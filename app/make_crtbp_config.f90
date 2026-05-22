program make_crtbp_config
    use pod_global, only: DP, MAX_STRING_LEN
    implicit none

    character(len=MAX_STRING_LEN) :: json_path, output_path, arg_tmp
    character(len=:), allocatable :: content
    logical :: ok

    real(DP) :: LU, TU, mu, t_end
    real(DP) :: nominal_state(6), cov(6,6)
    character(len=10) :: method
    integer  :: n_samples, da_order
    real(DP) :: rel_tol, abs_tol, dt_min, dt_max
    integer  :: max_steps_int
    character(len=MAX_STRING_LEN) :: output_prefix

    real(DP) :: pos_factor, vel_factor, time_factor
    real(DP) :: pp_factor, pv_factor, vv_factor
    logical :: state_is_dim, cov_is_dim, time_is_dim
    integer :: i, j, unit

    ! ---- CLI parsing ----
    if (command_argument_count() < 2) then
        print *, 'Usage: make_crtbp_config --in <config.json> [--out <output.json>]'
        stop
    end if

    call get_command_argument(1, json_path)
    if (trim(json_path) /= '--in') then
        print *, 'Error: expected --in <config.json>'
        stop
    end if
    call get_command_argument(2, json_path)

    output_path = 'crtbp_nondim.json'
    if (command_argument_count() >= 4) then
        call get_command_argument(3, arg_tmp)
        if (trim(arg_tmp) == '--out') then
            call get_command_argument(4, output_path)
        end if
    end if

    ! ---- Read JSON file ----
    call read_file_content(json_path, content, ok)
    if (.not. ok) then
        print *, 'Error: cannot read file ', trim(json_path)
        stop
    end if

    ! ---- Parse scale factors ----
    call json_get_real(content, 'LU', LU, ok)
    if (.not. ok) then
        print *, 'Error: missing "LU" (length unit, e.g. distance between primaries in km)'
        stop
    end if

    call json_get_real(content, 'TU', TU, ok)
    if (.not. ok) then
        print *, 'Error: missing "TU" (time unit, e.g. 1/mean_motion in seconds)'
        stop
    end if

    print '(A,ES15.6)', 'LU = ', LU
    print '(A,ES15.6)', 'TU = ', TU

    ! ---- Parse required fields ----
    call json_get_real(content, 'mu', mu, ok)
    if (.not. ok) stop 'Error: missing "mu"'

    call json_get_array(content, 'initial_state', nominal_state, 6, ok)
    if (.not. ok) stop 'Error: missing "initial_state"'

    call json_get_matrix(content, 'initial_covariance', cov, 6, 6, ok)
    if (.not. ok) stop 'Error: missing "initial_covariance"'

    call json_get_real(content, 'propagation_time', t_end, ok)
    if (.not. ok) stop 'Error: missing "propagation_time"'

    call json_get_string(content, 'method', method, ok)
    if (.not. ok) method = 'MC'

    call json_get_int(content, 'num_samples', n_samples, ok)
    if (.not. ok) n_samples = 10000

    call json_get_int(content, 'da_order', da_order, ok)
    if (.not. ok) da_order = 2

    ! ---- Integrator defaults ----
    rel_tol = 1.0e-12_DP
    abs_tol = 1.0e-12_DP
    dt_min  = 1.0e-10_DP
    dt_max  = 0.1_DP
    max_steps_int = 100000

    call json_get_string(content, 'prefix', output_prefix, ok)
    if (.not. ok) output_prefix = './output/crtbp_uprop'

    deallocate(content)

    ! ---- Detect and convert dimensional inputs ----
    state_is_dim = (maxval(abs(nominal_state(1:3))) > 100.0_DP) .or. &
                   (maxval(abs(nominal_state(4:6))) > 10.0_DP)
    cov_is_dim  = maxval(abs(cov)) > 1.0_DP
    time_is_dim = t_end > 100.0_DP

    pos_factor = 1.0_DP / LU
    vel_factor = TU / LU
    time_factor = 1.0_DP / TU

    pp_factor = 1.0_DP / (LU * LU)
    pv_factor = TU / (LU * LU)
    vv_factor = (TU * TU) / (LU * LU)

    if (state_is_dim) then
        print *, 'Converting initial_state: dimensional -> non-dimensional'
        nominal_state(1:3) = nominal_state(1:3) * pos_factor
        nominal_state(4:6) = nominal_state(4:6) * vel_factor
    else
        print *, 'initial_state appears non-dimensional, leaving as-is'
    end if

    if (cov_is_dim) then
        print *, 'Converting initial_covariance: dimensional -> non-dimensional'
        do i = 1, 3
            do j = 1, 3
                cov(i,j) = cov(i,j) * pp_factor
            end do
        end do
        do i = 1, 3
            do j = 4, 6
                cov(i,j) = cov(i,j) * pv_factor
                cov(j,i) = cov(j,i) * pv_factor
            end do
        end do
        do i = 4, 6
            do j = 4, 6
                cov(i,j) = cov(i,j) * vv_factor
            end do
        end do
    else
        print *, 'initial_covariance appears non-dimensional, leaving as-is'
    end if

    if (time_is_dim) then
        print *, 'Converting propagation_time: dimensional -> non-dimensional'
        t_end = t_end * time_factor
    else
        print *, 'propagation_time appears non-dimensional, leaving as-is'
    end if

    ! ---- Validate non-dimensional input ----
    if (maxval(abs(nominal_state(1:3))) > 100.0_DP .or. &
        maxval(abs(nominal_state(4:6))) > 10.0_DP) then
        print *, ''
        print *, 'WARNING: Converted state still appears dimensional. Check LU/TU values.'
        print *, '  |pos| max = ', maxval(abs(nominal_state(1:3)))
        print *, '  |vel| max = ', maxval(abs(nominal_state(4:6)))
    end if

    ! ---- Write output JSON ----
    open(newunit=unit, file=trim(output_path), status='replace', action='write')
    write(unit, '(A)') '{'
    write(unit, '(A,F12.8,A)') '    "mu": ', mu, ','
    write(unit, '(A,*(F18.12,:,","))') &
        '    "initial_state": [', nominal_state(1), nominal_state(2), nominal_state(3), &
        nominal_state(4), nominal_state(5), nominal_state(6), '],'
    write(unit, '(A)') '    "initial_covariance": ['
    do i = 1, 6
        if (i < 6) then
            write(unit, '(A,*(ES22.14,:,","),A)') &
                '        [', cov(i,1), cov(i,2), cov(i,3), &
                cov(i,4), cov(i,5), cov(i,6), '],'
        else
            write(unit, '(A,*(ES22.14,:,","),A)') &
                '        [', cov(i,1), cov(i,2), cov(i,3), &
                cov(i,4), cov(i,5), cov(i,6), ']'
        end if
    end do
    write(unit, '(A)') '    ],'
    write(unit, '(A,F14.6,A)') '    "propagation_time": ', t_end, ','
    write(unit, '(A,A,A)') '    "method": "', trim(method), '",'
    write(unit, '(A,I0,A)') '    "num_samples": ', n_samples, ','
    write(unit, '(A,I0,A)') '    "da_order": ', da_order, ','
    write(unit, '(A)') '    "integrator": {'
    write(unit, '(A,ES12.4,A)') '        "rel_tol": ', rel_tol, ','
    write(unit, '(A,ES12.4,A)') '        "abs_tol": ', abs_tol, ','
    write(unit, '(A,ES12.4,A)') '        "min_step": ', dt_min, ','
    write(unit, '(A,F10.6,A)') '        "max_step": ', dt_max, ','
    write(unit, '(A,I0)') '        "max_steps": ', max_steps_int
    write(unit, '(A)') '    },'
    write(unit, '(A)') '    "output": {'
    write(unit, '(A,A,A)') '        "prefix": "', trim(output_prefix), '"'
    write(unit, '(A)') '    }'
    write(unit, '(A)') '}'
    close(unit)

    print *, ''
    print *, 'Non-dimensional config written to: ', trim(output_path)

contains

    ! ================================================================
    !  JSON PARSER (minimal, embedded)
    ! ================================================================

    subroutine read_file_content(filename, content, ok)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable, intent(out) :: content
        logical, intent(out) :: ok
        character(len=4096) :: line
        integer :: unit, ios
        character(len=:), allocatable :: temp

        ok = .false.
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        allocate(character(len=0) :: content)
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            temp = content // trim(adjustl(line))
            deallocate(content)
            content = temp
        end do
        close(unit)
        if (len(content) > 0) ok = .true.
    end subroutine read_file_content

    function find_key_pos(content, key) result(pos)
        character(len=*), intent(in) :: content, key
        integer :: pos
        character(len=:), allocatable :: search
        search = '"' // trim(key) // '"'
        pos = index(content, search)
    end function find_key_pos

    subroutine json_get_real(content, key, value, found)
        character(len=*), intent(in) :: content, key
        real(DP), intent(out) :: value
        logical, intent(out) :: found
        integer :: p, p_end, ios
        character(len=256) :: token

        found = .false.
        value = 0.0_DP
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, ':')
        if (p == 0) return
        p = p + 1
        p = skip_whitespace(content, p)
        p_end = scan_forward(content, p, ',}]')
        if (p_end == 0) p_end = len(content)
        token = content(p:p_end-1)
        token = trim(adjustl(token))
        read(token, *, iostat=ios) value
        if (ios == 0) found = .true.
    end subroutine json_get_real

    subroutine json_get_int(content, key, value, found)
        character(len=*), intent(in) :: content, key
        integer, intent(out) :: value
        logical, intent(out) :: found
        integer :: p, p_end, ios
        character(len=256) :: token

        found = .false.
        value = 0
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, ':')
        if (p == 0) return
        p = p + 1
        p = skip_whitespace(content, p)
        p_end = scan_forward(content, p, ',}]')
        if (p_end == 0) p_end = len(content)
        token = content(p:p_end-1)
        token = trim(adjustl(token))
        read(token, *, iostat=ios) value
        if (ios == 0) found = .true.
    end subroutine json_get_int

    subroutine json_get_string(content, key, value, found)
        character(len=*), intent(in) :: content, key
        character(len=*), intent(out) :: value
        logical, intent(out) :: found
        integer :: p, p_start, p_end

        found = .false.
        value = ''
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, ':')
        if (p == 0) return
        p = p + 1
        p = skip_whitespace(content, p)
        if (p > len(content)) return
        if (content(p:p) /= '"') return
        p_start = p + 1
        p_end = scan_forward(content, p_start, '"')
        if (p_end == 0) return
        value = content(p_start:p_end-1)
        found = .true.
    end subroutine json_get_string

    subroutine json_get_array(content, key, arr, n, found)
        character(len=*), intent(in) :: content, key
        real(DP), intent(out) :: arr(n)
        integer, intent(in) :: n
        logical, intent(out) :: found
        integer :: p, p_end, i, ios
        character(len=4096) :: segment

        found = .false.
        arr = 0.0_DP
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, '[')
        if (p == 0) return
        p = p + 1
        p_end = find_matching_bracket(content, p, '[', ']')
        if (p_end == 0) return
        segment = content(p:p_end-1)
        do i = 1, n
            p = skip_whitespace(segment, 1)
            p_end = scan_forward(segment, p, ',]')
            if (p_end == 0) p_end = len_trim(segment) + 1
            read(segment(p:p_end-1), *, iostat=ios) arr(i)
            if (ios /= 0) return
            if (i < n .and. p_end <= len_trim(segment)) &
                segment = segment(p_end+1:)
        end do
        found = .true.
    end subroutine json_get_array

    subroutine json_get_matrix(content, key, mat, rows, cols, found)
        character(len=*), intent(in) :: content, key
        real(DP), intent(out) :: mat(rows, cols)
        integer, intent(in) :: rows, cols
        logical, intent(out) :: found
        integer :: p, p_end, i, j, ios, p_seg
        character(len=4096) :: outer_seg, row_seg

        found = .false.
        mat = 0.0_DP
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, '[')
        if (p == 0) return
        p = p + 1
        p_end = find_matching_bracket(content, p, '[', ']')
        if (p_end == 0) return
        outer_seg = content(p:p_end-1)

        do i = 1, rows
            p = scan_forward(outer_seg, 1, '[')
            if (p == 0) return
            p = p + 1
            p_end = find_matching_bracket(outer_seg, p, '[', ']')
            if (p_end == 0) return
            row_seg = outer_seg(p:p_end-1)
            do j = 1, cols
                p_seg = skip_whitespace(row_seg, 1)
                p_end = scan_forward(row_seg, p_seg, ',]')
                if (p_end == 0) p_end = len_trim(row_seg) + 1
                read(row_seg(p_seg:p_end-1), *, iostat=ios) mat(i, j)
                if (ios /= 0) return
                if (j < cols .and. p_end <= len_trim(row_seg)) &
                    row_seg = row_seg(p_end+1:)
            end do
            p_end = scan_forward(outer_seg, p, ']')
            if (p_end == 0) return
            outer_seg = outer_seg(p_end+1:)
        end do
        found = .true.
    end subroutine json_get_matrix

    ! ---- Utility functions ----

    function scan_forward(str, start, ch) result(pos)
        character(len=*), intent(in) :: str
        integer, intent(in) :: start
        character, intent(in) :: ch
        integer :: pos
        pos = index(str(start:), ch)
        if (pos > 0) pos = pos + start - 1
    end function scan_forward

    function skip_whitespace(str, start) result(pos)
        character(len=*), intent(in) :: str
        integer, intent(in) :: start
        integer :: pos
        pos = start
        do while (pos <= len(str))
            if (str(pos:pos) /= ' ' .and. str(pos:pos) /= achar(9) .and. &
                str(pos:pos) /= achar(10) .and. str(pos:pos) /= achar(13)) exit
            pos = pos + 1
        end do
    end function skip_whitespace

    function find_matching_bracket(str, open_pos, open_ch, close_ch) result(close_pos)
        character(len=*), intent(in) :: str
        integer, intent(in) :: open_pos
        character, intent(in) :: open_ch, close_ch
        integer :: close_pos
        integer :: depth, i

        depth = 1
        close_pos = 0
        do i = open_pos + 1, len(str)
            if (str(i:i) == open_ch) depth = depth + 1
            if (str(i:i) == close_ch) depth = depth - 1
            if (depth == 0) then
                close_pos = i
                return
            end if
        end do
    end function find_matching_bracket

end program make_crtbp_config
