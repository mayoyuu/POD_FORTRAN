program run_CRTBP_uprop
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_uq_crtbp_mc_module, only: crtbp_mc_propagate
    use pod_uq_crtbp_da_module, only: crtbp_da_propagate
    implicit none

    character(len=MAX_STRING_LEN) :: json_path
    character(len=:), allocatable :: content
    character(len=:), allocatable :: obj_content
    logical :: ok

    real(DP) :: mu, t_end
    real(DP) :: nominal_state(6)
    real(DP) :: cov(6,6)
    character(len=10) :: method
    integer  :: n_samples, da_order
    real(DP) :: rel_tol, abs_tol, dt_min, dt_max
    integer  :: max_steps_int
    character(len=MAX_STRING_LEN) :: output_prefix

    real(DP), allocatable :: final_samples(:,:)
    real(DP) :: final_mean(6), final_cov(6,6), propagated_ref(6)

    integer :: i, j, file_unit

    ! ---- CLI parsing ----
    if (command_argument_count() < 2) then
        print *, 'Usage: run_CRTBP_uprop --in <config.json>'
        stop
    end if

    call get_command_argument(1, json_path)
    if (trim(json_path) /= '--in') then
        print *, 'Error: expected --in <config.json>'
        stop
    end if
    call get_command_argument(2, json_path)

    ! ---- Read JSON file ----
    call read_file_content(json_path, content, ok)
    if (.not. ok) then
        print *, 'Error: cannot read file ', trim(json_path)
        stop
    end if

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
    if (.not. ok) stop 'Error: missing "method"'

    call json_get_int(content, 'num_samples', n_samples, ok)
    if (.not. ok) n_samples = 10000

    call json_get_int(content, 'da_order', da_order, ok)
    if (.not. ok) da_order = 2

    ! ---- Parse nested integrator object ----
    call json_get_object(content, 'integrator', obj_content, ok)
    if (ok) then
        call json_get_real(obj_content, 'rel_tol', rel_tol, ok)
        if (.not. ok) rel_tol = 1.0e-12_DP
        call json_get_real(obj_content, 'abs_tol', abs_tol, ok)
        if (.not. ok) abs_tol = 1.0e-12_DP
        call json_get_real(obj_content, 'min_step', dt_min, ok)
        if (.not. ok) dt_min = 1.0e-10_DP
        call json_get_real(obj_content, 'max_step', dt_max, ok)
        if (.not. ok) dt_max = 0.1_DP
        call json_get_int(obj_content, 'max_steps', max_steps_int, ok)
        if (.not. ok) max_steps_int = 100000
        deallocate(obj_content)
    else
        rel_tol = 1.0e-12_DP
        abs_tol = 1.0e-12_DP
        dt_min  = 1.0e-10_DP
        dt_max  = 0.1_DP
        max_steps_int = 100000
    end if

    ! ---- Parse output prefix ----
    call json_get_object(content, 'output', obj_content, ok)
    if (ok) then
        call json_get_string(obj_content, 'prefix', output_prefix, ok)
        if (.not. ok) output_prefix = './output/crtbp_uprop'
        deallocate(obj_content)
    else
        output_prefix = './output/crtbp_uprop'
    end if

    deallocate(content)

    ! ---- Validate non-dimensional input ----
    if (maxval(abs(nominal_state(1:3))) > 100.0_DP .or. &
        maxval(abs(nominal_state(4:6))) > 10.0_DP) then
        print *, 'ERROR: Input state appears to be dimensional.'
        print *, 'CRTBP expects non-dimensional units (positions ~O(1), velocities ~O(1)).'
        print *, 'Position magnitudes found: ', maxval(abs(nominal_state(1:3)))
        print *, 'Velocity magnitudes found: ', maxval(abs(nominal_state(4:6)))
        print *, 'If input is truly non-dimensional, adjust the threshold in source code.'
        stop
    end if

    ! ---- Print config summary ----
    print *, '========================================'
    print *, '  CRTBP Uncertainty Propagation'
    print *, '========================================'
    print '(A,F8.6)',  '  mu              = ', mu
    print '(A,6F10.6)', '  initial state   = ', nominal_state
    print '(A,F12.4)',  '  propagation T   = ', t_end
    print '(A,A)',      '  method          = ', trim(method)
    if (trim(method) == 'MC') then
        print '(A,I0)', '  num_samples     = ', n_samples
    else
        print '(A,I0)', '  da_order        = ', da_order
        print '(A,I0)', '  num_samples     = ', n_samples
    end if
    print '(A,ES10.3)', '  rel_tol         = ', rel_tol
    print '(A,ES10.3)', '  abs_tol         = ', abs_tol
    print '(A,A)',      '  output prefix   = ', trim(output_prefix)
    print *, '========================================'

    ! ---- Route to MC or DA ----
    if (trim(method) == 'MC') then
        call crtbp_mc_propagate(nominal_state, cov, mu, t_end, &
            n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
            final_samples, final_mean, final_cov, .true.)
    else if (trim(method) == 'DA') then
        call crtbp_da_propagate(nominal_state, cov, mu, t_end, &
            da_order, n_samples, rel_tol, abs_tol, dt_min, dt_max, max_steps_int, &
            final_samples, final_mean, final_cov, propagated_ref, .true.)
    else
        print *, 'Error: method must be "MC" or "DA", got: ', trim(method)
        stop
    end if

    ! ---- Write output CSV files ----
    call write_csv_particles(trim(output_prefix) // '_particles.csv', final_samples)
    call write_csv_stats(trim(output_prefix) // '_stats.csv', final_mean, final_cov)

    print *, 'Results saved to:'
    print *, '  ', trim(output_prefix) // '_particles.csv'
    print *, '  ', trim(output_prefix) // '_stats.csv'

    deallocate(final_samples)

contains

    ! ================================================================
    !  JSON PARSER (embedded, minimal)
    ! ================================================================

    subroutine read_file_content(filename, content, ok)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable, intent(out) :: content
        logical, intent(out) :: ok
        character(len=4096) :: line
        integer :: unit, ios, total_len
        character(len=:), allocatable :: temp

        ok = .false.
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) return

        allocate(character(len=0) :: content)
        total_len = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            temp = content // trim(adjustl(line))
            deallocate(content)
            content = temp
        end do
        close(unit)
        total_len = len(content)
        if (total_len > 0) ok = .true.
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

    subroutine json_get_object(content, key, obj_content, found)
        character(len=*), intent(in) :: content, key
        character(len=:), allocatable, intent(out) :: obj_content
        logical, intent(out) :: found
        integer :: p, p_end

        found = .false.
        allocate(character(len=0) :: obj_content)
        p = find_key_pos(content, key)
        if (p == 0) return
        p = p + len('"' // trim(key) // '"')
        p = scan_forward(content, p, '{')
        if (p == 0) return
        p_end = find_matching_bracket(content, p, '{', '}')
        if (p_end == 0) return
        obj_content = content(p:p_end)
        found = .true.
    end subroutine json_get_object

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

    ! ---- Output writers ----

    subroutine write_csv_particles(filename, samples)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: samples(:,:)
        integer :: unit, i

        open(newunit=unit, file=trim(filename), status='replace', action='write')
        write(unit, '(A)') 'x,y,z,vx,vy,vz'
        do i = 1, size(samples, 2)
            write(unit, '(*(ES22.14, :, ","))') samples(:, i)
        end do
        close(unit)
    end subroutine write_csv_particles

    subroutine write_csv_stats(filename, mean_vec, cov_mat)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: mean_vec(6), cov_mat(6,6)
        integer :: unit, i

        open(newunit=unit, file=trim(filename), status='replace', action='write')
        write(unit, '(A)') '# Mean'
        write(unit, '(*(ES22.14, :, ","))') mean_vec(:)
        write(unit, '(A)') '# Covariance Matrix'
        do i = 1, 6
            write(unit, '(*(ES22.14, :, ","))') cov_mat(i, :)
        end do
        close(unit)
    end subroutine write_csv_stats

end program run_CRTBP_uprop
