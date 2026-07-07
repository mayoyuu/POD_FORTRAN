!> @file test_emdac_gap_scatter.f90
!> @brief EMDAC 长间隔散点测试：捕获第一个长时间观测间隔前后的粒子散布
!>
!> 处理 L1Halo-2 案例 (n=1,3,5)，在第一个长时间间隔（obs 10→11，
!> 约 10.5 天）前后采样 100k 粒子，保存为 .txt 散点文件和 .json GMM 文件。
program test_emdac_gap_scatter
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_filter_emdac_module, only: emdac_filter
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_obs_io_module, only: obs_record, preload_observations, &
                                  station_record, preload_stations, find_station_by_id
    use pod_measurement_base_module, only: observation_station
    use pod_data_format_module, only: load_initial_opm, write_json_opm
    use pod_basicmath_module, only: PI

    implicit none

    character(len=MAX_STRING_LEN) :: config_file, obs_file, site_json_file
    character(len=MAX_STRING_LEN) :: initial_json_file, output_prefix, arg_str
    integer :: i, num_args

    type(emdac_filter) :: my_filter
    type(observation_station) :: current_station
    type(obs_record), allocatable :: obs_list(:)
    type(station_record), allocatable :: station_list(:)

    real(DP) :: initial_mean(6), initial_cov(6,6), noise_Q(6,6), noise_R(2,2)
    real(DP) :: et_current, et_obs, y_meas(2), dt
    real(DP) :: final_mean(6), final_cov(6,6)
    type(uq_gmm_state_type) :: post_gap_gmm

    real(DP), allocatable :: particles(:, :)
    real(DP), allocatable :: prop_particles(:, :)
    integer :: obs_count, n_particles_actual

    ! 算法参数默认值
    integer :: opt_particles = 100000
    integer :: opt_da_order = 4
    integer :: opt_em_max_iter = 50
    real(DP) :: opt_em_tol = 1.0e-4_DP
    integer :: n_components = 5
    real(DP), parameter :: sigma_a = 1.0e-11_DP

    ! 默认路径
    config_file       = 'config/config.txt'
    obs_file          = 'OBS/L2Halo-2/L2Halo-2_supp_single_R91_1h.obs'
    site_json_file    = 'config/site.json'
    initial_json_file = 'OPM/L1Halo-2/L1Halo-2_init.opm.json'
    output_prefix     = 'output/L1Halo-2_emdac_gap'

    ! 解析命令行参数
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))
            case ('-config')
                call get_command_argument(i+1, arg_str); config_file = trim(arg_str); i = i + 1
            case ('-obs')
                call get_command_argument(i+1, arg_str); obs_file = trim(arg_str); i = i + 1
            case ('-init')
                call get_command_argument(i+1, arg_str); initial_json_file = trim(arg_str); i = i + 1
            case ('-site')
                call get_command_argument(i+1, arg_str); site_json_file = trim(arg_str); i = i + 1
            case ('-out')
                call get_command_argument(i+1, arg_str); output_prefix = trim(arg_str); i = i + 1
            case ('-p')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_particles; i = i + 1
            case ('-o')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_da_order; i = i + 1
            case ('-n')
                call get_command_argument(i+1, arg_str); read(arg_str, *) n_components; i = i + 1
            case ('-iter')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_em_max_iter; i = i + 1
            case ('-tol')
                call get_command_argument(i+1, arg_str); read(arg_str, *) opt_em_tol; i = i + 1
        end select
        i = i + 1
    end do

    if (len_trim(initial_json_file) == 0) then
        write(*,*) '[ERROR] -init <opm.json> is required'
        error stop
    end if

    write(*,*) '=================================================='
    write(*,*) '  EMDAC 长间隔散点测试'
    write(*,*) '=================================================='
    write(*,*) '初始 OPM    : ', trim(initial_json_file)
    write(*,*) '观测文件    : ', trim(obs_file)
    write(*,*) '粒子总数    : ', opt_particles
    write(*,*) 'GMM 分量数  : ', n_components
    write(*,*) 'DA 阶数     : ', opt_da_order
    write(*,*) '输出前缀    : ', trim(output_prefix)
    write(*,*) '--------------------------------------------------'

    ! 1. 初始化引擎
    call pod_engine_init(trim(config_file))

    ! 2. 加载初始状态（含 GMM）
    call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)

    ! 3. 初始化滤波器
    call my_filter%init(et_current, initial_mean, initial_cov, n_components, opt_da_order, &
                        n_part=opt_particles, opt_em_max_iter=opt_em_max_iter, opt_em_tol=opt_em_tol)
    call my_filter%set_n_particles(opt_particles)

    ! 4. 加载观测与测站
    call preload_observations(obs_file, obs_list)
    call preload_stations(site_json_file, station_list)
    write(*,*) '观测总数: ', size(obs_list)

    ! 5. 测量噪声
    noise_R = 0.0_DP
    noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2
    noise_R(2,2) = noise_R(1,1)

    ! 6. 处理观测 1 到 11（跨越第一个长时间间隔）
    do obs_count = 1, min(11, size(obs_list))
        et_obs = obs_list(obs_count)%et
        y_meas(1) = obs_list(obs_count)%ra
        y_meas(2) = obs_list(obs_count)%dec
        current_station = find_station_by_id(obs_list(obs_count)%station_id, station_list)

        ! 时间更新
        dt = et_obs - et_current
        noise_Q = 0.0_DP
        do i = 1, 3
            noise_Q(i,i)     = (dt**4 / 4.0_DP) * sigma_a**2
            noise_Q(i+3,i+3) = dt**2 * sigma_a**2
        end do

        write(*,'(A,I0,A,F12.1,A)') '  Obs #', obs_count, '  dt=', dt, 's'
        call my_filter%time_update(et_obs, noise_Q)

        ! 测量更新
        call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)

        if (obs_count == 10) then
            ! 测量更新后：采样后验散点供 UT 传播
            allocate(particles(6, opt_particles))
            call my_filter%sample_particles_from_gmm(particles)
            call write_scatter_txt(trim(output_prefix) // '_before_gap.txt', particles, opt_particles)
            deallocate(particles)
            write(*,*) '  [后验散点] before_gap 已保存 (obs #10, 测量更新后)'
        end if

        if (obs_count == 11) then
            ! 测量更新后：获取传播粒子作为「间隔后」散点 + GMM
            call my_filter%get_propagated_particles(prop_particles)
            n_particles_actual = size(prop_particles, 2)
            call write_scatter_txt(trim(output_prefix) // '_after_gap.txt', &
                                   prop_particles, n_particles_actual)
            deallocate(prop_particles)
            write(*,*) '  [散点] after_gap 已保存 (obs #11, 时间更新后, ', n_particles_actual, ' 粒子)'

            ! 保存 GMM 到 OPM JSON
            call my_filter%get_current_gmm(post_gap_gmm)
            call my_filter%get_current_state(final_mean)
            call my_filter%get_current_cov(final_cov)
            call write_json_opm(trim(output_prefix) // '_gap_gmm', &
                               final_mean, final_cov, post_gap_gmm, 0.0_DP, "DRO", et_obs)
            write(*,*) '  [GMM]  gap_gmm.opm.json 已保存'
            exit
        end if
    end do

    write(*,*) '=================================================='
    write(*,*) '  完成！'
    write(*,*) '=================================================='

contains

    !> 将粒子散布写入 rand_list 格式的文本文件
    !> 格式：行号 <TAB> x,y,z,vx,vy,vz,
    subroutine write_scatter_txt(filename, samples, n_particles)
        character(len=*), intent(in) :: filename
        real(DP), intent(in) :: samples(:,:)  ! (6, n_particles)
        integer, intent(in) :: n_particles
        integer :: u, ios, p
        character(len=256) :: full_filename

        full_filename = trim(filename)
        open(newunit=u, file=full_filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[ERROR] 无法创建文件: ', trim(full_filename)
            error stop
        end if

        do p = 1, n_particles
            write(u, '(I0,A,6(ES25.16E3,A))') &
                p, char(9), &
                samples(1,p), ',', samples(2,p), ',', samples(3,p), ',', &
                samples(4,p), ',', samples(5,p), ',', samples(6,p), ','
        end do

        close(u)
        write(*,*) '  已写入 ', n_particles, ' 行到: ', trim(full_filename)
    end subroutine write_scatter_txt

end program test_emdac_gap_scatter
