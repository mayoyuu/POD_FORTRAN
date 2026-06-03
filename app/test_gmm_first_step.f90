!> @file test_gmm_first_step.f90
!> @brief 保存第一次测量更新前后的 GMM 快照，用于调试 EMDAC 测量更新过程
program test_gmm_first_step
    use pod_global, only: DP, MAX_STRING_LEN
    use pod_engine_module, only: pod_engine_init
    use pod_filter_emdac_module, only: emdac_filter
    use pod_uq_gmm_state_module, only: uq_gmm_state_type
    use pod_obs_io_module, only: obs_record, preload_observations, &
                                  station_record, preload_stations, find_station_by_id
    use pod_measurement_base_module, only: observation_station
    use pod_data_format_module, only: load_initial_opm
    use pod_basicmath_module, only: PI
    use pod_spice, only: et2utc

    implicit none

    character(len=MAX_STRING_LEN) :: config_file, obs_file, site_json_file
    character(len=MAX_STRING_LEN) :: initial_json_file, output_gmm_file, arg_str
    integer :: i, num_args

    type(emdac_filter) :: my_filter
    type(observation_station) :: current_station
    type(obs_record), allocatable :: obs_list(:)
    type(station_record), allocatable :: station_list(:)

    real(DP) :: initial_mean(6), initial_cov(6,6), noise_Q(6,6), noise_R(2,2)
    real(DP) :: et_current, et_obs, y_meas(2), dt
    real(DP), parameter :: sigma_a = 1.0e-11_DP
    real(DP) :: global_means(2, 6), global_covs(2, 6, 6)

    type(uq_gmm_state_type) :: gmm_snapshots(2)
    real(DP) :: gmm_epochs(2)
    character(len=64) :: label(2)

    ! 算法参数默认值
    integer :: opt_particles = 10000
    integer :: opt_da_order = 4
    integer :: opt_em_max_iter = 50
    real(DP) :: opt_em_tol = 1.0e-4_DP
    integer :: n_components = 5

    ! 默认路径
    config_file       = 'config/dummy_test_config.txt'
    obs_file          = 'OBS/L1Halo-1/L1Halo-1_supp_single_R91_1h.obs'
    site_json_file    = 'config/site.json'
    initial_json_file = 'OPM/L1Halo-1/L1Halo-1_init.opm.json'
    output_gmm_file   = 'OPM/L1Halo-1/L1Halo-1_supp_single_R91_1h_emdac_n5'

    ! 解析命令行参数
    num_args = command_argument_count()
    i = 1
    do while (i <= num_args)
        call get_command_argument(i, arg_str)
        select case (trim(arg_str))
            case ('-obs')
                call get_command_argument(i+1, arg_str); obs_file = trim(arg_str); i = i + 1
            case ('-init')
                call get_command_argument(i+1, arg_str); initial_json_file = trim(arg_str); i = i + 1
            case ('-out')
                call get_command_argument(i+1, arg_str); output_gmm_file = trim(arg_str); i = i + 1
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

    write(*,*) '=================================================='
    write(*,*) '  EMDAC 首步 GMM 快照调试工具'
    write(*,*) '=================================================='
    write(*,*) '粒子总数     : ', opt_particles
    write(*,*) 'GMM 分量数   : ', n_components
    write(*,*) 'DA 阶数      : ', opt_da_order
    write(*,*) '输出文件     : ', trim(output_gmm_file)
    write(*,*) '--------------------------------------------------'

    ! 1. 初始化引擎
    call pod_engine_init(trim(config_file))

    ! 2. 加载初始状态
    call load_initial_opm(initial_json_file, et_current, initial_mean, initial_cov)

    ! 3. 初始化滤波器
    call my_filter%init(et_current, initial_mean, initial_cov, n_components, opt_da_order, &
                        n_part=opt_particles, opt_em_max_iter=opt_em_max_iter, opt_em_tol=opt_em_tol)

    ! 4. 加载观测与测站
    call preload_observations(obs_file, obs_list)
    call preload_stations(site_json_file, station_list)

    ! 5. 测量噪声
    noise_R = 0.0_DP
    noise_R(1,1) = (0.1_DP * PI / 180.0_DP / 3600.0_DP)**2
    noise_R(2,2) = noise_R(1,1)

    ! 6. 只处理第一条观测
    et_obs = obs_list(1)%et
    y_meas(1) = obs_list(1)%ra
    y_meas(2) = obs_list(1)%dec
    current_station = find_station_by_id(obs_list(1)%station_id, station_list)

    write(*,*) '  观测历元: ', et_obs
    write(*,*) '  观测值  : ', y_meas
    write(*,*) '  测站    : ', trim(current_station%name)

    ! 7. 时间更新
    dt = et_obs - et_current
    noise_Q = 0.0_DP
    do i = 1, 3
        noise_Q(i,i)     = (dt**4 / 4.0_DP) * sigma_a**2
        noise_Q(i+3,i+3) = dt**2 * sigma_a**2
    end do

    write(*,*) '  时间步长 dt = ', dt, ' s'
    write(*,*) '  执行时间更新...'
    call my_filter%time_update(et_obs, noise_Q)

    ! 8. 保存测量更新前的 GMM
    call my_filter%get_current_gmm(gmm_snapshots(1))
    gmm_epochs(1) = et_obs
    label(1) = 'AFTER_TIME_UPDATE'

    write(*,*) '  [快照1] 时间更新后 GMM 已保存'
    write(*,*) '    分量数: ', gmm_snapshots(1)%n_components
    do i = 1, gmm_snapshots(1)%n_components
        write(*,*) '    分量 ', i, ' 权重: ', gmm_snapshots(1)%components(i)%weight
    end do

    call my_filter%get_current_state(global_means(1,:))
    call my_filter%get_current_cov(global_covs(1,:,:))

    ! 9. 执行测量更新
    write(*,*) '  执行测量更新...'
    call my_filter%measurement_update(y_meas, noise_R, et_obs, current_station)

    ! 10. 保存测量更新后的 GMM
    call my_filter%get_current_gmm(gmm_snapshots(2))
    gmm_epochs(2) = et_obs
    label(2) = 'AFTER_MEASUREMENT_UPDATE'

    write(*,*) '  [快照2] 测量更新后 GMM 已保存'
    write(*,*) '    分量数: ', gmm_snapshots(2)%n_components
    do i = 1, gmm_snapshots(2)%n_components
        write(*,*) '    分量 ', i, ' 权重: ', gmm_snapshots(2)%components(i)%weight, &
                   ' 均值前3维: ', gmm_snapshots(2)%components(i)%mean(1:3)
    end do

    call my_filter%get_current_state(global_means(2,:))
    call my_filter%get_current_cov(global_covs(2,:,:))

    ! 11. 写入文件
    write(*,*) '  写入 GMM 快照...'
    call write_gmm_snapshots_with_labels(output_gmm_file, gmm_epochs, gmm_snapshots, &
                                         global_means, global_covs, label, 2)

    write(*,*) '=================================================='
    write(*,*) '  完成！输出: ', trim(output_gmm_file) // '.gmms.json'
    write(*,*) '=================================================='

contains

    !> 带标签的 GMM 快照写入（与 data_format 中类似，但多了 label 字段）
    subroutine write_gmm_snapshots_with_labels(filename, epochs, gmms, &
                                               global_means, global_covs, labels, n_steps)
        character(len=*), intent(in) :: filename
        real(DP), intent(in)         :: epochs(:)
        type(uq_gmm_state_type), intent(in) :: gmms(:)
        real(DP), intent(in)         :: global_means(:,:), global_covs(:,:,:)
        character(len=*), intent(in) :: labels(:)
        integer, intent(in)          :: n_steps

        integer :: u, ios, s, i, j
        character(len=64) :: epoch_str
        character(len=256) :: full_filename

        full_filename = trim(filename) // ".gmms.json"

        open(newunit=u, file=full_filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[错误] 无法创建文件: ', trim(full_filename)
            return
        end if

        write(u, '(A)') '{'
        write(u, '(A)') '    "STEPS": ['

        do s = 1, n_steps
            call et2utc(epochs(s), 'ISOC', 3, epoch_str)

            write(u, '(A)') '        {'
            write(u, '(A,I0,A)')          '            "INDEX": ', s, ','
            write(u, '(A,A,A)')           '            "LABEL": "', trim(labels(s)), '",'
            write(u, '(A,A,A)')           '            "EPOCH": "', trim(epoch_str), '",'
            write(u, '(A,F16.4,A)')       '            "ET": ', epochs(s), ','
            write(u, '(A)')               '            "GLOBAL": {'

            write(u, '(A,5(ES22.15,", "),ES22.15,A)') &
                '                "MEAN": [', global_means(s, 1:6), '],'
            write(u, '(A)')               '                "COV": ['
            do j = 1, 6
                if (j < 6) then
                    write(u, '(A,5(ES22.15,", "),ES22.15,A)') &
                        '                    [', global_covs(s, j, 1:6), '],'
                else
                    write(u, '(A,5(ES22.15,", "),ES22.15,A)') &
                        '                    [', global_covs(s, j, 1:6), ']'
                end if
            end do
            write(u, '(A)')               '                ]'
            write(u, '(A)')               '            },'

            write(u, '(A,I0,A)')          '            "N_COMPONENTS": ', gmms(s)%n_components, ','
            write(u, '(A)')               '            "COMPONENTS": ['

            do i = 1, gmms(s)%n_components
                write(u, '(A)') '                {'
                write(u, '(A,I0,A)')      '                    "INDEX": ', i, ','
                write(u, '(A,ES22.15,A)') '                    "WEIGHT": ', gmms(s)%components(i)%weight, ','
                write(u, '(A,5(ES22.15,", "),ES22.15,A)') &
                    '                    "MEAN": [', gmms(s)%components(i)%mean, '],'
                write(u, '(A)')           '                    "COV": ['
                do j = 1, 6
                    if (j < 6) then
                        write(u, '(A,5(ES22.15,", "),ES22.15,A)') &
                            '                        [', gmms(s)%components(i)%cov(j, 1:6), '],'
                    else
                        write(u, '(A,5(ES22.15,", "),ES22.15,A)') &
                            '                        [', gmms(s)%components(i)%cov(j, 1:6), ']'
                    end if
                end do
                write(u, '(A)') '                    ]'

                if (i < gmms(s)%n_components) then
                    write(u, '(A)') '                },'
                else
                    write(u, '(A)') '                }'
                end if
            end do

            write(u, '(A)') '            ]'

            if (s < n_steps) then
                write(u, '(A)') '        },'
            else
                write(u, '(A)') '        }'
            end if
        end do

        write(u, '(A)') '    ]'
        write(u, '(A)') '}'

        close(u)
    end subroutine write_gmm_snapshots_with_labels

end program test_gmm_first_step
