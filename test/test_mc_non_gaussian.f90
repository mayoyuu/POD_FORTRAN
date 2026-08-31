! !> @file test_mc_non_gaussian.f90
! !> @brief MC 误差传播：检测初始高斯分布在 CRTBP 中传播后何时不再符合高斯
! !>
! !> 用法: fpm run test_mc_non_gaussian -- -init OPM/DRO-1/DRO-1_init.opm.json -days 30 -step 2
! !> 每个程序 run 一个轨道，通过脚本依次调用。
! program test_mc_non_gaussian
!     use pod_global, only: DP, MAX_STRING_LEN
!     use pod_engine_module, only: pod_engine_init
!     use pod_orbit_propagation, only: propagate_orbit, orbit_state, propagation_result
!     use pod_integrator_module, only: METHOD_RKF78
!     use pod_data_format_module, only: load_initial_opm
!     use pod_random_module, only: init_random_seed, generate_multivariate_normal
!     use pod_basicmath_module, only: dpotrf
!     implicit none

!     character(len=MAX_STRING_LEN) :: init_file, config_file, output_file, arg_str
!     real(DP) :: et0, init_state(6), init_cov(6, 6)
!     real(DP), allocatable :: particles(:, :)
!     integer :: i, n_particles, max_days, step_days, day, n_steps
!     integer :: u_out, ios, arg_i
!     real(DP) :: current_epoch, dt_step
!     type(orbit_state) :: orb0
!     type(propagation_result) :: prop_result
!     logical :: is_gaussian
!     integer :: first_non_gaussian_day

!     ! 默认参数
!     config_file   = 'config/config.txt'
!     init_file     = ''
!     max_days      = 30
!     step_days     = 2
!     n_particles   = 10000
!     output_file   = ''

!     ! 命令行参数
!     arg_i = 1
!     do while (arg_i <= command_argument_count())
!         call get_command_argument(arg_i, arg_str)
!         select case (trim(arg_str))
!             case ('-config')
!                 call get_command_argument(arg_i+1, arg_str); config_file = trim(arg_str); arg_i = arg_i + 1
!             case ('-init')
!                 call get_command_argument(arg_i+1, arg_str); init_file = trim(arg_str); arg_i = arg_i + 1
!             case ('-days')
!                 call get_command_argument(arg_i+1, arg_str); read(arg_str, *) max_days; arg_i = arg_i + 1
!             case ('-step')
!                 call get_command_argument(arg_i+1, arg_str); read(arg_str, *) step_days; arg_i = arg_i + 1
!             case ('-n')
!                 call get_command_argument(arg_i+1, arg_str); read(arg_str, *) n_particles; arg_i = arg_i + 1
!             case ('-out')
!                 call get_command_argument(arg_i+1, arg_str); output_file = trim(arg_str); arg_i = arg_i + 1
!         end select
!         arg_i = arg_i + 1
!     end do

!     if (len_trim(init_file) == 0) then
!         write(*,*) '用法: test_mc_non_gaussian -init <OPM.json> [-days 30] [-step 2] [-n 10000] [-out result.txt]'
!         error stop
!     end if

!     ! 自动生成输出文件名：与 init 文件同目录
!     if (len_trim(output_file) == 0) then
!         output_file = trim(init_file(1:len_trim(init_file)-5)) // '_mc_non_gaussian.txt'
!     end if

!     write(*,*) '=================================================='
!     write(*,*) '  MC 非高斯性检测'
!     write(*,*) '=================================================='
!     write(*,*) 'OPM 文件  : ', trim(init_file)
!     write(*,*) '总天数    : ', max_days
!     write(*,*) '检测步长  : ', step_days, ' 天'
!     write(*,*) '粒子数    : ', n_particles
!     write(*,*) '输出文件  : ', trim(output_file)

!     ! 1. 初始化引擎
!     call pod_engine_init(trim(config_file))

!     ! 2. 加载初始状态
!     call load_initial_opm(trim(init_file), et0, init_state, init_cov)
!     write(*,*) '初始历元  : ', et0
!     write(*,*) '初始协方差迹: ', init_cov(1,1)+init_cov(2,2)+init_cov(3,3)+ &
!                                         init_cov(4,4)+init_cov(5,5)+init_cov(6,6)

!     ! 3. 生成粒子
!     allocate(particles(6, n_particles))
!     call init_random_seed(.true.)
!     call generate_multivariate_normal(init_state, init_cov, particles)
!     write(*,*) '已生成 ', n_particles, ' 个粒子'

!     ! 4. 打开输出文件
!     open(newunit=u_out, file=trim(output_file), status='replace', action='write', iostat=ios)
!     if (ios /= 0) then
!         write(*,*) '[ERROR] 无法创建: ', trim(output_file)
!         error stop
!     end if
!     write(u_out, '(A)') '# day  skew_x  skew_y  skew_z  skew_vx  skew_vy  skew_vz' &
!         // '  ekurt_x  ekurt_y  ekurt_z  ekurt_vx  ekurt_vy  ekurt_vz' &
!         // '  md_mean  md_std  md_min  md_max  is_gaussian'
!     close(u_out)

!     ! 5. Day 0 高斯性检查
!     write(*,*) ''
!     write(*,*) '--- Day 0 ---'
!     call check_and_log(0, particles, n_particles, output_file, is_gaussian)

!     ! 6. 逐段传播与检测
!     n_steps = max_days / step_days
!     dt_step = real(step_days, DP) * 86400.0_DP
!     current_epoch = et0
!     first_non_gaussian_day = -1

!     do day = 1, n_steps
!         write(*,*) ''
!         write(*,*) '--- Propagating to Day ', day * step_days, ' ---'

!         !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
!         do i = 1, n_particles
!             orb0%epoch = current_epoch
!             orb0%state = particles(:, i)
!             call propagate_orbit(orb0, dt_step, METHOD_RKF78, prop_result)
!             particles(:, i) = prop_result%states(prop_result%n_steps, :)
!         end do
!         !$OMP END PARALLEL DO

!         current_epoch = current_epoch + dt_step

!         call check_and_log(day * step_days, particles, n_particles, output_file, is_gaussian)

!         if (.not. is_gaussian .and. first_non_gaussian_day < 0) then
!             first_non_gaussian_day = day * step_days
!         end if
!     end do

!     deallocate(particles)

!     write(*,*) ''
!     write(*,*) '=================================================='
!     if (first_non_gaussian_day > 0) then
!         write(*,*) '  首次检测到非高斯: Day ', first_non_gaussian_day
!     else
!         write(*,*) '  在 ', max_days, ' 天内未检测到明显非高斯'
!     end if
!     write(*,*) '  结果已保存至: ', trim(output_file)
!     write(*,*) '=================================================='

! contains

!     subroutine check_and_log(day_num, samples, n, out_file, is_gauss)
!         integer, intent(in) :: day_num, n
!         real(DP), intent(in) :: samples(:, :)
!         character(len=*), intent(in) :: out_file
!         logical, intent(out) :: is_gauss

!         real(DP) :: mean_vec(6), cov_mat(6, 6)
!         real(DP) :: skew(6), ekurt(6), md_mean, md_std, md_min, md_max
!         integer :: u_out, ios

!         call compute_mean_cov(samples, n, mean_vec, cov_mat)
!         call compute_gaussianity(samples, n, mean_vec, cov_mat, skew, ekurt, &
!                                   md_mean, md_std, md_min, md_max)

!         ! 判定标准（md_mean 是数学恒等式，始终≈6，不参与判定）：
!         !   马氏距离 std 显著偏离 sqrt(12)≈3.46，或 max 远超 50
!         !   任意维度偏度或峰度显著偏离 0
!         is_gauss = (md_std < 7.0_DP) .and. (md_max < 50.0_DP) .and. &
!                    all(abs(skew) < 0.5_DP) .and. all(abs(ekurt) < 1.0_DP)

!         write(*, '(A,I3,A,L1)') '  Day ', day_num, '  高斯: ', is_gauss
!         write(*, '(A,6F8.3)')   '    偏度      : ', skew
!         write(*, '(A,6F8.3)')   '    超额峰度  : ', ekurt
!         write(*, '(A,F8.3,A,F8.3,A,F8.3,A,F8.3)') &
!             '    马氏距离  : mean=', md_mean, ' std=', md_std, ' min=', md_min, ' max=', md_max, &
!             '  (mean≈6恒成立，不参与判定)'

!         open(newunit=u_out, file=trim(out_file), status='old', position='append', &
!              action='write', iostat=ios)
!         if (ios /= 0) return
!         write(u_out, '(I4,14F9.4,F9.3,F9.3,L2)') day_num, &
!             skew(1), skew(2), skew(3), skew(4), skew(5), skew(6), &
!             ekurt(1), ekurt(2), ekurt(3), ekurt(4), ekurt(5), ekurt(6), &
!             md_mean, md_std, md_min, md_max, is_gauss
!         close(u_out)
!     end subroutine check_and_log

!     subroutine compute_mean_cov(samples, n, mean_vec, cov_mat)
!         real(DP), intent(in) :: samples(:, :)
!         integer, intent(in) :: n
!         real(DP), intent(out) :: mean_vec(6), cov_mat(6, 6)
!         integer :: i, j

!         mean_vec = 0.0_DP
!         do i = 1, n
!             mean_vec = mean_vec + samples(:, i)
!         end do
!         mean_vec = mean_vec / real(n, DP)

!         cov_mat = 0.0_DP
!         do i = 1, n
!             do j = 1, 6
!                 cov_mat(:, j) = cov_mat(:, j) &
!                     + (samples(:, i) - mean_vec) * (samples(j, i) - mean_vec(j))
!             end do
!         end do
!         cov_mat = cov_mat / real(n - 1, DP)
!     end subroutine compute_mean_cov

!     subroutine compute_gaussianity(samples, n, mean_vec, cov_mat, skew, ekurt, &
!                                     md_mean, md_std, md_min, md_max)
!         real(DP), intent(in) :: samples(:, :), mean_vec(6), cov_mat(6, 6)
!         integer, intent(in) :: n
!         real(DP), intent(out) :: skew(6), ekurt(6), md_mean, md_std, md_min, md_max
!         real(DP) :: diff(6), m2(6), m3(6), m4(6), L(6, 6), y(6), md
!         integer :: i, j, k, info

!         skew = 0.0_DP; ekurt = 0.0_DP; m2 = 0.0_DP; m3 = 0.0_DP; m4 = 0.0_DP
!         do i = 1, n
!             diff = samples(:, i) - mean_vec
!             do j = 1, 6
!                 m2(j) = m2(j) + diff(j)**2
!                 m3(j) = m3(j) + diff(j)**3
!                 m4(j) = m4(j) + diff(j)**4
!             end do
!         end do
!         do j = 1, 6
!             skew(j)  = (m3(j) / n) / (m2(j) / n)**1.5_DP
!             ekurt(j) = (m4(j) / n) / (m2(j) / n)**2 - 3.0_DP
!         end do

!         ! Cholesky 分解
!         L = cov_mat
!         call dpotrf('L', 6, L, 6, info)
!         if (info /= 0) then
!             md_mean = -1.0_DP; md_std = -1.0_DP; md_min = -1.0_DP; md_max = -1.0_DP
!             return
!         end if
!         do i = 1, 6
!             L(1:i-1, i) = 0.0_DP
!         end do

!         md_mean = 0.0_DP; md_std = 0.0_DP
!         md_min = huge(1.0_DP); md_max = 0.0_DP
!         do i = 1, n
!             diff = samples(:, i) - mean_vec
!             do j = 1, 6
!                 y(j) = diff(j)
!                 do k = 1, j - 1
!                     y(j) = y(j) - L(j, k) * y(k)
!                 end do
!                 y(j) = y(j) / L(j, j)
!             end do
!             md = dot_product(y, y)
!             md_mean = md_mean + md
!             md_std  = md_std  + md * md
!             if (md < md_min) md_min = md
!             if (md > md_max) md_max = md
!         end do
!         md_mean = md_mean / n
!         md_std  = sqrt(max(md_std / n - md_mean**2, 0.0_DP))
!     end subroutine compute_gaussianity

! end program test_mc_non_gaussian
