!--------------------------------------------------------------------------------------------------------------
!> STT 对称张量模块
!>
!> 提供 STT 计算所需的对称张量数据结构和组合算法：
!>   - 对称索引压缩映射 (词典序)
!>   - 整数划分生成
!>   - Faà di Bruno 组合引擎
!>   - 非线性矩计算
!>
!> ## 数学背景
!>
!> 对称张量 Φ_{i, k₁...kₚ} 的指标 (k₁,...,kₚ) 在排列下不变。
!> 利用可重复组合将 C(6+p-1, p) 个独立分量压缩存储。
!>
!> ## 依赖
!> - pod_global: DP 精度类型
!--------------------------------------------------------------------------------------------------------------
module pod_stt_tensor_module
    use pod_global, only: DP
    implicit none
    private

    ! =========================================================================
    ! 公共常量
    ! =========================================================================
    integer, parameter, public :: STT_MAX_ORDER = 6
    integer, parameter, public :: STT_DIM = 6

    ! =========================================================================
    ! 模块级存储 (初始化后只读)
    ! =========================================================================
    integer, allocatable, save :: binom_table(:,:)
    integer, save :: stt_sizes(0:STT_MAX_ORDER)
    logical, save :: indexing_initialized = .false.

    ! =========================================================================
    ! 公开接口
    ! =========================================================================
    public :: STT_MAX_ORDER, STT_DIM, stt_sizes
    public :: init_stt_indexing, stt_order_p_size, indexing_initialized
    public :: tuple_to_sym_index, sym_index_to_tuple
    public :: generate_partitions, factorial
    public :: compute_stt_rhs, compute_stt_moments
    public :: stt_store_type
    public :: gen_block_assignments, expand_sym_index

    ! =========================================================================
    ! 辅助类型: 划分存储
    ! =========================================================================
    type, public :: partition_list_type
        integer :: n_parts
        integer, allocatable :: counts(:,:)   ! (n_parts, 6) — zero-padded
        integer, allocatable :: k(:)          ! (n_parts) — actual number of blocks
        integer, allocatable :: mults(:)      ! (n_parts) — combinatorial multiplicity
    end type partition_list_type

    ! =========================================================================
    ! 辅助类型: 块指配存储
    ! =========================================================================
    type, public :: block_assign_type
        integer :: n_assign
        integer :: k
        integer, allocatable :: blocks(:,:)   ! (n_assign, k) — compressed block indices
    end type block_assign_type

    ! =========================================================================
    ! STT 存储类型
    ! =========================================================================
    type, public :: stt_store_type
        integer :: order = 0
        real(DP), allocatable :: stt(:,:,:)   ! (6, max_size, 0:order)
    contains
        procedure :: init => stt_store_init
        procedure :: get  => stt_store_get
        procedure :: set  => stt_store_set
        procedure :: get_ptr => stt_store_get_ptr
        procedure :: destroy => stt_store_destroy
    end type stt_store_type

contains

    ! =========================================================================
    ! 阶乘 (查表)
    ! =========================================================================
    integer function factorial(n) result(res)
        integer, intent(in) :: n
        integer, parameter :: ft(0:12) = [1, 1, 2, 6, 24, 120, 720, 5040, 40320, &
                                           362880, 3628800, 39916800, 479001600]
        if (n < 0) then
            res = 1
        else if (n <= 12) then
            res = ft(n)
        else
            res = ft(12)
        end if
    end function factorial

    ! =========================================================================
    ! 二项式系数 (仅在建表时使用)
    ! =========================================================================
    integer function binom_direct(n, k) result(res)
        integer, intent(in) :: n, k
        integer :: i, kk
        if (k < 0 .or. k > n) then
            res = 0
            return
        end if
        kk = min(k, n - k)
        res = 1
        do i = 1, kk
            res = res * (n - i + 1) / i
        end do
    end function binom_direct

    ! =========================================================================
    ! 初始化索引系统 (调用一次)
    ! =========================================================================
    subroutine init_stt_indexing(max_order)
        integer, intent(in) :: max_order
        integer :: n, k, p

        if (indexing_initialized) return

        allocate(binom_table(0:STT_DIM + STT_MAX_ORDER, 0:STT_DIM + STT_MAX_ORDER))
        do n = 0, STT_DIM + STT_MAX_ORDER
            do k = 0, n
                binom_table(n, k) = binom_direct(n, k)
            end do
        end do

        stt_sizes(0) = 1
        do p = 1, STT_MAX_ORDER
            stt_sizes(p) = binom_table(STT_DIM + p - 1, p)
        end do

        indexing_initialized = .true.
    end subroutine init_stt_indexing

    ! =========================================================================
    ! 查询阶数 p 对应的对称张量大小
    ! =========================================================================
    integer function stt_order_p_size(p) result(sz)
        integer, intent(in) :: p
        if (p >= 0 .and. p <= STT_MAX_ORDER) then
            sz = stt_sizes(p)
        else
            sz = 0
        end if
    end function stt_order_p_size

    ! =========================================================================
    ! 词典序: 非降序元组 → 压缩索引 (1-based)
    !
    ! 对于 tup(1) ≤ tup(2) ≤ ... ≤ tup(p)，元素 ∈ {1..6}
    ! 压缩索引 = 1 + Σ_{j=1}^{p} Σ_{t=prev}^{tup(j)-1} C(6 + p - j - t, p - j)
    ! =========================================================================
    integer function tuple_to_sym_index(tup, p) result(idx)
        integer, intent(in) :: tup(:), p
        integer :: j, t, prev

        idx = 1
        prev = 1
        do j = 1, p
            do t = prev, tup(j) - 1
                idx = idx + binom_table(STT_DIM + p - j - t, p - j)
            end do
            prev = tup(j)
        end do
    end function tuple_to_sym_index

    ! =========================================================================
    ! 逆映射: 压缩索引 → 非降序元组
    ! =========================================================================
    subroutine sym_index_to_tuple(idx, p, tup)
        integer, intent(in)  :: idx, p
        integer, intent(out) :: tup(p)
        integer :: j, t, rem, cum, cnt

        rem = idx - 1
        cum = 1
        do j = 1, p
            do t = cum, STT_DIM
                cnt = binom_table(STT_DIM + p - j - t, p - j)
                if (rem < cnt) then
                    tup(j) = t
                    cum = t
                    exit
                end if
                rem = rem - cnt
            end do
        end do
    end subroutine sym_index_to_tuple

    ! =========================================================================
    ! 展开压缩索引为全排列元组集合
    !
    ! 给定压缩索引对应的对称元组，生成所有 6^p 种排列中属于该对称类的元组。
    ! 返回: expanded_tuples(n_perm, p) — 所有排列
    ! =========================================================================
    subroutine expand_sym_index(idx, p, expanded_tuples, n_perm)
        integer, intent(in)  :: idx, p
        integer, allocatable, intent(out) :: expanded_tuples(:,:)
        integer, intent(out) :: n_perm
        integer :: tup(p), i, j
        integer, parameter :: MAX_PERM = 6**6  ! worst case, will rarely reach

        call sym_index_to_tuple(idx, p, tup)

        ! 生成所有非降序排列 = 所有保持可重复性的重排
        ! 对于对称元组，所有不同的赋值对应于把每个不同值分配到可用位置上
        ! 简化处理: 生成所有 6^p 排列，再选出映射到同一压缩索引的
        ! 这对于 p ≤ 6 是可接受的 (6^6 = 46656)
        allocate(expanded_tuples(MAX_PERM, p))
        n_perm = 0
        call gen_all_assignments_rec(tup, p, 1, expanded_tuples, n_perm)
    end subroutine expand_sym_index

    recursive subroutine gen_all_assignments_rec(tup, p, pos, result, count)
        integer, intent(in)    :: tup(p), p, pos
        integer, intent(inout) :: result(:,:)
        integer, intent(inout) :: count
        integer, save :: current(p)
        integer :: val, i

        if (pos > p) then
            count = count + 1
            result(count, 1:p) = current(1:p)
            return
        end if

        do val = 1, STT_DIM
            current(pos) = val
            call gen_all_assignments_rec(tup, p, pos + 1, result, count)
        end do
    end subroutine gen_all_assignments_rec

    ! =========================================================================
    ! 生成整数 p 的所有划分 (partition)
    !
    ! 输出: plist — 划分列表
    !   counts(i, 1:k(i)) = (λ₁, λ₂, ..., λₖ) with λ₁ ≥ λ₂ ≥ ... ≥ λₖ
    !   mults(i) = p! / (Π λⱼ! · Π (cnt_of_each_equal_λ)!)
    ! =========================================================================
    subroutine generate_partitions(p, plist)
        integer, intent(in) :: p
        type(partition_list_type), intent(out) :: plist
        integer, parameter :: MAX_PARTS = 32
        integer :: stack(6)
        integer :: i, d, cnt, depth, part_val, rem

        plist%n_parts = 0
        allocate(plist%counts(MAX_PARTS, 6), plist%k(MAX_PARTS), plist%mults(MAX_PARTS))
        plist%counts = 0
        plist%k = 0
        plist%mults = 0

        ! 使用迭代栈模拟递归，避免 save 变量问题
        call rec_part(p, p, 0, stack, plist)

    contains
        recursive subroutine rec_part(rem, max_part, depth, stack, plist)
            integer, intent(in)    :: rem, max_part, depth
            integer, intent(inout) :: stack(6)
            type(partition_list_type), intent(inout) :: plist
            integer :: part_val, d, cnt

            if (rem == 0) then
                plist%n_parts = plist%n_parts + 1
                plist%k(plist%n_parts) = depth
                plist%counts(plist%n_parts, 1:depth) = stack(1:depth)

                ! 多重性: p! / (Π λᵢ! · Π (cnt_of_eq_λ)!)
                plist%mults(plist%n_parts) = factorial(p)
                do d = 1, depth
                    plist%mults(plist%n_parts) = &
                        plist%mults(plist%n_parts) / factorial(stack(d))
                end do
                ! 除以相同划分值的排列数
                d = 1
                do while (d <= depth)
                    cnt = 1
                    do while (d + cnt <= depth .and. stack(d + cnt) == stack(d))
                        cnt = cnt + 1
                    end do
                    plist%mults(plist%n_parts) = &
                        plist%mults(plist%n_parts) / factorial(cnt)
                    d = d + cnt
                end do
                return
            end if

            do part_val = min(max_part, rem), 1, -1
                stack(depth + 1) = part_val
                call rec_part(rem - part_val, part_val, depth + 1, stack, plist)
            end do
        end subroutine rec_part
    end subroutine generate_partitions

    ! =========================================================================
    ! 生成块指配: 将元组的 p 个指标分配到 k 个块中
    !
    ! 给定一个非降序元组 tup(p) 和划分 (λ₁,...,λₖ),
    ! 生成所有不同的块指配。每个指配产生 k 个压缩索引 (B₁,...,Bₖ)。
    !
    ! 方法: 先将 tup 中相同值的元素合并为 (value, count) 对，
    !       然后枚举每个值在 k 个块中的分布方式。
    ! =========================================================================
    subroutine gen_block_assignments(tup, p, part_sizes, k, blocks_list, n_assign)
        integer, intent(in)  :: tup(p), p
        integer, intent(in)  :: part_sizes(k)
        integer, intent(in)  :: k
        integer, allocatable, intent(out) :: blocks_list(:,:)
        integer, intent(out) :: n_assign

        integer :: values(6), counts(6), uniq_n
        integer :: i, j, v, prev
        integer, allocatable :: dist(:,:)    ! (uniq_n, k) — 每个值的分布
        integer, allocatable :: block_tuples(:,:,:)  ! (k, p, temporary)
        integer, allocatable :: buf(:,:)
        integer, parameter :: MAX_ASSIGN = 10000
        integer :: buf_size

        allocate(buf(MAX_ASSIGN, k))
        buf = 0
        buf_size = 0

        ! 1. 压缩 tup 为 (value, count) 对
        uniq_n = 0
        i = 1
        do while (i <= p)
            uniq_n = uniq_n + 1
            values(uniq_n) = tup(i)
            counts(uniq_n) = 1
            do while (i + counts(uniq_n) <= p .and. tup(i + counts(uniq_n)) == tup(i))
                counts(uniq_n) = counts(uniq_n) + 1
            end do
            i = i + counts(uniq_n)
        end do

        ! 2. 为每个 (value, count) 枚举在 k 块中的分布
        !    使用递推: 对当前值，枚举将其 count 个元素分配到 k 块的方案
        !    (每个块获得的元素数必须非负且总和 = count)
        allocate(dist(uniq_n, k))
        allocate(block_tuples(k, p, 1))  ! placeholder, 将在递归中填充
        call enum_distributions(1)

        ! 3. 复制结果
        n_assign = buf_size
        allocate(blocks_list(n_assign, k))
        blocks_list(1:n_assign, 1:k) = buf(1:n_assign, 1:k)

        deallocate(dist, block_tuples, buf)

    contains
        recursive subroutine enum_distributions(vi)
            integer, intent(in) :: vi
            integer :: d1, d2, d3, d4, d5, d6
            integer :: c, blk, total

            if (vi > uniq_n) then
                ! 所有值已分配完毕，计算各块的压缩索引
                call build_block_indices()
                return
            end if

            c = counts(vi)

            if (k == 1) then
                dist(vi, 1) = c
                call enum_distributions(vi + 1)
            else if (k == 2) then
                do d1 = 0, c
                    dist(vi, 1) = d1
                    dist(vi, 2) = c - d1
                    call enum_distributions(vi + 1)
                end do
            else if (k == 3) then
                do d1 = 0, c
                    do d2 = 0, c - d1
                        dist(vi, 1) = d1
                        dist(vi, 2) = d2
                        dist(vi, 3) = c - d1 - d2
                        call enum_distributions(vi + 1)
                    end do
                end do
            else if (k == 4) then
                do d1 = 0, c
                    do d2 = 0, c - d1
                        do d3 = 0, c - d1 - d2
                            dist(vi, 1) = d1
                            dist(vi, 2) = d2
                            dist(vi, 3) = d3
                            dist(vi, 4) = c - d1 - d2 - d3
                            call enum_distributions(vi + 1)
                        end do
                    end do
                end do
            else if (k == 5) then
                do d1 = 0, c
                    do d2 = 0, c - d1
                        do d3 = 0, c - d1 - d2
                            do d4 = 0, c - d1 - d2 - d3
                                dist(vi, 1) = d1
                                dist(vi, 2) = d2
                                dist(vi, 3) = d3
                                dist(vi, 4) = d4
                                dist(vi, 5) = c - d1 - d2 - d3 - d4
                                call enum_distributions(vi + 1)
                            end do
                        end do
                    end do
                end do
            else if (k == 6) then
                do d1 = 0, c
                    do d2 = 0, c - d1
                        do d3 = 0, c - d1 - d2
                            do d4 = 0, c - d1 - d2 - d3
                                do d5 = 0, c - d1 - d2 - d3 - d4
                                    dist(vi, 1) = d1
                                    dist(vi, 2) = d2
                                    dist(vi, 3) = d3
                                    dist(vi, 4) = d4
                                    dist(vi, 5) = d5
                                    dist(vi, 6) = c - d1 - d2 - d3 - d4 - d5
                                    call enum_distributions(vi + 1)
                                end do
                            end do
                        end do
                    end do
                end do
            end if
        end subroutine enum_distributions

        subroutine build_block_indices()
            integer :: blk, vi, v, di, idx
            integer :: block_tup(p)
            integer :: block_order, bp, blk_size
            integer :: block_idx

            ! 检查每个块的分配是否符合其应有的大小
            do blk = 1, k
                blk_size = 0
                do vi = 1, uniq_n
                    blk_size = blk_size + dist(vi, blk)
                end do
                if (blk_size /= part_sizes(blk)) return  ! 分配与划分不匹配
            end do

            ! 为每个块构建压缩索引
            buf_size = buf_size + 1
            if (buf_size > MAX_ASSIGN) then
                ! 扩展缓冲区
                return
            end if

            do blk = 1, k
                ! 收集块 blk 中的指标
                bp = 0
                do vi = 1, uniq_n
                    v = values(vi)
                    do di = 1, dist(vi, blk)
                        bp = bp + 1
                        block_tup(bp) = v
                    end do
                end do
                block_order = part_sizes(blk)
                if (block_order > 0) then
                    block_idx = tuple_to_sym_index(block_tup(1:block_order), block_order)
                    buf(buf_size, blk) = block_idx
                else
                    buf(buf_size, blk) = 0
                end if
            end do
        end subroutine build_block_indices
    end subroutine gen_block_assignments

    ! =========================================================================
    ! Faà di Bruno α-收缩核心
    !
    ! 对于给定的 k 个块，计算:
    !   Σ_{α₁,...,αₖ ∈ {1..6}^k} f*_{i, α₁...αₖ} · Π_{b=1}^{k} Φ_{α_b, B_b, λ_b}
    !
    ! 使用显式嵌套循环优化 k=1..6 的情况。
    ! =========================================================================
    real(DP) function contract_alpha(i_comp, k, block_indices, block_orders, &
                                      f_tensors, max_f_order, stt_store) result(res)
        integer, intent(in) :: i_comp, k
        integer, intent(in) :: block_indices(k)
        integer, intent(in) :: block_orders(k)
        real(DP), intent(in) :: f_tensors(:, :, 0:)
        integer, intent(in) :: max_f_order
        type(stt_store_type), intent(in) :: stt_store
        integer :: a1, a2, a3, a4, a5, a6
        integer :: fidx
        real(DP) :: f_val, prod

        res = 0.0_DP

        ! f* 的压缩索引: 将 α 元组压缩
        ! 注意: f* 的指标 α₁...αₖ 不排序 (因为它们是求和的虚拟指标)
        ! 但 f_tensors 存储的是对称压缩形式
        ! 对于 f* 的 RHS 使用，我们需要 f*_{i, α₁...αₖ} 其中指标已排序
        ! 由于 f_tensors 是压缩对称存储的，需要先排序 α 再查询

        select case (k)
        case (1)
            do a1 = 1, 6
                fidx = tuple_to_sym_index([a1], 1)
                f_val = f_tensors(i_comp, fidx, 1)
                prod = f_val * stt_store%get(a1, block_indices(1), block_orders(1))
                res = res + prod
            end do
        case (2)
            do a1 = 1, 6
                do a2 = 1, 6
                    fidx = tuple_to_sym_index(sort2(a1, a2), 2)
                    f_val = f_tensors(i_comp, fidx, 2)
                    prod = f_val &
                        * stt_store%get(a1, block_indices(1), block_orders(1)) &
                        * stt_store%get(a2, block_indices(2), block_orders(2))
                    res = res + prod
                end do
            end do
        case (3)
            do a1 = 1, 6
                do a2 = 1, 6
                    do a3 = 1, 6
                        fidx = tuple_to_sym_index(sort3(a1, a2, a3), 3)
                        f_val = f_tensors(i_comp, fidx, 3)
                        prod = f_val &
                            * stt_store%get(a1, block_indices(1), block_orders(1)) &
                            * stt_store%get(a2, block_indices(2), block_orders(2)) &
                            * stt_store%get(a3, block_indices(3), block_orders(3))
                        res = res + prod
                    end do
                end do
            end do
        case (4)
            do a1 = 1, 6
                do a2 = 1, 6
                    do a3 = 1, 6
                        do a4 = 1, 6
                            fidx = tuple_to_sym_index(sort4(a1, a2, a3, a4), 4)
                            f_val = f_tensors(i_comp, fidx, 4)
                            prod = f_val &
                                * stt_store%get(a1, block_indices(1), block_orders(1)) &
                                * stt_store%get(a2, block_indices(2), block_orders(2)) &
                                * stt_store%get(a3, block_indices(3), block_orders(3)) &
                                * stt_store%get(a4, block_indices(4), block_orders(4))
                            res = res + prod
                        end do
                    end do
                end do
            end do
        case (5)
            do a1 = 1, 6
                do a2 = 1, 6
                    do a3 = 1, 6
                        do a4 = 1, 6
                            do a5 = 1, 6
                                fidx = tuple_to_sym_index(sort5(a1, a2, a3, a4, a5), 5)
                                f_val = f_tensors(i_comp, fidx, 5)
                                prod = f_val &
                                    * stt_store%get(a1, block_indices(1), block_orders(1)) &
                                    * stt_store%get(a2, block_indices(2), block_orders(2)) &
                                    * stt_store%get(a3, block_indices(3), block_orders(3)) &
                                    * stt_store%get(a4, block_indices(4), block_orders(4)) &
                                    * stt_store%get(a5, block_indices(5), block_orders(5))
                                res = res + prod
                            end do
                        end do
                    end do
                end do
            end do
        case (6)
            do a1 = 1, 6
                do a2 = 1, 6
                    do a3 = 1, 6
                        do a4 = 1, 6
                            do a5 = 1, 6
                                do a6 = 1, 6
                                    fidx = tuple_to_sym_index( &
                                        sort6(a1, a2, a3, a4, a5, a6), 6)
                                    f_val = f_tensors(i_comp, fidx, 6)
                                    prod = f_val &
                                        * stt_store%get(a1, block_indices(1), block_orders(1)) &
                                        * stt_store%get(a2, block_indices(2), block_orders(2)) &
                                        * stt_store%get(a3, block_indices(3), block_orders(3)) &
                                        * stt_store%get(a4, block_indices(4), block_orders(4)) &
                                        * stt_store%get(a5, block_indices(5), block_orders(5)) &
                                        * stt_store%get(a6, block_indices(6), block_orders(6))
                                    res = res + prod
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end select
    end function contract_alpha

    ! ---- 内联排序辅助函数 (2~6 个元素) ----

    function sort2(a, b) result(r)
        integer, intent(in) :: a, b
        integer :: r(2)
        if (a <= b) then; r = [a, b]; else; r = [b, a]; end if
    end function

    function sort3(a, b, c) result(r)
        integer, intent(in) :: a, b, c
        integer :: r(3)
        r = [a, b, c]
        call isort3(r)
    end function

    function sort4(a, b, c, d) result(r)
        integer, intent(in) :: a, b, c, d
        integer :: r(4)
        r = [a, b, c, d]
        call isort4(r)
    end function

    function sort5(a, b, c, d, e) result(r)
        integer, intent(in) :: a, b, c, d, e
        integer :: r(5)
        r = [a, b, c, d, e]
        call isort5(r)
    end function

    function sort6(a, b, c, d, e, f) result(r)
        integer, intent(in) :: a, b, c, d, e, f
        integer :: r(6)
        r = [a, b, c, d, e, f]
        call isort6(r)
    end function

    subroutine isort3(a)
        integer, intent(inout) :: a(3)
        integer :: t
        if (a(1) > a(2)) then; t = a(1); a(1) = a(2); a(2) = t; end if
        if (a(2) > a(3)) then; t = a(2); a(2) = a(3); a(3) = t; end if
        if (a(1) > a(2)) then; t = a(1); a(1) = a(2); a(2) = t; end if
    end subroutine

    subroutine isort4(a)
        integer, intent(inout) :: a(4)
        integer :: t, i, j
        do i = 1, 3
            do j = i + 1, 4
                if (a(i) > a(j)) then; t = a(i); a(i) = a(j); a(j) = t; end if
            end do
        end do
    end subroutine

    subroutine isort5(a)
        integer, intent(inout) :: a(5)
        integer :: t, i, j
        do i = 1, 4
            do j = i + 1, 5
                if (a(i) > a(j)) then; t = a(i); a(i) = a(j); a(j) = t; end if
            end do
        end do
    end subroutine

    subroutine isort6(a)
        integer, intent(inout) :: a(6)
        integer :: t, i, j
        do i = 1, 5
            do j = i + 1, 6
                if (a(i) > a(j)) then; t = a(i); a(i) = a(j); a(j) = t; end if
            end do
        end do
    end subroutine

    ! =========================================================================
    ! compute_stt_rhs — Faà di Bruno 右端项计算
    !
    ! 给定分量 i、阶数 p、压缩索引对应的元组 tup(1:p)、
    ! 力场导数张量 f_tensors 和当前 STT 存储，
    ! 计算 dΦ_{i, tup} / dt。
    !
    ! 算法:
    !   1. 均匀项: Σ_α f*_{i,α} · Φ_{α, tup, p}
    !   2. 划分项: 对每个 λ ⊢ p (k ≥ 2):
    !        对每个块指配:
    !          dphi += contract_alpha(i, k, [B₁..Bₖ], [λ₁..λₖ], f*, stt)
    ! =========================================================================
    subroutine compute_stt_rhs(i_comp, p, tup, f_tensors, max_f_order, stt_store, dphi)
        integer, intent(in) :: i_comp, p
        integer, intent(in) :: tup(p)
        real(DP), intent(in) :: f_tensors(:, :, 0:)
        integer, intent(in) :: max_f_order
        type(stt_store_type), intent(in) :: stt_store
        real(DP), intent(out) :: dphi

        type(partition_list_type) :: plist
        integer, allocatable :: blocks_list(:,:)
        integer :: n_assign
        integer :: ip, ia
        integer :: tup_idx
        integer :: alpha, k_blk
        real(DP) :: term

        dphi = 0.0_DP

        ! ---- 1. 均匀项: Σ_α f*_{i,α} · Φ_{α, tup_idx, p} ----
        tup_idx = tuple_to_sym_index(tup, p)
        do alpha = 1, 6
            dphi = dphi + f_tensors(i_comp, alpha, 1) &
                        * stt_store%get(alpha, tup_idx, p)
        end do

        ! ---- 2. 非均匀项: 对所有划分 λ ⊢ p (k ≥ 2) ----
        if (p >= 2) then
            call generate_partitions(p, plist)

            do ip = 1, plist%n_parts
                k_blk = plist%k(ip)
                if (k_blk < 2) cycle  ! 均匀项已处理

                call gen_block_assignments(tup, p, plist%counts(ip, 1:k_blk), &
                                           k_blk, blocks_list, n_assign)

                do ia = 1, n_assign
                    term = contract_alpha(i_comp, k_blk, &
                                          blocks_list(ia, 1:k_blk), &
                                          plist%counts(ip, 1:k_blk), &
                                          f_tensors, max_f_order, stt_store)
                    dphi = dphi + term
                end do

                deallocate(blocks_list)
            end do
        end if
    end subroutine compute_stt_rhs

    ! =========================================================================
    ! compute_stt_moments — 占位 (Task 3 实现)
    ! =========================================================================
    subroutine compute_stt_moments(x_star, stt_store, P0, order, mean_out, cov_out)
        real(DP), intent(in) :: x_star(6), P0(6,6)
        type(stt_store_type), intent(in) :: stt_store
        integer, intent(in) :: order
        real(DP), intent(out) :: mean_out(6), cov_out(6,6)
        mean_out = x_star
        cov_out = 0.0_DP
    end subroutine compute_stt_moments

    ! =========================================================================
    ! STT Store Type 方法实现
    ! =========================================================================
    subroutine stt_store_init(this, order)
        class(stt_store_type), intent(inout) :: this
        integer, intent(in) :: order
        integer :: max_sz, p
        if (order < 0 .or. order > STT_MAX_ORDER) return
        this%order = order
        max_sz = 0
        do p = 1, order
            max_sz = max(max_sz, stt_sizes(p))
        end do
        max_sz = max(max_sz, 1)
        allocate(this%stt(6, max_sz, 0:order))
        this%stt = 0.0_DP
    end subroutine stt_store_init

    subroutine stt_store_destroy(this)
        class(stt_store_type), intent(inout) :: this
        if (allocated(this%stt)) deallocate(this%stt)
        this%order = 0
    end subroutine stt_store_destroy

    real(DP) function stt_store_get(this, icomp, cidx, ord) result(val)
        class(stt_store_type), intent(in) :: this
        integer, intent(in) :: icomp, cidx
        integer, intent(in), optional :: ord
        integer :: o
        o = 1
        if (present(ord)) o = ord
        if (allocated(this%stt) .and. o <= this%order .and. o >= 0) then
            if (cidx >= 1 .and. cidx <= size(this%stt, 2)) then
                val = this%stt(icomp, cidx, o)
            else
                val = 0.0_DP
            end if
        else
            val = 0.0_DP
        end if
    end function stt_store_get

    subroutine stt_store_set(this, icomp, cidx, ord, val)
        class(stt_store_type), intent(inout) :: this
        integer, intent(in) :: icomp, cidx, ord
        real(DP), intent(in) :: val
        if (allocated(this%stt) .and. ord <= this%order .and. ord >= 0) then
            if (cidx >= 1 .and. cidx <= size(this%stt, 2)) then
                this%stt(icomp, cidx, ord) = val
            end if
        end if
    end subroutine stt_store_set

    function stt_store_get_ptr(this, ord) result(ptr)
        class(stt_store_type), intent(in), target :: this
        integer, intent(in) :: ord
        real(DP), pointer :: ptr(:,:)
        if (allocated(this%stt) .and. ord <= this%order) then
            ptr => this%stt(:,:,ord)
        else
            nullify(ptr)
        end if
    end function stt_store_get_ptr

end module pod_stt_tensor_module
