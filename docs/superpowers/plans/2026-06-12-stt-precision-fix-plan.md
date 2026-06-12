# STT Precision Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix two STT bugs: (1) missing multinomial multiplicity weights in tensor contraction causing ~2.1e-4 mean mismatch, (2) slow DIAGNOSTIC path blocking the optimized order-based RHS computation.

**Architecture:** Two-file change: `pod_stt_tensor_module.f90` gets multiplicity tracking in `gen_block_assignments_v2` and `compute_stt_rhs_order`; `pod_uq_prop_stt_module.f90` replaces the diagnostic `compute_stt_rhs` loop with the batch-optimized `compute_stt_rhs_order`.

**Tech Stack:** Fortran 2008+, no new dependencies

---

### Task 1: Add multinomial multiplicity to gen_block_assignments_v2

**Files:**
- Modify: `POD_Fortran/src/lib/math/pod_stt_tensor_module.f90:997-1002` (signature)
- Modify: `POD_Fortran/src/lib/math/pod_stt_tensor_module.f90:1051` (insert multiplicity calc)

- [ ] **Step 1: Add mults_out parameter to gen_block_assignments_v2 signature**

Replace the subroutine header:

```fortran
    subroutine gen_block_assignments_v2(tup, p, part_sizes, k, blocks_out, n_assign)
        integer, intent(in)  :: tup(p), p
        integer, intent(in)  :: part_sizes(:)
        integer, intent(in)  :: k
        integer, intent(out) :: blocks_out(:,:)  ! (:, k), 调用者提供足够空间
        integer, intent(out) :: n_assign
```

with:

```fortran
    subroutine gen_block_assignments_v2(tup, p, part_sizes, k, blocks_out, n_assign, mults_out)
        integer, intent(in)  :: tup(p), p
        integer, intent(in)  :: part_sizes(:)
        integer, intent(in)  :: k
        integer, intent(out) :: blocks_out(:,:)  ! (:, k), 调用者提供足够空间
        integer, intent(out) :: n_assign
        integer, intent(out) :: mults_out(:)     ! multiplicity per assignment
```

- [ ] **Step 2: Insert multiplicity calculation in enum_dist_v2**

After line 1051 (`bpos = bpos + 1`), before the block-index construction comment, insert:

```fortran
                bpos = bpos + 1
                ! 计算多重性: Πᵢ cnts(i)! / Π_{blk} gba_dist(i, blk)!
                block
                    integer :: mult, fac_i, mb
                    mult = 1
                    do i = 1, uniq_n
                        fac_i = factorial(cnts(i))
                        do mb = 1, kk
                            fac_i = fac_i / factorial(gba_dist(i, mb))
                        end do
                        mult = mult * fac_i
                    end do
                    mults_out(bpos) = mult
                end block
```

Note: `i` variable is already declared in the host `gen_block_assignments_v2` scope (line 1005). The `block` construct scopes the new locals `mult`, `fac_i`, `mb`.

- [ ] **Step 3: Commit**

```bash
git add POD_Fortran/src/lib/math/pod_stt_tensor_module.f90
git commit -m "fix: add multinomial multiplicity to gen_block_assignments_v2

Add mults_out(:) output parameter storing the combinatorial
multiplicity for each block assignment, computed as the product
of multinomial coefficients over unique tuple values.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: Use multiplicity weights in compute_stt_rhs_order

**Files:**
- Modify: `POD_Fortran/src/lib/math/pod_stt_tensor_module.f90:936-989`

- [ ] **Step 1: Allocate blk_mults array alongside blk_assigns**

Replace lines 936-943:
```fortran
        integer, allocatable :: blk_assigns(:,:,:)  ! (partition, n_assign, k)
        integer, allocatable :: blk_counts(:)       ! (partition) n_assign
        integer, parameter :: MAX_ASSIGN_PER_PART = 2000

        np = stt_sizes(p)

        associate(pl => part_cache(p))
            allocate(blk_assigns(pl%n_parts, MAX_ASSIGN_PER_PART, p))
            allocate(blk_counts(pl%n_parts))
```

with:
```fortran
        integer, allocatable :: blk_assigns(:,:,:)  ! (partition, n_assign, k)
        integer, allocatable :: blk_mults(:,:)      ! (partition, n_assign) multiplicity
        integer, allocatable :: blk_counts(:)       ! (partition) n_assign
        integer, parameter :: MAX_ASSIGN_PER_PART = 2000

        np = stt_sizes(p)

        associate(pl => part_cache(p))
            allocate(blk_assigns(pl%n_parts, MAX_ASSIGN_PER_PART, p))
            allocate(blk_mults(pl%n_parts, MAX_ASSIGN_PER_PART))
            allocate(blk_counts(pl%n_parts))
```

- [ ] **Step 2: Pass blk_mults to gen_block_assignments_v2**

Replace lines 954-957:
```fortran
                    call gen_block_assignments_v2(tup, p, &
                        pl%counts(ip, 1:k_blk), k_blk, &
                        blk_assigns(ip, 1:MAX_ASSIGN_PER_PART, 1:k_blk), &
                        blk_counts(ip))
```

with:
```fortran
                    call gen_block_assignments_v2(tup, p, &
                        pl%counts(ip, 1:k_blk), k_blk, &
                        blk_assigns(ip, 1:MAX_ASSIGN_PER_PART, 1:k_blk), &
                        blk_counts(ip), &
                        blk_mults(ip, 1:MAX_ASSIGN_PER_PART))
```

- [ ] **Step 3: Multiply term by multiplicity in accumulation**

Replace line 981:
```fortran
                            dphi = dphi + term / real(eq_fac, DP)
```

with:
```fortran
                            dphi = dphi + term * real(blk_mults(ip, ia), DP) / real(eq_fac, DP)
```

- [ ] **Step 4: Deallocate blk_mults**

Replace line 989:
```fortran
            deallocate(blk_assigns, blk_counts)
```

with:
```fortran
            deallocate(blk_assigns, blk_mults, blk_counts)
```

- [ ] **Step 5: Commit**

```bash
git add POD_Fortran/src/lib/math/pod_stt_tensor_module.f90
git commit -m "fix: apply multinomial multiplicity weights in compute_stt_rhs_order

Multiply each non-uniform contraction term by the block assignment's
combinatorial multiplicity to correct under-counting in the Faa di Bruno
tensor contraction. This fixes the mean mismatch (~2.1e-4) and the
STT covariance diagonal deviation in component 5.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Remove DIAGNOSTIC path, enable optimized RHS in compute_stt_aug_derivatives

**Files:**
- Modify: `POD_Fortran/src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90:137-195`

- [ ] **Step 1: Remove DIAGNOSTIC local variables**

Replace lines 137-139:
```fortran
        ! DIAGNOSTIC: 旧版 RHS 路径临时变量
        type(stt_store_type) :: stt_tmp
        integer :: pp, np2, cidx2, tup2(6)
```

with nothing (delete all 3 lines plus the blank line before them). The resulting code should have just the blank line before `mu_ = mu`.

- [ ] **Step 2: Replace DIAGNOSTIC block with optimized loop**

Replace lines 169-195:
```fortran
        ! ---- 5. 高阶 STT RHS (order 2..max_order) ----
        ! DIAGNOSTIC: 使用旧版 compute_stt_rhs 逐分量路径
        ! 从 y(7:) 构建 stt_store (所有阶数 1..order)
        call stt_tmp%init(order)
        do pp = 1, order
            np2 = stt_sizes(pp)
            do i_comp = 1, 6
                do cidx2 = 1, np2
                    call stt_tmp%set(i_comp, cidx2, pp, &
                        y(7 + stt_flat_offset(i_comp, cidx2, pp) - 1))
                end do
            end do
        end do

        ! 使用旧版 compute_stt_rhs 逐分量计算 dydt
        do p = 2, order
            np2 = stt_sizes(p)
            do i_comp = 1, 6
                do cidx2 = 1, np2
                    call sym_index_to_tuple(cidx2, p, tup2(1:p))
                    call compute_stt_rhs(i_comp, p, tup2(1:p), &
                        f_tensors, order, stt_tmp, dydt(pos))
                    pos = pos + 1
                end do
            end do
        end do

        call stt_tmp%destroy()
```

with:
```fortran
        ! ---- 5. 高阶 STT RHS (order 2..max_order) ----
        do p = 2, order
            call compute_stt_rhs_order(p, y(7:), f_tensors, order, dydt, pos)
            pos = pos + 6 * stt_sizes(p)
        end do
```

Note: `pos`, `p`, `f_tensors`, `order`, `dydt` already exist as local variables in the subroutine (declared at line 135). No new variable declarations needed.

- [ ] **Step 3: Verify no unused imports**

The `stt_store_type` import at line 27 is still used by `stt_propagate` (line 236, 293, 321) and `stt_propagate_deviates` (line 590, 626, 639). The `compute_stt_rhs` import at line 25 is no longer called directly in this module, but `compute_stt_rhs_order` is already imported at line 26. No import changes needed — the compiler will silently accept the unused `compute_stt_rhs` import, or it can be removed from the `only` list if warnings are a concern.

- [ ] **Step 4: Commit**

```bash
git add POD_Fortran/src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90
git commit -m "fix: replace DIAGNOSTIC STT RHS path with optimized compute_stt_rhs_order

Remove the slow per-component-per-tuple diagnostic loop that packed
STT data into a stt_store object each step and replace with the batch
order-based RHS computation that shares block assignments across all
6 components and uses direct flat-array access.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: Rebuild and run teststt

- [ ] **Step 1: Rebuild the project**

```bash
cd POD_Fortran && make
```

- [ ] **Step 2: Run teststt**

```bash
cd POD_Fortran && ./teststt
```

Expected: all three checks pass:
```
PASS: mean matches (max diff < threshold)
PASS: STT covariance matches DA (rel err < 5%)
PASS: DA covariance matches MC (rel err < 10%)
```

- [ ] **Step 3: Verify specific improvements**

Check that:
- Mean max diff drops from ~2.1e-4 to well below 1e-6
- STT cov diag component 5 (~3.04e-7 before) rises to match DA value (~5.19e-7)
- No NaN or overflow in any output
