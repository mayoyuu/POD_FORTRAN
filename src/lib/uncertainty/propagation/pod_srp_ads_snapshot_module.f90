!> Reproducible sparse ADS patch snapshots at selected checkpoints.
module pod_srp_ads_snapshot_module
    use pod_global, only: DP
    use pod_ads_split_module, only: manifold_type, sh_center, sh_width
    use pod_uq_hfem_ads_module, only: hfem_ads_history_type
    implicit none
    private
    public :: write_ads_patch_snapshot
contains
    subroutine write_ads_patch_snapshot(history, prefix, use_previous, status, message)
        type(hfem_ads_history_type), target, intent(in) :: history
        character(len=*), intent(in) :: prefix
        logical, intent(in) :: use_previous
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(manifold_type), pointer :: domain
        real(DP) :: time_seconds
        real(DP), allocatable :: center(:), width(:)
        integer :: umeta, ucoeff, uhist, umap, io, p, k, j, powers(7), nsplit
        character(len=64) :: label

        status=0
        message=''
        if (.not.history%initialized) then
            status=1
            message='ADS history is not initialized'
            return
        end if
        if (use_previous) then
            if (.not.history%previous_available) then
                status=2
                message='previous checkpoint is unavailable'
                return
            end if
            domain=>history%previous_domain
            time_seconds=history%previous_time
        else
            domain=>history%domain
            time_seconds=history%current_time
        end if
        open(newunit=umeta,file=trim(prefix)//'_meta.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 900
        open(newunit=ucoeff,file=trim(prefix)//'_coeff.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 901
        open(newunit=uhist,file=trim(prefix)//'_splits.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 902
        open(newunit=umap,file=trim(prefix)//'_map.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 903
        write(umeta,'(A)') 'patch_id,time_seconds,center1,center2,center3,center4,center5,center6,center7,'// &
            'width1,width2,width3,width4,width5,width6,width7,split_count'
        write(ucoeff,'(A)') 'patch_id,component,power1,power2,power3,power4,power5,power6,power7,coefficient'
        write(uhist,'(A)') 'patch_id,sequence,split_direction'
        write(umap,'(A)') 'key,i,j,value'
        write(umap,'(A)') 'state_unit,0,0,km_and_km_per_s'
        write(umap,'(A)') 'polynomial_variables,0,0,local_patch_coordinates_minus_one_to_one'
        write(umap,'(A)') 'coordinate_mode,0,0,whitened'
        write(umap,'(A,I0)') 'da_order,0,0,',history%options%da_order
        write(umap,'(A,ES25.16E3)') 'domain_sigma,0,0,',history%coordinate_map%domain_sigma
        write(umap,'(A,ES25.16E3)') 'eta_sigma,0,0,',history%coordinate_map%srp_sigma
        do j=1,6
            write(umap,'(A,I0,A,ES25.16E3)') 'mean,',j,',0,',history%coordinate_map%mean6(j)
            do k=1,6
                write(umap,'(A,I0,A,I0,A,ES25.16E3)') 'basis,',j,',',k,',',history%coordinate_map%basis6(j,k)
            end do
        end do
        do p=1,domain%n_patches
            center=sh_center(domain%patches(p)%history)
            width=sh_width(domain%patches(p)%history)
            nsplit=0
            if(allocated(domain%patches(p)%history%entries)) nsplit=size(domain%patches(p)%history%entries)
            write(umeta,'(I0,",",ES25.16E3,14(",",ES25.16E3),",",I0)') &
                p,time_seconds,center,width,nsplit
            do j=1,nsplit
                write(uhist,'(I0,",",I0,",",I0)') p,j,domain%patches(p)%history%entries(j)
            end do
            do k=1,6
                powers=0
                call emit_degree(0,1,history%options%da_order)
            end do
        end do
        close(umap)
        close(uhist)
        close(ucoeff)
        close(umeta)
        return
900     status=10
        message='cannot open patch metadata output'
        return
901     close(umeta)
        status=11
        message='cannot open patch coefficients output'
        return
902     close(ucoeff)
        close(umeta)
        status=12
        message='cannot open patch split output'
        return
903     close(uhist)
        close(ucoeff)
        close(umeta)
        status=13
        message='cannot open patch map output'
        return
    contains
        recursive subroutine emit_degree(total,index,max_degree)
            integer,intent(in) :: total,index,max_degree
            integer :: exponent
            real(DP) :: coeff
            if(index==8) then
                if(total>max_degree) return
                coeff=domain%patches(p)%da_vec%elements(k)%get_coeff(powers)
                if(coeff/=0.0_DP) write(ucoeff,'(I0,",",I0,7(",",I0),",",ES25.16E3)') &
                    p,k,powers,coeff
                return
            end if
            do exponent=0,max_degree-total
                powers(index)=exponent
                call emit_degree(total+exponent,index+1,max_degree)
            end do
        end subroutine emit_degree
    end subroutine write_ads_patch_snapshot
end module pod_srp_ads_snapshot_module
