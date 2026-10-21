!> Reproducible sparse ADS patch snapshots at selected checkpoints.
module pod_srp_ads_snapshot_module
    use pod_global, only: DP
    use pod_ads_split_module, only: manifold_type, sh_center, sh_width
    use pod_uq_hfem_ads_module, only: hfem_ads_history_type
    use pod_uq_ads_coordinates_module, only: ADS_PARAM_NONE,ADS_PARAM_SRP_ETA,ADS_PARAM_SRP_MODES
    implicit none
    private
    public :: write_ads_patch_snapshot
contains
    subroutine write_ads_patch_snapshot(history, prefix, use_previous, status, message, &
            schema_version,mode_csv_sha256,mode_metadata_sha256)
        type(hfem_ads_history_type), target, intent(in) :: history
        character(len=*), intent(in) :: prefix
        logical, intent(in) :: use_previous
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer,intent(in),optional :: schema_version
        character(len=*),intent(in),optional :: mode_csv_sha256,mode_metadata_sha256
        type(manifold_type), pointer :: domain
        real(DP) :: time_seconds
        real(DP), allocatable :: center(:), width(:)
        real(DP),allocatable :: saved_center(:),saved_width(:)
        integer,allocatable :: powers(:)
        integer :: umeta, ucoeff, uhist, umap, io, p, k, j, nsplit, nvars,columns,version
        character(len=64) :: label
        character(len=:),allocatable :: meta_header,coeff_header
        character(len=2),parameter :: roles(8)=[character(len=2) :: 'x','y','z','vx','vy','vz','s','u']

        status=0
        message=''
        version=1
        if(present(schema_version)) version=schema_version
        if(version/=1 .and. version/=2) then
            status=3; message='unsupported snapshot schema'; return
        end if
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
        nvars=history%da_variables
        if(version==1 .and. nvars>7) then
            status=3; message='8D snapshots require schema v2'; return
        end if
        columns=7
        if(version==2) columns=nvars
        allocate(saved_center(columns),saved_width(columns),powers(columns))
        meta_header='patch_id,time_seconds'; coeff_header='patch_id,component'
        do j=1,columns
            write(label,'(I0)') j
            meta_header=meta_header//',center'//trim(label)
            coeff_header=coeff_header//',power'//trim(label)
        end do
        do j=1,columns
            write(label,'(I0)') j
            meta_header=meta_header//',width'//trim(label)
        end do
        meta_header=meta_header//',split_count'; coeff_header=coeff_header//',coefficient'
        open(newunit=umeta,file=trim(prefix)//'_meta.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 900
        open(newunit=ucoeff,file=trim(prefix)//'_coeff.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 901
        open(newunit=uhist,file=trim(prefix)//'_splits.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 902
        open(newunit=umap,file=trim(prefix)//'_map.csv',status='replace',action='write',iostat=io)
        if(io/=0) goto 903
        write(umeta,'(A)') meta_header
        write(ucoeff,'(A)') coeff_header
        write(uhist,'(A)') 'patch_id,sequence,split_direction'
        write(umap,'(A)') 'key,i,j,value'
        write(umap,'(A)') 'state_unit,0,0,km_and_km_per_s'
        write(umap,'(A)') 'polynomial_variables,0,0,local_patch_coordinates_minus_one_to_one'
        write(umap,'(A)') 'coordinate_mode,0,0,whitened'
        if(version==2) then
            write(umap,'(A)') 'schema_version,0,0,2'
            write(umap,'(A,I0)') 'n_variables,0,0,',history%coordinate_map%n_variables
            write(umap,'(A,I0)') 'da_variables,0,0,',nvars
            select case(history%coordinate_map%parameterization)
            case(ADS_PARAM_NONE); label='NONE'
            case(ADS_PARAM_SRP_ETA); label='SRP_ETA'
            case(ADS_PARAM_SRP_MODES); label='SRP_MODES'
            end select
            write(umap,'(A)') 'parameterization,0,0,'//trim(label)
            do j=1,nvars
                label=roles(j)
                if(j>history%coordinate_map%n_variables) label='inactive'
                if(j==7 .and. history%coordinate_map%parameterization==ADS_PARAM_SRP_ETA) label='eta'
                write(umap,'(A,I0,A)') 'variable_role,',j,',0,'//trim(label)
            end do
            if(present(mode_csv_sha256)) write(umap,'(A)') 'mode_csv_sha256,0,0,'//trim(mode_csv_sha256)
            if(present(mode_metadata_sha256)) write(umap,'(A)') 'mode_metadata_sha256,0,0,'//trim(mode_metadata_sha256)
        end if
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
            saved_center=0.0_DP; saved_width=0.0_DP
            saved_center(1:nvars)=center; saved_width(1:nvars)=width
            nsplit=0
            if(allocated(domain%patches(p)%history%entries)) nsplit=size(domain%patches(p)%history%entries)
            write(umeta,'(I0,",",ES25.16E3,",")',advance='no') p,time_seconds
            write(umeta,'(*(ES25.16E3,","))',advance='no') saved_center,saved_width
            write(umeta,'(I0)') nsplit
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
            if(index==nvars+1) then
                if(total>max_degree) return
                coeff=domain%patches(p)%da_vec%elements(k)%get_coeff(powers(1:nvars))
                if(coeff/=0.0_DP) then
                    write(ucoeff,'(I0,",",I0,",")',advance='no') p,k
                    write(ucoeff,'(*(I0,","))',advance='no') powers
                    write(ucoeff,'(ES25.16E3)') coeff
                end if
                return
            end if
            do exponent=0,max_degree-total
                powers(index)=exponent
                call emit_degree(total+exponent,index+1,max_degree)
            end do
        end subroutine emit_degree
    end subroutine write_ads_patch_snapshot
end module pod_srp_ads_snapshot_module
