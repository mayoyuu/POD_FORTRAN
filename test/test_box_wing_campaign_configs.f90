program test_box_wing_campaign_configs
    use pod_global, only: DP
    use pod_config, only: config, set_default_config, load_config, validate_config
    implicit none

    call load_and_check('config/config_box_wing_earth.txt','earth')
    call load_and_check('config/config_box_wing_moon.txt','moon')
    write(*,'(a)') 'PASS: box-wing campaign configuration files.'

contains

    subroutine load_and_check(path,expected_mode)
        character(len=*),intent(in) :: path,expected_mode
        call set_default_config()
        call load_config(path)
        call assert_true(trim(config%srp_model)=='box_wing','box-wing model')
        call assert_true(trim(config%srp_attitude_mode)==expected_mode,'attitude mode')
        call assert_scalar(config%srp_mass_kg,1200.0_DP,'mass')
        call assert_close(config%srp_box_dimensions_m,[2.0_DP,3.0_DP,4.0_DP],'box dimensions')
        call assert_close(config%srp_box_optical,[0.30_DP,0.40_DP,0.30_DP],'box optical')
        call assert_scalar(config%srp_array_total_area_m2,24.0_DP,'array area')
        call assert_true(trim(config%srp_array_tracking_mode)=='single_axis','array mode')
        call assert_close(config%srp_array_hinge_axis_body,[0.0_DP,1.0_DP,0.0_DP],'hinge axis')
        call assert_close(config%srp_array_reference_normal_body,[1.0_DP,0.0_DP,0.0_DP], &
            'reference normal')
        call assert_close(config%srp_array_front_optical,[0.85_DP,0.08_DP,0.07_DP], &
            'front optical')
        call assert_close(config%srp_array_back_optical,[0.60_DP,0.20_DP,0.20_DP], &
            'back optical')
        call assert_close(config%srp_primary_axis_body,[0.0_DP,0.0_DP,1.0_DP],'primary axis')
        call assert_close(config%srp_secondary_axis_body,[0.0_DP,1.0_DP,0.0_DP], &
            'secondary axis')
        call assert_true(trim(config%srp_roll_reference)=='orbit_normal','roll reference')
        call assert_scalar(config%srp_pressure_1au_n_m2, &
            1367.0_DP/299792458.0_DP,'solar pressure')
        call assert_scalar(config%srp_scale_da_span,0.0_DP,'scale span')
        call assert_close(config%srp_attitude_bias_span_arcsec,[0.0_DP,0.0_DP,0.0_DP], &
            'attitude span')
        call assert_scalar(config%srp_array_angle_span_deg,0.0_DP,'array span')
        call assert_true(validate_config(),'valid configuration')
    end subroutine load_and_check

    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_scalar(actual,expected,label)
        real(DP),intent(in) :: actual,expected
        character(len=*),intent(in) :: label
        if(abs(actual-expected)>1.0e-14_DP*max(1.0_DP,abs(expected))) then
            write(*,'(a,2(1x,es24.16))') 'FAIL: '//trim(label),actual,expected
            error stop 1
        end if
    end subroutine assert_scalar

    subroutine assert_close(actual,expected,label)
        real(DP),intent(in) :: actual(3),expected(3)
        character(len=*),intent(in) :: label
        if(maxval(abs(actual-expected))>1.0e-14_DP) then
            write(*,'(a,6(1x,es24.16))') 'FAIL: '//trim(label),actual,expected
            error stop 1
        end if
    end subroutine assert_close
end program test_box_wing_campaign_configs
