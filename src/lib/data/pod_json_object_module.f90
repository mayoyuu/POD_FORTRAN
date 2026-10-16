!> Small strict JSON reader for typed object members (no substring searches).
module pod_json_object_module
    use pod_global, only: DP
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: json_member, json_string, json_number, json_boolean
contains
    subroutine whitespace(text,p)
        character(len=*),intent(in) :: text
        integer,intent(inout) :: p
        do while(p<=len(text))
            if(index(' '//achar(9)//achar(10)//achar(13),text(p:p))==0) exit
            p=p+1
        end do
    end subroutine

    subroutine string_token(text,p,value,ok)
        character(len=*),intent(in) :: text
        integer,intent(inout) :: p
        character(len=:),allocatable,intent(out) :: value
        logical,intent(out) :: ok
        character :: c
        integer :: k,n,ios
        value=''; ok=.false.
        if(p>len(text)) return
        if(text(p:p)/='"') return
        p=p+1
        do while(p<=len(text))
            c=text(p:p); p=p+1
            if(c=='"') then
                ok=.true.; return
            end if
            if(iachar(c)<32) return
            if(c==achar(92)) then
                if(p>len(text)) return
                c=text(p:p); p=p+1
                select case(c)
                case('"','/',achar(92))
                case('b'); c=achar(8)
                case('f'); c=achar(12)
                case('n'); c=achar(10)
                case('r'); c=achar(13)
                case('t'); c=achar(9)
                case('u')
                    if(p+3>len(text)) return
                    do k=p,p+3
                        if(index('0123456789abcdefABCDEF',text(k:k))==0) return
                    end do
                    read(text(p:p+3),'(Z4)',iostat=ios) n
                    if(ios/=0) return
                    ! Non-ASCII escapes are valid JSON, but cannot equal ASCII identity keys.
                    c=achar(min(n,127)); p=p+4
                case default
                    return
                end select
            end if
            value=value//c
        end do
    end subroutine

    recursive subroutine scan_value(text,p,ok,depth)
        character(len=*),intent(in) :: text
        integer,intent(inout) :: p
        logical,intent(out) :: ok
        integer,intent(in) :: depth
        character(len=:),allocatable :: key,keys(:)
        character :: closing
        integer :: first,i,nkeys
        ok=.false.
        call whitespace(text,p)
        if(p>len(text) .or. depth>32) return
        select case(text(p:p))
        case('{','[')
            closing=']'
            if(text(p:p)=='{') closing='}'
            p=p+1; call whitespace(text,p)
            if(p>len(text)) return
            if(text(p:p)==closing) then
                p=p+1; ok=.true.; return
            end if
            allocate(character(len=len(text)) :: keys(0)); nkeys=0
            do
                if(closing=='}') then
                    call string_token(text,p,key,ok)
                    if(.not.ok) return
                    do i=1,nkeys
                        if(trim(keys(i))==key) then
                            ok=.false.; return
                        end if
                    end do
                    keys=[character(len=len(text)) :: keys,key]; nkeys=nkeys+1
                    call whitespace(text,p)
                    ok=.false.
                    if(p>len(text)) return
                    if(text(p:p)/=':') return
                    p=p+1
                end if
                call scan_value(text,p,ok,depth+1)
                if(.not.ok) return
                call whitespace(text,p); ok=.false.
                if(p>len(text)) return
                if(text(p:p)==closing) then
                    p=p+1; ok=.true.; return
                end if
                if(text(p:p)/=',') return
                p=p+1; call whitespace(text,p)
            end do
        case('"')
            call string_token(text,p,key,ok)
        case('t')
            if(p+3>len(text)) return
            if(text(p:p+3)/='true') return
            p=p+4; ok=.true.
        case('f')
            if(p+4>len(text)) return
            if(text(p:p+4)/='false') return
            p=p+5; ok=.true.
        case('n')
            if(p+3>len(text)) return
            if(text(p:p+3)/='null') return
            p=p+4; ok=.true.
        case default
            first=p
            if(text(p:p)=='-') p=p+1
            if(p>len(text)) return
            if(text(p:p)=='0') then
                p=p+1
            else
                if(index('123456789',text(p:p))==0) return
                call digits(text,p)
            end if
            if(p<=len(text)) then
                if(text(p:p)=='.') then
                    p=p+1; i=p; call digits(text,p)
                    if(p==i) return
                end if
            end if
            if(p<=len(text)) then
                if(index('eE',text(p:p))/=0) then
                    p=p+1
                    if(p<=len(text)) then
                        if(index('+-',text(p:p))/=0) p=p+1
                    end if
                    i=p; call digits(text,p)
                    if(p==i) return
                end if
            end if
            ok=p>first
        end select
    end subroutine

    subroutine digits(text,p)
        character(len=*),intent(in) :: text
        integer,intent(inout) :: p
        do while(p<=len(text))
            if(index('0123456789',text(p:p))==0) exit
            p=p+1
        end do
    end subroutine

    subroutine json_member(text,name,value,ok)
        character(len=*),intent(in) :: text,name
        character(len=:),allocatable,intent(out) :: value
        logical,intent(out) :: ok
        character(len=:),allocatable :: key
        integer :: p,first,last
        logical :: valid,found
        value=''; ok=.false.; p=1
        call scan_value(text,p,valid,0)
        if(.not.valid) return
        call whitespace(text,p)
        if(p/=len(text)+1) return
        p=1; call whitespace(text,p)
        if(text(p:p)/='{') return
        p=p+1; found=.false.
        do
            call whitespace(text,p)
            if(text(p:p)=='}') exit
            call string_token(text,p,key,valid)
            call whitespace(text,p); p=p+1
            call whitespace(text,p); first=p
            call scan_value(text,p,valid,1); last=p-1
            if(key==name) then
                value=text(first:last); found=.true.
            end if
            call whitespace(text,p)
            if(text(p:p)=='}') exit
            p=p+1
        end do
        ok=found
    end subroutine

    subroutine json_string(text,name,value,ok)
        character(len=*),intent(in) :: text,name
        character(len=:),allocatable,intent(out) :: value
        logical,intent(out) :: ok
        character(len=:),allocatable :: raw
        integer :: p
        call json_member(text,name,raw,ok)
        value=''
        if(.not.ok) return
        p=1; call string_token(raw,p,value,ok)
        if(p/=len(raw)+1) ok=.false.
    end subroutine

    subroutine json_number(text,name,value,ok)
        character(len=*),intent(in) :: text,name
        real(DP),intent(out) :: value
        logical,intent(out) :: ok
        character(len=:),allocatable :: raw
        integer :: ios
        call json_member(text,name,raw,ok)
        value=0.0_DP
        if(.not.ok) return
        if(index('-0123456789',raw(1:1))==0) then
            ok=.false.; return
        end if
        read(raw,*,iostat=ios) value
        ok=ios==0 .and. ieee_is_finite(value)
    end subroutine

    subroutine json_boolean(text,name,value,ok)
        character(len=*),intent(in) :: text,name
        logical,intent(out) :: value,ok
        character(len=:),allocatable :: raw
        call json_member(text,name,raw,ok)
        value=.false.
        if(.not.ok) return
        ok=raw=='true' .or. raw=='false'
        value=raw=='true'
    end subroutine
end module pod_json_object_module
