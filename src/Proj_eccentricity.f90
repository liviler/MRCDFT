Module Eccentricity
    use Constants,only: r64
    implicit none
    
contains

    subroutine calculate_Eccentri_n(n,iphi,it,Eccentri,pEccentri)
        !--------------------------------------------------
        !  calculate <q1| E_n |q2>/<q1|R|q2>
        !--------------------------------------------------
        use Constants,only: itx
        use Globals, only: BS, mix
        integer :: n,iphi,it
        integer :: ifg1,ifg2,ifg3,ifg4,m1,m2,m3,m4
        complex(r64) :: Eccentri,pEccentri,Eccentri_1B,pEccentri_1B,Eccentri_2B,pEccentri_2B
        real(r64) :: e_1B,e_2B

        Eccentri_1B = 0.d0
        pEccentri_1B = 0.d0
        Eccentri_2B = 0.d0
        pEccentri_2B = 0.d0
        do ifg1 = 1,2
        do m1 = 1,BS%HO_sph%iasp(1,ifg1)
            do ifg2 = 1,2
            do m2 = 1,BS%HO_sph%iasp(1,ifg2)
                ! 1 Body
                if(ifg1==ifg2) then 
                    call eccentricity_matrix_element_one_body(ifg1,m1,m2,2,e_1B)
                    Eccentri_1B = Eccentri_1B + e_1B*mix%rho_mm(m2,m1,ifg1,iphi,it)
                    pEccentri_1B = pEccentri_1B + e_1B*mix%rho_mm(m2,m1,ifg1,iphi,it)
                end if 
                do ifg3 = 1,2
                do m3 = 1,BS%HO_sph%iasp(1,ifg3)
                    do ifg4 = 1,2
                    do m4 = 1,BS%HO_sph%iasp(1,ifg4)
                        ! 2 Body
                        if(ifg1==ifg2 .and. ifg2==ifg3 .and. ifg3==ifg4) then 
                            call eccentricity_matrix_element_two_body(ifg1,m1,m2,m3,m4,n,e_2B)
                            Eccentri_2B = Eccentri_2B+ 1.d0/4.d0*e_2B*( &
                                          mix%rho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%rho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                        - mix%rho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%rho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                        + mix%kappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%kappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it))
                            pEccentri_2B = pEccentri_2B + 1.d0/4.d0*e_2B*( &
                                          mix%prho_mm(m4,m1,indexfg(ifg4,ifg1),iphi,it)*mix%prho_mm(m3,m2,indexfg(ifg3,ifg2),iphi,it) &
                                        - mix%prho_mm(m3,m1,indexfg(ifg3,ifg1),iphi,it)*mix%prho_mm(m4,m2,indexfg(ifg4,ifg2),iphi,it) &
                                        + mix%pkappa01c_mm(m1,m2,indexfg(ifg1,ifg2),iphi,it)*mix%pkappa10_mm(m4,m3,indexfg(ifg4,ifg3),iphi,it))
                        end if 
                    end do 
                    end do               
                end do 
                end do 
            end do 
            end do 
        end do 
        end do 
        Eccentri = Eccentri_1B + Eccentri_2B
        pEccentri = pEccentri_1B + pEccentri_2B
    end subroutine

    subroutine eccentricity_matrix_element_one_body(ifg,m1,m2,n,EME1B)
        !-----------------------------------------------------------
        !   e^{n}_{m1 m2} = \sum_m f^n_{m1 m} f^{-n}_{m m2} 
        !-----------------------------------------------------------
        use Globals, only: BS
        integer :: ifg,m1,m2,n
        real(r64) :: EME1B
        integer :: ndsp,i0sp,nr1,nl1,nj1,nm1,nr2,nl2,nj2,nm2,m,nr,nl,nj,nm
        EME1B =0.d0
        ndsp = BS%HO_sph%idsp(1,ifg)
        i0sp = BS%HO_sph%iasp(1,ifg)

        nr1 = BS%HO_sph%nljm(i0sp+m1,1)
        nl1 = BS%HO_sph%nljm(i0sp+m1,2)
        nj1 = BS%HO_sph%nljm(i0sp+m1,3)
        nm1 = BS%HO_sph%nljm(i0sp+m1,4)

        nr2 = BS%HO_sph%nljm(i0sp+m2,1)
        nl2 = BS%HO_sph%nljm(i0sp+m2,2)
        nj2 = BS%HO_sph%nljm(i0sp+m2,3)
        nm2 = BS%HO_sph%nljm(i0sp+m2,4)

        do m = 1, ndsp
            nr= BS%HO_sph%nljm(i0sp+m,1) ! n_r
            nl= BS%HO_sph%nljm(i0sp+m,2) ! l
            nj= BS%HO_sph%nljm(i0sp+m,3) ! j +1/2
            nm= BS%HO_sph%nljm(i0sp+m,4) ! m_j + 1/2
            EME1B = EME1B + f_n(nr1,nl1,nj1,nm1,nr,nl,nj,nm,n)*f_n(nr,nl,nj,nm,nr2,nl2,nj2,nm2,-n)
        end do
    end subroutine

    subroutine eccentricity_matrix_element_two_body(ifg,m1,m2,m3,m4,n,EME2B)
        !-------------------------------------------------------------------------------
        !   e^{n}_{m1 m2 m3 m4} = 2(f^n_{m1 m3}f^{-n}_{m2 m4} - f^{-n}_{m1 m4}f^n_{m2 m3} 
        !-------------------------------------------------------------------------------
        use Globals, only: BS
        integer :: ifg,m1,m2,m3,m4,n
        real(r64) :: EME2B
        integer :: ndsp,i0sp,nr1,nl1,nj1,nm1,nr2,nl2,nj2,nm2,nr3,nl3,nj3,nm3,nr4,nl4,nj4,nm4
        ndsp = BS%HO_sph%idsp(1,ifg)
        i0sp = BS%HO_sph%iasp(1,ifg)

        nr1 = BS%HO_sph%nljm(i0sp+m1,1)
        nl1 = BS%HO_sph%nljm(i0sp+m1,2)
        nj1 = BS%HO_sph%nljm(i0sp+m1,3)
        nm1 = BS%HO_sph%nljm(i0sp+m1,4)

        nr2 = BS%HO_sph%nljm(i0sp+m2,1)
        nl2 = BS%HO_sph%nljm(i0sp+m2,2)
        nj2 = BS%HO_sph%nljm(i0sp+m2,3)
        nm2 = BS%HO_sph%nljm(i0sp+m2,4)

        nr3 = BS%HO_sph%nljm(i0sp+m3,1)
        nl3 = BS%HO_sph%nljm(i0sp+m3,2)
        nj3 = BS%HO_sph%nljm(i0sp+m3,3)
        nm3 = BS%HO_sph%nljm(i0sp+m3,4)

        nr4 = BS%HO_sph%nljm(i0sp+m4,1)
        nl4 = BS%HO_sph%nljm(i0sp+m4,2)
        nj4 = BS%HO_sph%nljm(i0sp+m4,3)
        nm4 = BS%HO_sph%nljm(i0sp+m4,4)

        EME2B = 2.d0*(f_n(nr1,nl1,nj1,nm1,nr3,nl3,nj3,nm3,n)*f_n(nr2,nl2,nj2,nm2,nr4,nl4,nj4,nm4,-n)&
                     -f_n(nr1,nl1,nj1,nm1,nr4,nl4,nj4,nm4,-n)*f_n(nr2,nl2,nj2,nm2,nr3,nl3,nj3,nm3,n))
        
    end subroutine

    double precision function f_n(n1,l1,j1,m1,n2,l2,j2,m2,n)
        !-------------------------------------------
        !   calculate <n1 l1 j1 m1| F_n | n1 l2 j2 m2>
        !   where F_n =  r^{|n|} * Y_{|n|n}
        !--------------------------------------------
        use EM, only: multipole_matrix_elements
        integer :: n1,l1,j1,m1,n2,l2,j2,m2,n
        real(r64) :: mpme
        f_n = 0.d0
        call multipole_matrix_elements(n1,l1,j1,m1,abs(n),n,n2,l2,j2,m2,mpme)
        f_n =  mpme
    end function

    integer function indexfg(ifg1,ifg2)
        integer, intent(in) :: ifg1, ifg2
        indexfg = ifg1
        ! if(ifg1==1.and.ifg2==1) then
        !     indexfg = 1 ! ++
        ! else if(ifg1==2.and.ifg2==2) then
        !     indexfg = 2 ! --
        ! else if(ifg1==1.and.ifg2==2) then
        !     indexfg = 3 ! +-
        ! else if(ifg1==2.and.ifg2==1) then
        !     indexfg = 4 ! -+
        ! else 
        !     stop 'wrong ifg1 and ifg2'
        ! end if 
    end function
end Module Eccentricity