Module Eccentricity
    use Constants,only: r64
    implicit none
    
contains
    
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

        nr3 = BS%HO_sph%nljm(i0sp+m1,1)
        nl3 = BS%HO_sph%nljm(i0sp+m1,2)
        nj3 = BS%HO_sph%nljm(i0sp+m1,3)
        nm3 = BS%HO_sph%nljm(i0sp+m1,4)

        nr4 = BS%HO_sph%nljm(i0sp+m2,1)
        nl4 = BS%HO_sph%nljm(i0sp+m2,2)
        nj4 = BS%HO_sph%nljm(i0sp+m2,3)
        nm4 = BS%HO_sph%nljm(i0sp+m2,4)

        EME2B = 2.d0*(f_n(nr1,nl1,nj1,nm1,nr3,nl3,nj3,nm3,n)*f_n(nr2,nl2,nj2,nm2,nr4,nl4,nj4,nm4,-n)&
                     -f_n(nr1,nl1,nj1,nm1,nr4,nl4,nj4,nm4,-n)*f_n(nr2,nl2,nj2,nm2,nr3,nl3,nj3,nm3,n))
        
    end subroutine

    double precision function f_n(n1,l1,j1,m1,n2,l2,j2,m2,n)
        !-------------------------------------------
        !   calculate <n1 l1 j1 m1| F_n | n1 l2 j2 m2>
        !   where F_n = 1/c_n * r^{|n|} * Y_{|n|n}
        !   while n > 0:  c_n = (-1)^n / (2^n n!) \sqrt((2n+1)!/(4pi))
        !   while n < 0 : c_n = 1 / (2^n n!) \sqrt((2n+1)!/(4pi))
        !--------------------------------------------
        use Constants, only: pi
        use Globals, only: gfv
        use EM, only: multipole_matrix_elements
        integer :: n1,l1,j1,m1,n2,l2,j2,m2,n
        real(r64) :: c_n, mpme
        f_n = 0.d0
        if (n>0) then
            c_n = (-1)**n / (2**n * gfv%fak(n)) * sqrt(gfv%fak(2*n+1)/(4*pi)) 
        else 
            c_n = 1.d0 / (2**n * gfv%fak(n)) * sqrt(gfv%fak(2*n+1)/(4*pi))
        end if 
        c_n = 1.d0
        call multipole_matrix_elements(n1,l1,j1,m1,abs(n),n,n2,l2,j2,m2,mpme)
        f_n= 1.d0/c_n * mpme
    end function
end Module Eccentricity