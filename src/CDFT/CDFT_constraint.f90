!==============================================================================!
! MODULE Constraint                                                            !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
Module Constraints
use Constants, only: r64,zero,one
use Globals, only: constraint
implicit none
contains

subroutine select_the_constraint_deformation(factor_index,option)
    use Globals, only: constraint, Factors
    integer :: factor_index, option
    if(option==0) then 
        ! deformation from b23.dat file
        write(*,*) 'deformation from b23.dat.'
    else 
        ! deformation from factor file
        constraint%length = 1
        constraint%betac(1) = Factors%beta2(factor_index)
        constraint%bet3c(1) = 0.d0
        write(*,*) 'deformation from factor file.'
    end if

end subroutine

subroutine select_the_constraint_deformation_space(factor_index1,factor_index2,option)
    use Globals, only: constraint,constraint1,constraint2,Factors
    integer :: factor_index1,factor_index2,option
    if(.not. allocated(constraint1%betac)) allocate(constraint1%betac(size(constraint%betac)))
    if(.not. allocated(constraint1%bet3c)) allocate(constraint1%bet3c(size(constraint%bet3c)))
    if(.not. allocated(constraint1%clam2)) allocate(constraint1%clam2(size(constraint%clam2)))
    if(.not. allocated(constraint2%betac)) allocate(constraint2%betac(size(constraint%betac)))
    if(.not. allocated(constraint2%bet3c)) allocate(constraint2%bet3c(size(constraint%bet3c)))
    if(.not. allocated(constraint2%clam2)) allocate(constraint2%clam2(size(constraint%clam2)))
    if(option==0) then ! from b23.dat
        constraint1%length = constraint%length
        constraint1%betac = constraint%betac
        constraint1%bet3c = constraint%bet3c
        constraint2%length = constraint%length
        constraint2%betac = constraint%betac
        constraint2%bet3c = constraint%bet3c
    else ! from factor file
        constraint1%length = 1
        constraint1%betac(1) = Factors%beta2(factor_index1)
        constraint1%bet3c(1) = 0.d0
        constraint2%length = 1
        constraint2%betac(1) = Factors%beta2(factor_index2)
        constraint2%bet3c(1) = 0.d0
    end if
end subroutine

subroutine set_constraint_parameters
    !----------------------------------------------------------------
    ! `constraint%c1x constraint%c2x constraint%c3x` will be updated 
    !  after each iteration.
    !----------------------------------------------------------------
    use Constants, only: radius_r0,third,pi,nghl,ngh,two
    use Globals, only: input_par,nucleus_attributes,gauss
    integer :: ih,il,ii
    real(r64) :: R00,tmp_fac,tmp_c1,tmp_c2,tmp_c3
    ! set input parameter 
    constraint%icstr = input_par%constraint_icstr ! ! quadratic constraint (no 0; beta2 1; b2+b3 2)
    constraint%cspr  = input_par%constraint_cspr * nucleus_attributes%mass_number
    constraint%cmax  = input_par%constraint_cmax * nucleus_attributes%mass_number

    ! initial constraining field
    constraint%c1xold = zero
    constraint%c2xold = zero
    constraint%c3xold = zero
    constraint%c1x = zero
    constraint%c2x = zero
    constraint%c3x = zero

    R00 = radius_r0 * nucleus_attributes%mass_number**third
    tmp_fac = 4*pi/3/nucleus_attributes%mass_number
    tmp_c1 = tmp_fac*dsqrt(3/(4*pi)) /R00
    tmp_c2 = tmp_fac*dsqrt(5/(16*pi))/(R00**2)
    tmp_c3 = tmp_fac*dsqrt(7/(16*pi))/(R00**3)
    do ii =1,nghl
        il = 1 + (ii-1)/ngh
        ih = ii - (il-1)*ngh
        constraint%vc(ii,1) = tmp_c1 * gauss%zb(ih)
        constraint%vc(ii,2) = tmp_c2 * (two*gauss%zb(ih)**2 - gauss%rb(il)**2) ! Quardrupole: (4pi/(3A))* \sqrt(5/16pi)*(2z^2-x^2-y^2)/(R_0^2)
        constraint%vc(ii,3) = tmp_c3 * (two*gauss%zb(ih)**3 - 3.d0*gauss%zb(ih)*gauss%rb(il)**2) ! Octupole: (4pi/(3A))* \sqrt(7/16pi)*(2z^3-3z(x^2-y^2)^2)(R_0^3)
    enddo
end subroutine set_constraint_parameters

subroutine calculate_constraint_potential_coefficients(ifPrint)
    !----------------------------------------------------------------
    !    calculation of the constrainting potential for axial case
    ! 
    !   Note: 
    !   1) Constrain the dipole deformation to zero to keep the 
    !       center of mass at the origin.
    !----------------------------------------------------------------
    logical,intent(in),optional :: ifPrint
    real(r64) :: d1,dc1x,d2,dc2x,d3,dc3x
    if (constraint%icstr.eq.0) then
        d1   = constraint%calq1-zero
        dc1x = d1/(one/constraint%cspr+abs(d1)/constraint%cmax)
        constraint%c1x  = constraint%c1x + dc1x
    else if (constraint%icstr.eq.1) then
        d1   = constraint%calq1-zero
        d2   = constraint%calq2 - constraint%betac(constraint%index)
        dc1x = d1/(one/constraint%cspr+abs(d1)/constraint%cmax)
        dc2x = d2/(one/constraint%cspr+abs(d2)/constraint%cmax)
        constraint%c1x  = constraint%c1x + dc1x
        constraint%c2x  = constraint%c2x + dc2x
    elseif(constraint%icstr.eq.2) then
        d1   = constraint%calq1-zero
        d2   = constraint%calq2 - constraint%betac(constraint%index)
        d3   = constraint%calq3 - constraint%bet3c(constraint%index)
        dc1x = d1/(one/constraint%cspr+abs(d1)/constraint%cmax)
        dc2x = d2/(one/constraint%cspr+abs(d2)/constraint%cmax)
        dc3x = d3/(one/constraint%cspr+abs(d3)/constraint%cmax)
        constraint%c1x  = constraint%c1x + dc1x
        constraint%c2x  = constraint%c2x + dc2x
        constraint%c3x  = constraint%c3x + dc3x
    else if(constraint%icstr.eq.3) then
        d1   = constraint%calq1-zero
        d3   = constraint%calq3 - constraint%bet3c(constraint%index)
        dc1x = d1/(one/constraint%cspr+abs(d1)/constraint%cmax)
        dc3x = d3/(one/constraint%cspr+abs(d3)/constraint%cmax)
        constraint%c1x  = constraint%c1x + dc1x
        constraint%c3x  = constraint%c3x + dc3x
    endif
    
    if(ifPrint) call print_cstrpot
    contains
    subroutine print_cstrpot
        use Globals,only: outputfile
        write(outputfile%u_outputf,*) '*************************BEGIN print_cstrpot ********************'
        write(outputfile%u_outputf,*) 'c1x', constraint%c1x
        write(outputfile%u_outputf,*) 'c2x', constraint%c2x
        write(outputfile%u_outputf,*) 'c3x', constraint%c3x
        write(outputfile%u_outputf,*) '*************************END print_cstrpot **********************'
    end subroutine print_cstrpot
end subroutine calculate_constraint_potential_coefficients

END Module Constraints