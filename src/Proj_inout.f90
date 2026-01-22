!==============================================================================!
! MODULE Proj_Inout                                                            !
!                                                                              !
! This module contains functions and routines for reading and writing files.   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Proj_Inout
use Globals, only: outputfile
use Constants, only: i16,r64,u_start,pi,ngl,OUTPUT_PATH, Jmax_max
use CDFT_Inout, only: file_path_para,set_output_filename, int2str, adjust_left
implicit none
integer, private :: u_pko = u_start + 11

contains

subroutine read_Proj_configuration(ifPrint)
    use Globals, only: input_par,pko_option
    integer :: i,is
    logical,intent(in),optional :: ifPrint
    character(len=*), parameter ::  format1 = "(10x,2f10.4)", &
                                    format2 = "(10x,i5)", &
                                    format3 = "(10x,2i5)" 

    open(u_pko, file=file_path_para, status='old')
    ! skip CDFT parameters
    do i = 1, 29
        read(u_pko,'(A)', iostat=is) 
        if (is /= 0) then
            print *, "Error reading file"
            stop
        end if
    end do

    read(u_pko, format2) input_par%ProjectionType
    read(u_pko, format2) input_par%AMPType
    read(u_pko, format2) input_par%PNPType
    read(u_pko, format2) input_par%Kernel_Symmetry
    read(u_pko, format3) input_par%q1_start, input_par%q1_end
    read(u_pko, format3) input_par%q2_start, input_par%q2_end
    read(u_pko, format2) input_par%Jmax
    read(u_pko, format2) input_par%icm
    read(u_pko, format2) input_par%nphi
    read(u_pko, format2) input_par%Euler_Symmetry
    read(u_pko, format2) input_par%nalpha
    read(u_pko, format2) input_par%nbeta
    read(u_pko, format2) input_par%ngamma
    read(u_pko, format2) input_par%TDType
    read(u_pko, format2) input_par%lambda_max
    call set_pko_parameters
    if(ifPrint) call printParameters
    contains
    subroutine set_pko_parameters
        use Globals, only: input_par,pko_option,gcm_space,constraint,TDs
        ! set projection type ( 0 : no 1: RMF+AMP 2: only AMP)
        pko_option%ProjectionType = input_par%ProjectionType
        if(pko_option%ProjectionType > 2 .or. pko_option%ProjectionType < 0) stop 'ProjectionType wrong!'

        ! set AMP type ((0) no (1) 1DAMP (2) 3DAMP)
        pko_option%AMPtype = input_par%AMPType
        if(input_par%AMPType.ne.0 .and. input_par%AMPType.ne.1 .and. input_par%AMPType.ne.2) stop 'AMPType wrong!'
        if( input_par%AMPType==2) stop '(3DAMP) Not yet verified !'

        ! set PNP type (0: no PNP; 1: PNP)
        pko_option%PNPtype = input_par%PNPType
        if(input_par%PNPType.ne.0 .and. input_par%PNPType.ne.1 ) stop 'PNPType wrong!'

        ! Kernel Symmetry  (0: All ; 1: Triangular Matrix ; 2: Diagonal elements only)
        pko_option%Kernel_Symmetry = input_par%Kernel_Symmetry
        if(pko_option%Kernel_Symmetry > 2 .or. pko_option%Kernel_Symmetry < 0) stop 'Kernel_Symmetry wrong!'

        ! Euler angles Symmetry (0: no, 1: Axially, 2: D2)
        pko_option%Euler_Symmetry = input_par%Euler_Symmetry
        if(pko_option%Euler_Symmetry .ne. 0 .and. pko_option%Euler_Symmetry .ne. 1 .and. pko_option%Euler_Symmetry .ne. 2) stop 'Euler_Symmetry wrong!'
        if(input_par%AMPType == 1  .and. (pko_option%Euler_Symmetry .ne. 0 .and. pko_option%Euler_Symmetry .ne. 1)) stop '1DAMP: Euler_Symmetry should be 0 or 1!'
        if(input_par%AMPType == 2  .and. (pko_option%Euler_Symmetry .ne. 2)) stop '3DAMP: Euler_Symmetry should 2 !'

        ! center of mass correction (1: average  ; 2: HO approximation)
        pko_option%icm = input_par%icm
        ! Set Norm overlap calculation method
        pko_option%ihf = 3

        ! set GCM space
        gcm_space%Jmin = 0
        gcm_space%Jmax = input_par%Jmax
        if(gcm_space%Jmax > Jmax_max) stop 'Jmax_max too small!'
        gcm_space%Jstep = 1 
        gcm_space%q1_start = input_par%q1_start
        gcm_space%q1_end = input_par%q1_end
        if(gcm_space%q1_end == -1) gcm_space%q1_end = constraint%length
        gcm_space%q2_start = input_par%q2_start
        gcm_space%q2_end = input_par%q2_end
        if(gcm_space%q2_end == -1) gcm_space%q2_end = constraint%length
        if(gcm_space%q1_start < 1) stop 'q1_start less than 1 !'
        if(gcm_space%q1_end < gcm_space%q1_start) stop 'q1_end less than q1_start !'
        if(gcm_space%q1_end > constraint%length) stop 'q1_end greater than max length.'
        if(gcm_space%q2_start < 1) stop 'q2_start less than 1 !'
        if(gcm_space%q2_end < gcm_space%q2_start) stop 'q2_end less than q2_start !'
        if(gcm_space%q2_end > constraint%length) stop 'q2_end greater than max length.'

        ! set TD type
        pko_option%TDType = input_par%TDType
        ! max lambda of 1B transition density
        TDs%lambda_max = input_par%lambda_max
    end subroutine
    subroutine printParameters
        use Globals, only: input_par,pko_option,gcm_space,TDs
        character(len=5) :: AMP_char, PNP_char
        if(pko_option%AMPtype==0) then 
            AMP_char = 'noAMP'
        else if (pko_option%AMPtype==1) then
            AMP_char = '1DAMP'
        else if (pko_option%AMPtype==2) then
            AMP_char = '3DAMP'
        end if 
        if(pko_option%PNPtype==0) then 
            PNP_char = 'noPNP'
        else if (pko_option%PNPtype==1) then
            PNP_char = 'PNP'
        end if 

        write(*,'(5x,A)') AMP_char//'  +  '//PNP_char//'+  '//'PP :'
        if(pko_option%AMPtype /= 0) then
            write(*,"(5x,a,':   ',3(i2,a))") adjust_left('Number euler angles',35),input_par%nalpha,' (nalpha),  ',input_par%nbeta,' (nbeta),  ',input_par%ngamma,' (ngamma)'
            if(input_par%Euler_Symmetry==0) then
                write(*,"(5x,a,':   ',a)") adjust_left('Symmetry of euler angles',35),'no'
            else if(input_par%Euler_Symmetry==1) then
                write(*,"(5x,a,':   ',a)") adjust_left('Symmetry of euler angles',35),'Axially'
            else if(input_par%Euler_Symmetry==2) then
                write(*,"(5x,a,':   ',a)") adjust_left('Symmetry of euler angles',35),'D2'
            end if 
        end if 
        if(pko_option%PNPtype /= 0) then
            write(*,"(5x,a,':   ',i2,a)") adjust_left('Number gauge angles',35),input_par%nphi,' (nphi)'
        end if

        if(pko_option%Kernel_Symmetry==0) then 
            write(*,"(5x,a,':   ',a)") adjust_left('Kernels',35),'All kernels'
        else if(pko_option%Kernel_Symmetry==1) then
            write(*,"(5x,a,':   ',a)") adjust_left('Kernels',35),'Upper triangular kernels'
        else if(pko_option%Kernel_Symmetry==2) then
            write(*,"(5x,a,':   ',a)") adjust_left('Kernels',35),'Diagonal kernels'
        end if 
        write(*,"(5x,a,':   ',a,i3,a,i3,a)") adjust_left('Quadratic constraint q1 range',35), '[',gcm_space%q1_start,',',gcm_space%q1_end,' ]'
        write(*,"(5x,a,':   ',a,i3,a,i3,a)") adjust_left('Quadratic constraint q2 range',35), '[',gcm_space%q2_start,',',gcm_space%q2_end,' ]'

        write(*,"(5x,a,':   ',i3)") adjust_left('Maximal J value',35), gcm_space%Jmax
        write(*,"(5x,a,':   ',i3)") adjust_left('Maximal lambda value(1BTD)',35),TDs%lambda_max

        if(pko_option%ihf == 1) then 
            write(*,"(5x,a,':   ',a)") adjust_left('Norm overlap formula',35),'sqrt(det(D) det(R))'
        else if(pko_option%ihf == 2) then
            write(*,"(5x,a,':   ',a)") adjust_left('Norm overlap formula',35),'Robledo (2009) formula'
        else if(pko_option%ihf == 3) then
            write(*,"(5x,a,':   ',a)") adjust_left('Norm overlap formula',35),'Bertsch & Robledo (2011) formula'
        end if 
        write(*,"(a)") '=========================================================================================='
    end subroutine
end subroutine

subroutine read_wavefuntion_files(q1,q2)
    use Constants, only: nb_max 
    use Globals, only: wf1,wf2,constraint
    integer :: q1,q2,nb,it,temp,ib
    integer, dimension(nb_max,2) :: kd
    call set_output_filename(constraint%betac(q1),constraint%bet3c(q1))
    open(outputfile%u_outputwf,file=outputfile%outputwf,form='unformatted',status='unknown')      
    read(outputfile%u_outputwf) !dirac%ka
    read(outputfile%u_outputwf) kd
    read(outputfile%u_outputwf) wf1%v2
    read(outputfile%u_outputwf) wf1%v2d
    read(outputfile%u_outputwf) wf1%skk
    read(outputfile%u_outputwf) wf1%fg  
    read(outputfile%u_outputwf) nb
    read(outputfile%u_outputwf) wf1%ecm 
    read(outputfile%u_outputwf) !pairing%ibk
    close(outputfile%u_outputwf)
    do it=1,2
        temp = 0
        do ib=1,nb
            temp = temp + kd(ib,it)
        enddo
        wf1%nk(it) = temp
    enddo 
    call set_output_filename(constraint%betac(q2),constraint%bet3c(q2))
    open(outputfile%u_outputwf,file=outputfile%outputwf,form='unformatted',status='unknown')      
    read(outputfile%u_outputwf) !dirac%ka
    read(outputfile%u_outputwf) kd
    read(outputfile%u_outputwf) wf2%v2
    read(outputfile%u_outputwf) wf2%v2d
    read(outputfile%u_outputwf) wf2%skk
    read(outputfile%u_outputwf) wf2%fg  
    read(outputfile%u_outputwf) nb
    read(outputfile%u_outputwf) wf2%ecm 
    read(outputfile%u_outputwf) !pairing%ibk
    close(outputfile%u_outputwf)
    do it=1,2
        temp = 0
        do ib=1,nb
            temp = temp + kd(ib,it)
        enddo
        wf2%nk(it) = temp
    enddo
end subroutine

subroutine write_pko_output(q1,q2)
    use Globals, only: gcm_space,pko_option
    integer,intent(in) :: q1,q2
    call set_pko_output_filename(q1,q2,pko_option%AMPType)
    call write_kernels
    if (pko_option%TDType==1) then
        call write_reduced_1B_transition_density_matrix_elements(q1,q2)
        if(q1== gcm_space%q1_start .and. q2==gcm_space%q2_start) then 
            call write_reduced_1B_multipole_matrix_elements
        end if 
    end if 
end subroutine

subroutine set_pko_output_filename(q1,q2,AMPType)
    use Globals, only:constraint,BS,projection_mesh,nucleus_attributes
    integer :: q1,q2,AMPType,A
    real(r64) :: beta2_1, beta3_1,beta2_2,beta3_2,abs2c1,abs3c1,abs2c2,abs3c2
    character :: signb21,signb31,signb22,signb32
    integer(i16), dimension(6) :: name1,name2
    integer(i16) :: AMP,name_nf1,name_nf2,nphi_1,nphi_2,nbeta_1,nbeta_2
    A = nucleus_attributes%mass_number_int
    beta2_1 = constraint%betac(q1)
    beta3_1 = constraint%bet3c(q1)
    beta2_2 = constraint%betac(q2)
    beta3_2 = constraint%bet3c(q2)
    if(beta2_1 >= 0.d0) signb21 = '+'
    if(beta2_1 < 0.d0)  signb21 = '-'
    if(beta3_1 >= 0.d0) signb31 = '+'
    if(beta3_1 < 0.d0)  signb31 = '-'
    if(beta2_2 >= 0.d0) signb22 = '+'
    if(beta2_2 < 0.d0)  signb22 = '-'
    if(beta3_2 >= 0.d0) signb32 = '+'
    if(beta3_2 < 0.d0)  signb32 = '-'
    !------
    if(AMPType==0) AMP = 0 + 48
    if(AMPType==1) AMP = 1 + 48
    if(AMPType==3) AMP = 3 + 48
    !-------
    abs2c1 = abs(beta2_1)
    abs3c1 = abs(beta3_1)
    name1(1) = abs2c1 + 48 !In ASCII, character '0' start from 48. 
    name1(2) = mod(abs2c1*10,10.d0)+48
    name1(3) = mod(abs2c1*100,10.d0)+48
    name1(4) = abs3c1+48
    name1(5) = mod(abs3c1*10,10.d0)+48
    name1(6) = mod(abs3c1*100,10.d0)+48
    !------
    abs2c2 = abs(beta2_2)
    abs3c2 = abs(beta3_2)
    name2(1) = abs2c2 + 48 !In ASCII, character '0' start from 48. 
    name2(2) = mod(abs2c2*10,10.d0)+48
    name2(3) = mod(abs2c2*100,10.d0)+48
    name2(4) = abs3c2+48
    name2(5) = mod(abs3c2*10,10.d0)+48
    name2(6) = mod(abs3c2*100,10.d0)+48
    !-----
    name_nf1 = mod(BS%HO_sph%n0f/10,10) + 48
    name_nf2 = mod(BS%HO_sph%n0f,10) + 48
    !-----
    nphi_1 = mod(projection_mesh%nphi(1)/10,10) + 48
    nphi_2 = mod(projection_mesh%nphi(1),10) + 48
    nbeta_1 = mod(projection_mesh%nbeta/10,10) + 48
    nbeta_2 = mod(projection_mesh%nbeta,10) + 48
    outputfile%outputelem = OUTPUT_PATH//'kern.'//char(AMP)//'D' &
                        //'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2) &
                        //'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.elem'
    outputfile%outputTDME1B = OUTPUT_PATH//'TD1B.'//char(AMP)//'D' &
                        //'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2) &
                        //'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) &
                        //'_'//signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6))//'.dens'
    outputfile%outputTDME1B_c = OUTPUT_PATH//'TD1B.'//char(AMP)//'D' &
                        //'_eMax'//char(name_nf1)//char(name_nf2) &
                        //'.'//char(nphi_1)//char(nphi_2) &
                        //'.'//char(nbeta_1)//char(nbeta_2) &
                        //signb22//char(name2(1))//char(name2(2))//char(name2(3)) &
                        //signb32//char(name2(4))//char(name2(5))//char(name2(6)) &
                        //'_'//signb21//char(name1(1))//char(name1(2))//char(name1(3)) &
                        //signb31//char(name1(4))//char(name1(5))//char(name1(6)) //'.dens'
    outputfile%outputEMelem = OUTPUT_PATH//'EM'//'_A'//int2str(A) &
                        //'_eMax'//char(name_nf1)//char(name_nf2)//'.elem'       
end subroutine

subroutine write_kernels
    use Globals, only: gcm_space,kernels
    integer :: J,K1,K2,parity
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(len=*), parameter ::  format1 = "(3i5,4x,a)", &
                                    format2 = "(4e15.8)"
    open(outputfile%u_outputelem ,form='formatted',file=outputfile%outputelem)
        do J = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
            do K1 = -0,0
                do K2 = -0,0
                    ! In the axially symmetric case, the kernel is non-zero only when
                    ! the parity satisfies  Pi  = (-1)^J for  N_KK, H_KK, X_KK and E0_KK
                    ! the parity satisfies Pi_i = (-1)^J_i for  Q2_KK_12
                    if ((-1)**J == 1) then
                        parity = 1 ! +
                    else
                        parity = 2 ! -
                    end if
                    write(outputfile%u_outputelem,format1)  J,K1,K2,ParityChar(parity)
                    write(outputfile%u_outputelem,format2)  kernels%N_KK(J,K1,K2,parity), kernels%H_KK(J,K1,K2,parity)/kernels%N_KK(J,K1,K2,parity)
                    write(outputfile%u_outputelem,format2)  kernels%X_KK(J,K1,K2,1,parity)/kernels%N_KK(J,K1,K2,parity),kernels%X_KK(J,K1,K2,2,parity)/kernels%N_KK(J,K1,K2,parity)
                    ! ! proton part
                    write(outputfile%u_outputelem,format2)  kernels%Q2_KK_12(J,K1,K2,2,parity),kernels%Q2_KK_21(J,K1,K2,2,parity)
                    write(outputfile%u_outputelem,format2)  kernels%E0_KK(J,K1,K2,2,parity), kernels%E0_KK(J,K1,K2,2,parity)
                    ! ! neutron part
                    write(outputfile%u_outputelem,format2)  kernels%Q2_KK_12(J,K1,K2,1,parity),kernels%Q2_KK_21(J,K1,K2,1,parity)
                    write(outputfile%u_outputelem,format2)  kernels%E0_KK(J,K1,K2,1,parity), kernels%E0_KK(J,K1,K2,1,parity)
                end do
            end do
        end do 
    close(outputfile%u_outputelem)
end subroutine

subroutine write_reduced_1B_transition_density_matrix_elements(q1,q2)
    use Globals, only: outputfile,gcm_space,pko_option,TDs
    integer, intent(in) :: q1,q2
    integer :: J,Ji,Jf,lambda,Ki_start,Ki_end,Kf_start,Kf_end,Kf,Ki,ifg,a,b,Parity_i,Parity_f,iPi,iPf
    character(1), dimension(2) :: ParityChar = ['+', '-']
    character(1) :: Parity_f_c,Parity_i_c
    ! q1-q2
    open(outputfile%u_outputTDME1B ,form='formatted',file=outputfile%outputTDME1B)
    write(outputfile%u_outputTDME1B,*) "Pf Pi  Jf  Ji  l   Kf  Ki  ifg a   b    neutron            proton"
    ! q2-q1
    if (pko_option%Kernel_Symmetry == 1 .and. q1/=q2) then
        open(outputfile%u_outputTDME1B_c ,form='formatted',file=outputfile%outputTDME1B_c)
        write(outputfile%u_outputTDME1B_c,*) "Pf Pi  Jf  Ji  l   Kf  Ki  ifg a   b    neutron            proton"
    end if 
    do Ji = gcm_space%Jmin, gcm_space%Jmax, gcm_space%Jstep
        do Jf = Ji, Ji+TDs%lambda_max
            if(pko_option%AMPtype==0 .or. pko_option%AMPtype==1) then
                Ki_start = 0
                Ki_end = 0
                Kf_start = 0
                Kf_end = 0
            else
                Ki_start = -Ji
                Ki_end = Ji
                Kf_start = -Jf
                Kf_end = Jf
            end if
            do lambda = abs(Jf-Ji), min(TDs%lambda_max,Jf+Ji)
                do Kf = Kf_start, Kf_end
                    do Ki = Ki_start, Ki_end 
                        do ifg = 1, 2 
                            do a = 1, TDs%nlj_length(ifg)
                                do b = 1, TDs%nlj_length(ifg)
                                    Parity_i = (-1)**Ji ! In the axially symmetric case, the kernel is non-zero only when Parity_i * (-1)^J_i = 1
                                    iPi = (3-Parity_i)/2 ! +1: 1, -1: 2
                                    Parity_i_c = ParityChar(iPi)
                                    do iPf = 1,2
                                        Parity_f = (-1)**(iPf+1) ! 1: +1 , 2: -1
                                        Parity_f_c = ParityChar(iPf)
                                        write(outputfile%u_outputTDME1B,"(1x,a1,2x,a1,8i4,2x,2f18.14)") Parity_f_c,Parity_i_c,Jf,Ji,lambda,Kf,Ki,ifg,a,b, &
                                            Real(TDs%reduced_TDME1B(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,1)),&
                                            Real(TDs%reduced_TDME1B(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,2))
                                        if (pko_option%Kernel_Symmetry == 1 .and. q1/=q2) then
                                            write(outputfile%u_outputTDME1B_c,"(1x,a1,2x,a1,8i4,2x,2f18.14)")Parity_f_c,Parity_i_c,Jf,Ji,lambda,Kf,Ki,ifg,a,b, &
                                                Real(TDs%reduced_TDME1B_c(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,1)),&
                                                Real(TDs%reduced_TDME1B_c(Jf,Kf,iPf,lambda,Ji,Ki,iPi,ifg,a,b,2))
                                        end if 
                                    end do
                                end do 
                            end do 
                        end do
                    end do 
                end do
            end do
        end do 
    end do
    close(outputfile%u_outputTDME1B)
    if (pko_option%Kernel_Symmetry == 1 .and. q1/=q2) then 
        close(outputfile%u_outputTDME1B_c)
    end if 
end subroutine

subroutine  write_reduced_1B_multipole_matrix_elements
    use Globals, only: outputfile,BS,TDs,gcm_space
    use EM, only: reduced_multipole_matrix_elements,reduced_monopole_matrix_elements,rl_nl
    integer :: ifg,lambda_start,lambda_end,lambda,a,b,i0sp,a_index,nra,nla,nja,b_index,nrb,nlb,njb
    real(r64),allocatable,dimension(:) :: Ql_ab
    real(r64) :: monopole_ab
    character(len=200) :: header, temp
    open(outputfile%u_outputEMelem ,form='formatted',file=outputfile%outputEMelem)
    ! store the nlj of a/b
    write(outputfile%u_outputEMelem,*) "--------------------------------------------------------"
    write(outputfile%u_outputEMelem,*) " n l j of a/b "
    write(outputfile%u_outputEMelem,*) "--------------------------------------------------------"
    write(outputfile%u_outputEMelem,*) " ifg a/b  n   l   j"   
    do ifg =1, 2 
        do a = 1, TDs%nlj_length(ifg)
            i0sp = BS%HO_sph%iasp(1,ifg)
            a_index = TDs%nlj_index(a,ifg)
            nra= BS%HO_sph%nljm(i0sp+a_index,1) ! n_r
            nla= BS%HO_sph%nljm(i0sp+a_index,2) ! l
            nja= BS%HO_sph%nljm(i0sp+a_index,3) ! j +1/2
            write(outputfile%u_outputEMelem,"(5i4,'/2')") ifg, a, nra, nla, (nja*2-1)
        end do 
    end do

    ! store <a||Q_lambda||b>
    lambda_start = 1
    lambda_end = min(gcm_space%Jmax,TDs%lambda_max)
    header = '   a   b  sqrt(4pi)<a||r^2Y_0||b>'
    do lambda = lambda_start, lambda_end
        write(temp,'(A,I0,A)') '        <a||Q', lambda, '||b>'
        header = trim(header)//temp
    end do

    write(outputfile%u_outputEMelem,*) "------------------------------------------------------------------------"
    write(outputfile%u_outputEMelem,*) " reduced single-particle matrix element of the electromagnetic operator "
    write(outputfile%u_outputEMelem,*) "-------------------------------------------------------------------------"
    write(outputfile%u_outputEMelem,'(A)') trim(header)

    allocate(Ql_ab(lambda_start:lambda_end))
    ifg =  2 
    do a = 1, TDs%nlj_length(ifg)
        do b = 1, TDs%nlj_length(ifg)
            i0sp = BS%HO_sph%iasp(1,ifg)
            a_index = TDs%nlj_index(a,ifg)
            nra= BS%HO_sph%nljm(i0sp+a_index,1) ! n_r
            nla= BS%HO_sph%nljm(i0sp+a_index,2) ! l
            nja= BS%HO_sph%nljm(i0sp+a_index,3) ! j +1/2
            b_index = TDs%nlj_index(b,ifg)
            nrb= BS%HO_sph%nljm(i0sp+b_index,1) ! n_r
            nlb= BS%HO_sph%nljm(i0sp+b_index,2) ! l
            njb= BS%HO_sph%nljm(i0sp+b_index,3) ! j +1/2
            do lambda = lambda_start, lambda_end
                call reduced_multipole_matrix_elements(nra,nla,nja,lambda,nrb,nlb,njb,Ql_ab(lambda))
            end do
            call reduced_monopole_matrix_elements(nra,nla,nja,nrb,nlb,njb,monopole_ab)
            write(outputfile%u_outputEMelem,'(2i4,1x,f17.10,8x,15(1x,f17.10))') a,b, monopole_ab,(Ql_ab(lambda), lambda=lambda_start,lambda_end)
        end do 
    end do
end subroutine

END MODULE Proj_Inout