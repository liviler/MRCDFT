!==============================================================================!
! MODULE Proj                                                                  !
!                                                                              !
! This module calculates the                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine                                                                 !
!==============================================================================!
MODULE Proj
implicit none
contains
    subroutine Proj_Main
        use Globals, only: gcm_space,pko_option
        use Energy, only: calculate_sigma_nabla_Spherical
        use Basis, only: set_spherical_oscillator_wave_function
        use Proj_Inout, only: read_wavefuntion_files,write_pko_output
        use Mixed, only: determine_truncated_dimension
        use Kernel, only: set_projection_mesh_points,calculate_Kernel
        use Proj_Density, only: calculate_density_matrix_element
        use TD, only: calculate_reduced_transition_density_matrix_element
        implicit none
        integer :: q1,q2_start,q2_end,q2

        write(*,*)'PKO_Main: Start PKO calculations'
        call set_projection_mesh_points
        call set_spherical_oscillator_wave_function
        call calculate_sigma_nabla_Spherical
        do q1 = gcm_space%q1_start, gcm_space%q1_end
            if(pko_option%Kernel_Symmetry==0) then      ! All Kernels 
                q2_start = gcm_space%q2_start
                q2_end = gcm_space%q2_end
            else if(pko_option%Kernel_Symmetry==1) then ! Triangular Kernels
                q2_start = q1
                q2_end = gcm_space%q2_end
            else if(pko_option%Kernel_Symmetry==2) then ! Diagonal elements
                q2_start = q1
                q2_end = q1
            else
                stop 'Kernel_Symmetry should be 0, 1, 2'
            end if
            do q2 = q2_start, q2_end
                write(*,"(5x,'(q1,q2)=(',i3,',',i3,')')") q1, q2
                call read_wavefuntion_files(q1,q2)
                call determine_truncated_dimension
                call calculate_Kernel
                if (pko_option%DsType == 1) then
                    call calculate_density_matrix_element(q1,q2)
                end if 
                if (pko_option%TDType == 1) then
                    call calculate_reduced_transition_density_matrix_element(q1,q2)
                end if 
                call write_pko_output(q1,q2)
            end do
        end do

        write(*,*)'PKO_Main: PKO calculations completed'
    end subroutine
END MODULE Proj