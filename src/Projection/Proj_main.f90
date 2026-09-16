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
        use Tools, only: int2str
        use Globals, only: gcm_space,Proj_option,Proj_outputfile,MPI_Infor,Factors
        use Energy, only: calculate_sigma_nabla_Spherical
        use Basis, only: set_spherical_oscillator_wave_function
        use Proj_Inout, only: set_Proj_Expectation_filename,read_wavefuntion_files_with_factor,write_Proj_output_with_fact
        use Mixed, only: determine_truncated_dimension
        use Kernel, only: set_projection_mesh_points,calculate_Kernel
        use Proj_Density, only: calculate_density_matrix_element
        use TD, only: calculate_reduced_transition_density_matrix_element
        use Constraints, only: select_the_constraint_deformation_space
        use Forces, only: set_force_parameters_with_factor
        implicit none
        integer :: factor_index1,factor_index2,factor_index2_start,factor_index2_end,q1,q2,q2_start, q2_end
        integer :: task_id

        call set_projection_mesh_points
        call set_spherical_oscillator_wave_function
        call calculate_sigma_nabla_Spherical

        call set_Proj_Expectation_filename
        open(Proj_outputfile%u_outExpectation, file=trim(Proj_outputfile%outExpectation)//'_'//int2str(MPI_Infor%rank), status='unknown')

        task_id = -1
        do factor_index1 = Factors%start_index, Factors%end_index
            if(Proj_option%Kernel_Symmetry==0) then      ! All Kernels 
                factor_index2_start = Factors%start_index
                factor_index2_end = Factors%end_index
            else if(Proj_option%Kernel_Symmetry==1) then ! Triangular Kernels
                factor_index2_start = factor_index1
                factor_index2_end = Factors%end_index
            else if(Proj_option%Kernel_Symmetry==2) then ! Diagonal elements
                factor_index2_start = factor_index1
                factor_index2_end = factor_index1
            else
                stop 'Kernel_Symmetry should be 0, 1, 2'
            end if
            do factor_index2 = factor_index2_start, factor_index2_end
                task_id = task_id + 1
                if (mod(task_id, MPI_Infor%nprocs) /= MPI_Infor%rank) cycle ! MPI process distribution
                if(factor_index1==factor_index2) then 
                    call set_force_parameters_with_factor(factor_index1, .True.)
                else 
                    stop 'Correct H have not implement yet!'
                end if 
                call select_the_constraint_deformation_space(factor_index1,factor_index2,Factors%deformation_option)
                do q1 = gcm_space%q1_start, gcm_space%q1_end
                    if(Proj_option%Kernel_Symmetry==0) then      ! All Kernels 
                        q2_start = gcm_space%q2_start
                        q2_end = gcm_space%q2_end
                    else if(Proj_option%Kernel_Symmetry==1) then ! Triangular Kernels
                        if(factor_index1==factor_index2) then 
                            q2_start = q1
                            q2_end = gcm_space%q2_end
                        else 
                            q2_start = gcm_space%q2_start
                            q2_end = gcm_space%q2_end
                        end if 
                    else if(Proj_option%Kernel_Symmetry==2) then ! Diagonal elements
                        q2_start = q1
                        q2_end = q1
                    else
                        stop 'Kernel_Symmetry should be 0, 1, 2'
                    end if
                    do q2 = q2_start, q2_end
                        write(*,"(5x,'(f1,f2,q1,q2)=(',i4,',',i4,',',i3,',',i3,') >>>>>>>>>>')") factor_index1,factor_index2,q1, q2
                        call read_wavefuntion_files_with_factor(factor_index1,factor_index2,q1,q2)
                        call determine_truncated_dimension
                        call calculate_Kernel
                        ! if (Proj_option%DsType > 0) then
                        !     call calculate_density_matrix_element(q1,q2)
                        ! end if 
                        ! if (Proj_option%AMPtype > 0 .and. Proj_option%TDType > 0) then
                        !     call calculate_reduced_transition_density_matrix_element(q1,q2)
                        ! end if 
                        call write_Proj_output_with_fact(factor_index1,factor_index2,q1,q2)
                    end do
                end do
            end do
        end do
        close(Proj_outputfile%u_outExpectation) 
    end subroutine
END MODULE Proj