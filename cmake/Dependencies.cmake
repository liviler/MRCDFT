####### Dependence
# MPI
find_package(MPI REQUIRED COMPONENTS Fortran)

# OpenMP
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# MKL
set(MKL_LINK static)
set(MKL_INTERFACE lp64)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(MKL_THREADING gnu_thread)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(MKL_THREADING intel_thread)
else()
    set(MKL_THREADING sequential)
endif()
find_package(MKL CONFIG REQUIRED)