###### Fortran standard and Compiler version
#  compile with Fortran 2008 standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# set(CMAKE_Fortran_EXTENSIONS OFF)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0)
        message(FATAL_ERROR "GNU Fortran >= 9.0 is required.")
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19.0)
        message(FATAL_ERROR "Intel Fortran >= 19.0 is required.")
    endif()
else()
    message(WARNING "Untested Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()