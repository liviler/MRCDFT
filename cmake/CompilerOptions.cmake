####### Additional option
# Set compiler options for different build types
set(COMPILER_OPTIONS "")
set(LINKER_OPTIONS "")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(COMPILER_OPTIONS 
        -ffree-line-length-none
        # -fheap-arrays
        -cpp
        $<$<CONFIG:Debug>:-g;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow;-Wall>
    )
    # Windows gfortran stack size
    if(WIN32)
        set(LINKER_OPTIONS -Wl,--stack,100000000)
    endif()
    
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    if(WIN32) # Windows
        set(COMPILER_OPTIONS
            /free
            /heap-arrays
            /fpp
            $<$<CONFIG:Debug>:/debug:full;/check:all;/traceback;/warn:all>
        )
        set(LINKER_OPTIONS
            /libs:static
            /Qopenmp-link:static
            /STACK:100000000
        )
    else() # Linux/macOS
        set(COMPILER_OPTIONS
            -free
            -heap-arrays
            -fpp
            $<$<CONFIG:Debug>:-g;-check all;-traceback;-warn all>
        )
    endif()
else()
    set(COMPILER_OPTIONS $<$<CONFIG:Debug>:-g>)
endif()