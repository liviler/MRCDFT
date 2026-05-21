    You can compile the code using either gfortran or Intel Fortran (ifort), depending on the compiler available on your system.
    * Using CMake 
        ```bash
        cmake --preset mpi-gfortran
        cmake --build --preset mpi-gfortran
        ```
    * Using Make
        ```bash
        make mpif90
        ```
        Before running the make command, ensure that GNU Make and the selected Fortran compiler are properly installed and available in your environment.

    After successful compilation, the executable, the executable `MRCDFT` will be generated in the `bin/` directory.

3. Adding to Environment Variables :

    To run `MRCDFT` from any directory, you need to add the executable path to your system’s PATH environment variable.
    * Linux / macOS
        1) Open your shell configuration file (e.g. `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`).
        2) Add the following line (replace the path with the actual location of the bin/ directory):
            ```bash
            export PATH="/full/path/to/MR_CDFT_f90/bin:$PATH"
            ```
        3) Reload the configuration file:
            ``` bash
            source ~/.bashrc
            ```
        4) Verify the installation:
            ```bash
            which MRCDFT
            ```
    * Windows
        1) Open System Properties → Advanced system settings → Environment Variables.
        2) Under User variables or System variables, select Path and click Edit.
        3) Add the full path to the bin/ directory containing MRCDFT.exe.
        4) Open a new command prompt and test:
            ```cmd
                MRCDFT
            ```