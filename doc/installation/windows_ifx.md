# Installation Guide (Windows + ifx)

This guide describes how to build the project on Windows using CMake, Visual Studio, and Intel oneAPI.

## Prerequisites

Before building the project, please install the following dependencies:

1. **CMake**  
   Download and install CMake from:  
   https://cmake.org/download/

2. **Visual Studio**  
   Visual Studio is required to provide the Microsoft linker (link.exe) and Windows SDK, which ifx relies on for linking Fortran code into an executable.    If you have already installed Visual Studio, please verify that the version meets the requirements of the Intel oneAPI Toolkit.
   If not installed, download and install Visual Studio Community from:
   https://visualstudio.microsoft.com/zh-hans/downloads/

   During installation, make sure to select the following workload:
   - **Desktop development with C++**

   This workload provides the required MSVC compiler and build tools.

3. **Intel oneAPI Toolkit**  
   Intel oneAPI Toolkit provides the required Fortran compiler (`ifx`), Intel MPI Library, and Intel Math Kernel Library (MKL) needed for building the project. Download and install Intel oneAPI Toolkit from:  
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/oneapi-toolkit-download.html

>  [!CAUTION]
> * Please remember the installation paths of Visual Studio and Intel oneAPI Toolkit, as they will be needed during the later compilation process.
> * After installing CMake, the environment variables are typically configured automatically. Open a terminal and run `cmake --version` to verify it works. If the command is not recognized, manually add CMake's bin directory (e.g., `C:\Program Files\CMake\bin`) to your system's PATH.
> * After installing the Intel oneAPI Toolkit, the environment variable is typically configured automatically. Open a terminal and run `mpiexec --version` to verify `mpiexec` works. If the command is not recognized, manually add the bin directory where mpiexec is located (e.g., `C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin`) to your system PATH.

## Environment Setup

After installation, open a **Command Prompt** (cmd) terminal and set the Visual Studio and Intel oneAPI environments:

```cmd
call "C:\Program Files (x86)\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```
> [!IMPORTANT]
> The installation paths may differ depending on your local setup. Please replace the paths above with the actual locations of your `vcvars64.bat` and `setvars.bat` files.

These two commands set the required compilation environment variables for the current terminal session.
- `vcvars64.bat` sets the Visual Studio C/C++ build environment, including the MSVC compiler, linker, Windows SDK, and related build tools.
- `setvars.bat` sets the Intel oneAPI environment, including the Intel Fortran compiler (`ifx`/`ifort`), Intel MPI Library, Intel MKL, and other Intel development tools.



---

## Build Instructions
Use the **Command Prompt** terminal in which the Visual Studio and Intel oneAPI environments have already been set, and navigate to the root directory of the MRCDFT project.

First, verify that the compilation environment has been configured correctly by running:
```cmd
cmake --preset mpi-ifx
```
This command configures the project using the predefined CMake preset mpi-ifx. It detects the Intel Fortran compiler, MPI environment, MKL libraries, and generates the corresponding build files in the `build\mpi-ifx\` directory.

If an error occurs, it usually means that the current terminal environment is not configured correctly. In this case, please check whether Visual Studio and Intel oneAPI Toolkit have been installed properly, and make sure the environment initialization commands were executed successfully.

If the command completes successfully, the compilation environment is ready. Then build the project by running:

```cmd
cmake --build --preset mpi-ifx
```
This command compiles the source code and builds the executable files according to the configuration defined in the mpi-ifx preset.

After the build completes successfully, the executable file `MRCDFT.exe` will be generated in the project's `bin\` directory.

To run `MRCDFT.exe` from any directory in the command line, you need to add the directory containing `MRCDFT.exe` to your system `PATH` environment variable.

## Running a Test Calculation

After `MRCDFT.exe` has been successfully built, you can either continue using the Command Prompt with the environments set above, or open a new terminal and navigate to the root directory of the MRCDFT project.

#### Set the number of threads

If you are using a **Command Prompt (cmd)** terminal, run:

```cmd
set OMP_NUM_THREADS=4
set MKL_NUM_THREADS=4
```
These commands set the number of OpenMP and MKL threads used by each MPI process.

If you are using a **Power Shell** (PS) terminal, run:
```ps
$env:OMP_NUM_THREADS=4
$env:MKL_NUM_THREADS=4
```

#### Run with multiple processes
To run a test calculation for $^{22}\mathrm{Ne}$, execute the test calculation with:
```bash
cd examples/22Ne
mpiexec -np 2 ../../bin/MRCDFT -p 22Ne_para.dat -d 22Ne_b23.dat
```
This example will run with 2 MPI processes, with each process using 4 threads.
