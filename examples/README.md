# Examples :

## Visual Studio:
In the Visual Studio Project the solution is provided along with the source file and the 
dependencies. The solution was implemented in Microsoft Visual Studio 2019.Simply open the solution in 
Microsoft Visual Studio 2019

## CMake: 
The CMakeLists file is also provided to create the project with CMake. Armadillo comes with an implemented library of OPENBlas. In case that you have LAPACK or BLAS and want to use those, or you do not have these libraries, configure the "config.hpp" file in armadillo/include/armadillo_bits as per the instructions so that no errors will occur. After that, use CMake to make the .exe in the build directory of your choosing and run the .exe to get the testing results.