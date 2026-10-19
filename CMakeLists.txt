cmake_minimum_required(VERSION 3.10)

project(YSpec)
enable_language(Fortran)

option(YSPEC_USE_OPENMP "Parallelise the frequency loop using OpenMP" ON)
option(YSPEC_BUILD_TESTS "Build the regression tests" ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# FFTW3 is used for the transform to the time domain, through its
# Fortran 2003 interface (fftw3.f03)
find_path(FFTW3_INCLUDE_DIR fftw3.f03)
find_library(FFTW3_LIBRARY fftw3)
if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
  message(FATAL_ERROR "FFTW3 not found: install it (e.g. the libfftw3-dev "
                      "package), or set FFTW3_INCLUDE_DIR and FFTW3_LIBRARY")
endif()

add_library(YSpecModules src/module_kinds.f90
                         src/module_util.f90
                         src/module_spline.f90
                         src/module_model.f90
                         src/module_function.f90
                         src/module_start.f90
                         src/module_int.f90
                         src/module_solver.f90
                         src/module_fourier.f90)
target_include_directories(YSpecModules PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(YSpecModules PUBLIC ${FFTW3_LIBRARY})

if(YSPEC_USE_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(YSpecModules PUBLIC OpenMP::OpenMP_Fortran)
  else()
    message(STATUS "OpenMP not found: yspec will run serially")
  endif()
endif()

add_custom_target(CopyExamples
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/examples
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )

add_executable(yspec src/yspec.f90)
target_link_libraries(yspec YSpecModules)
add_dependencies(yspec CopyExamples)

if(YSPEC_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
