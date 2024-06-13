cmake_minimum_required(VERSION 3.5)

project(YSpec)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)


add_library(YSpecModules src/module_fourier.f90
                         src/module_function.f90
                         src/module_int.f90
           		         src/module_model.f90
			             src/module_solver.f90
			             src/module_spline.f90
			             src/module_start.f90
			             src/module_util.f90
			             src/nrtype.f90
			            src/nrutil.f90)







add_custom_target(CopyExamples        
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/examples
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples        
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )

add_executable(yspec src/yspec.f90)
target_link_libraries(yspec YSpecModules)
add_dependencies(yspec CopyExamples)

#add_custom_command(TARGET bar
#    # Run after all other rules within the target have been executed
#    POST_BUILD
#    COMMAND echo "executing a POST_BUILD command"
#    COMMENT "This command will be executed after building bar"
#    VERBATIM
#)


