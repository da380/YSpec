add_library(YSpecModules module_fourier.f90
                         module_function.f90
                         module_int.f90
           		 module_model.f90
			 module_solver.f90
			 module_spline.f90
			 module_start.f90
			 module_util.f90
			 nrtype.f90
			 nrutil.f90)


add_executable(yspec yspec.f90)
target_link_libraries(yspec YSpecModules)

		 
