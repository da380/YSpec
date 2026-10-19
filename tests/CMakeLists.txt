add_executable(yspec_compare compare.f90)

# Each case directory holds an input file and the reference outputs that
# it produced. Differences are measured relative to the peak amplitude of
# each component, so the tolerance allows for rounding in the 6-digit
# output but catches any change in the physics.
set(YSPEC_TEST_TOLERANCE 1e-4)

function(yspec_add_case name model)
  add_test(NAME ${name}
           COMMAND ${CMAKE_COMMAND}
                   -DYSPEC=$<TARGET_FILE:yspec>
                   -DCOMPARE=$<TARGET_FILE:yspec_compare>
                   -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${name}
                   -DMODEL=${PROJECT_SOURCE_DIR}/examples/${model}
                   -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${name}
                   -DTOL=${YSPEC_TEST_TOLERANCE}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/RunCase.cmake)
endfunction()

# full gravity, solid inner core, attenuation, displacement; receivers from
# 5 to 175 degrees and at the north pole
yspec_add_case(selfgrav_disp prem.200)

# Cowling, solid inner core, attenuation, acceleration with potential and
# tilt corrections; source in the low-Q zone and a receiver at the south pole
yspec_add_case(cowling_accel prem.200)

# no gravity, fluid inner core, no attenuation, velocity; receiver below the
# source
yspec_add_case(fic_nograv_vel prem.200.fic)

# radial modes only (lmax = 0) from an isotropic source, acceleration with the
# free-air correction, which is significant near 0S0
yspec_add_case(radial_accel prem.200)
