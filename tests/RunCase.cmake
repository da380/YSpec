# Runs one regression case: runs yspec on the case's input file in a fresh
# working directory, then compares each output file with the reference of
# the same name. Invoked by ctest with -P; see tests/CMakeLists.txt.

file(REMOVE_RECURSE ${WORK_DIR})
file(MAKE_DIRECTORY ${WORK_DIR})
file(COPY ${CASE_DIR}/input.txt ${MODEL} DESTINATION ${WORK_DIR})

execute_process(COMMAND ${YSPEC} input.txt
                WORKING_DIRECTORY ${WORK_DIR}
                OUTPUT_FILE ${WORK_DIR}/yspec.log
                RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "yspec exited with status ${status}; see ${WORK_DIR}/yspec.log")
endif()

file(GLOB references ${CASE_DIR}/reference/*)
if(NOT references)
  message(FATAL_ERROR "no reference files in ${CASE_DIR}/reference")
endif()

set(failed FALSE)
foreach(reference ${references})
  get_filename_component(name ${reference} NAME)
  execute_process(COMMAND ${COMPARE} ${WORK_DIR}/${name} ${reference} ${TOL}
                  RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    set(failed TRUE)
  endif()
endforeach()

if(failed)
  message(FATAL_ERROR "outputs differ from the references")
endif()
