# First we set up test macro
# In any subdirectory, you specify a test by
# run_test(test_name input)
#  - test_name is the name of your test (must be unique)
#  - input is the name of the input file
#
# For any .cfg file in the example directory, there should
# also be a .ref file which contains the expected output
# if the numdiff binary is found on your system, the example
# output is then compared to the reference output .ref, according
# to the configured rules in .cfg (by default identical files expected)
#
# See manual http://cds.cern.ch/record/1574594/

# Find numdiff binary..
find_program(NUMDIFF_BIN NAMES numdiff)

macro(run_test test_name input)
   # listing of different files..
   file(GLOB CFG_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.cfg")
   file(GLOB ALL_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*")
   string(REGEX REPLACE ".cfg$" "" TEST_OUTPUT "${CFG_FILES}")

   # copy all input files to working directory..
   file(COPY ${ALL_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

   # make sure all output files are removed (from previous runs)..
   file(REMOVE ${OLD_OUTPUT} ${TEST_OUTPUT})

   # run the test itself..
   add_test(NAME ${test_name} COMMAND dynac ${input})
   set_tests_properties(${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION "STOP ON KEY")

   if(NUMDIFF_BIN)
      # run numdiff comparison..
      add_test(NAME ${test_name}_numdiff COMMAND numdiff -b -l -t ${test_name} ${TEST_OUTPUT})
      set_tests_properties (${test_name}_numdiff PROPERTIES
          PASS_REGULAR_EXPRESSION ".*${test_name}.*PASS"
          DEPENDS ${test_name})
   endif()
endmacro()

add_subdirectory(egun)

