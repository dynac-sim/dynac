
if(CMAKE_Fortran_COMPILER_ID EQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic")
endif()
add_executable(dynac dynac.F)
set_target_properties(dynac PROPERTIES VERSION ${DYNAC_VERSION})
install(TARGETS dynac RUNTIME DESTINATION bin COMPONENT Runtime)
