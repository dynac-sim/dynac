
if(CMAKE_Fortran_COMPILER_ID EQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic")
endif()
add_executable(tw2dyn tw2dyn.f)
set_target_properties(tw2dyn PROPERTIES VERSION ${DYNAC_VERSION})
install(TARGETS tw2dyn RUNTIME DESTINATION bin COMPONENT Runtime)
