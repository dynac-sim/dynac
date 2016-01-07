
#gfortran dyndat.f -o dynplt
add_executable(dynplt dyndat.F)
install(TARGETS dynplt RUNTIME DESTINATION bin COMPONENT Runtime)
