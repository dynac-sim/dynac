
#gfortran dyndat.f -o dynplt
add_executable(dynplt dyndat.f)
install(TARGETS dynplt RUNTIME DESTINATION bin COMPONENT Runtime)
