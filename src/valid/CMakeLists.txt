set(LEVEL1_LIBS utils base)

FILE(GLOB sourcefiles "fdm/vpartial.f90")
add_executable(vpartial.x ${sourcefiles})
target_link_libraries(vpartial.x physics operators fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "fdm/vpartial3d.f90")
add_executable(vpartial3d.x ${sourcefiles})
target_link_libraries(vpartial3d.x physics operators fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "fdm/vintegral.f90")
add_executable(vintegral.x ${sourcefiles})
target_link_libraries(vintegral.x physics operators fdm ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "elliptic/velliptic.f90")
add_executable(velliptic.x ${sourcefiles})
target_link_libraries(velliptic.x physics operators fdm ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "burgers/*.f90")
# add_executable(vburgers.x ${sourcefiles})
# target_link_libraries(vburgers.x physics mappings operators thermodynamics filters fdm ibm ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "interpol/*.f90")
# add_executable(vinterpol.x ${sourcefiles})
# target_link_libraries(vinterpol.x operators filters fdm ibm ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "interpartial/*.f90")
# add_executable(vinterpartial.x ${sourcefiles})
# target_link_libraries(vinterpartial.x physics operators thermodynamics filters fdm ibm ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "integration/*.f90")
# add_executable(vintegration.x ${sourcefiles})
# target_link_libraries(vintegration.x ${LEVEL1_LIBS} ${LIBS})
