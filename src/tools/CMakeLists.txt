set(LEVEL1_LIBS utils base)

# Preprocessing
FILE(GLOB sourcefiles "preprocessing/inigrid/*.f90")
add_executable(inigrid.x ${sourcefiles})
target_link_libraries(inigrid.x operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "preprocessing/inirand/*.f90")
add_executable(inirand.x ${sourcefiles})
target_link_libraries(inirand.x physics operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "preprocessing/iniscal/*.f90")
add_executable(iniscal.x ${sourcefiles})
target_link_libraries(iniscal.x physics operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "preprocessing/iniflow/*.f90")
add_executable(iniflow.x ${sourcefiles})
target_link_libraries(iniflow.x physics mappings operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

# Simulation
FILE(GLOB sourcefiles "simulation/*.f90")
add_executable(dns.x ${sourcefiles})
target_link_libraries(dns.x statistics physics mappings operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

# Postprocessing
FILE(GLOB sourcefiles "postprocessing/visuals.f90")
add_executable(visuals.x ${sourcefiles})
target_link_libraries(visuals.x physics mappings operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

FILE(GLOB sourcefiles "postprocessing/averages*.f90")
add_executable(averages.x ${sourcefiles})
target_link_libraries(averages.x statistics physics mappings operators fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "statistics/spectra*.f90")
# add_executable(spectra.x ${sourcefiles})
# target_link_libraries(spectra.x statistics physics mappings operators filters fdm ibm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(spectra.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "statistics/pdfs*.f90")
# add_executable(pdfs.x ${sourcefiles})
# target_link_libraries(pdfs.x statistics physics mappings operators filters fdm ibm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(pdfs.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# # Transformation
# FILE(GLOB sourcefiles "transform/transgrid.f90")
# add_executable(transgrid.x ${sourcefiles})
# target_link_libraries(transgrid.x operators ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "transform/transfields.f90")
# add_executable(transfields.x ${sourcefiles})
# target_link_libraries(transfields.x statistics physics mappings operators filters fdm ibm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(transfields.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# # Particles
# FILE(GLOB sourcefiles "initialize/part/*.f90")
# add_executable(inipart.x ${sourcefiles})
# target_link_libraries(inipart.x statistics physics mappings operators particles filters fdm ibm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(inipart.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "initialize/part/particle_ini_trajec.f90")
# add_executable(l_ini_trajec.x ${sourcefiles})
# target_link_libraries(l_ini_trajec.x statistics physics mappings operators particles fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_ini_trajec.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "particles/particle_build_pdf.f90")
# add_executable(l_pdf.x ${sourcefiles})
# target_link_libraries(l_pdf.x statistics physics mappings operators particles filters fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_pdf.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "particles/particle_trajec.f90")
# add_executable(l_trajec.x ${sourcefiles})
# target_link_libraries(l_trajec.x statistics physics mappings operators particles filters fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_trajec.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# FILE(GLOB sourcefiles "particles/particle_pos_trajec*.f90")
# add_executable(l_pos_trajec.x ${sourcefiles})
# target_link_libraries(l_pos_trajec.x statistics physics mappings operators particles filters fdm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(l_pos_trajec.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# # Structure analysis
# FILE(GLOB sourcefiles "structure/apriori.f90")
# add_executable(apriori.x ${sourcefiles})
# target_link_libraries(apriori.x statistics physics mappings operators filters fdm ibm thermodynamics ${LEVEL1_LIBS} ${LIBS})
# if (PSFFT_LIB)
# set_target_properties(apriori.x PROPERTIES LINK_FLAGS "${PSFFT_LINK_FLAGS} ${LINK_FLAGS}")
# endif()

# # Cloud
# FILE(GLOB sourcefiles "cloud/state.f90")
# add_executable(state.x ${sourcefiles})
# target_link_libraries(state.x thermodynamics ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "cloud/smooth.f90")
# add_executable(smooth.x ${sourcefiles})
# target_link_libraries(smooth.x thermodynamics ${LEVEL1_LIBS} ${LIBS})

# FILE(GLOB sourcefiles "cloud/saturation.f90")
# add_executable(saturation.x ${sourcefiles})
# target_link_libraries(saturation.x thermodynamics ${LEVEL1_LIBS} ${LIBS})
