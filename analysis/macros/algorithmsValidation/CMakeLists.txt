message(STATUS "building algorithmsValidation ... ")


if    (ROOT_FOUND)
  ADD_EXECUTABLE(lightReconstruction            lightReconstruction.cc)
  TARGET_LINK_LIBRARIES(lightReconstruction     SuperFGD_Analysis ${ROOT_LIBRARIES})

  ADD_EXECUTABLE(tracking                    tracking.cc)
  TARGET_LINK_LIBRARIES(tracking     SuperFGD_Analysis ${ROOT_LIBRARIES})

endif (ROOT_FOUND)
