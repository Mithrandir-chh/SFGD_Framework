message(STATUS "building muonDecay ... ")


if    (ROOT_FOUND)
  ADD_EXECUTABLE(muonDecay                     muonDecay.cc)
  TARGET_LINK_LIBRARIES(muonDecay              SuperFGD_Analysis ${ROOT_LIBRARIES})

  



  
endif (ROOT_FOUND)
