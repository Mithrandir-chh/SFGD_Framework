message(STATUS "building pionReco ... ")


if    (ROOT_FOUND)
  ADD_EXECUTABLE(pionReco                     pionReco.cc)
  TARGET_LINK_LIBRARIES(pionReco              SuperFGD_Analysis ${ROOT_LIBRARIES})

  



  
endif (ROOT_FOUND)
