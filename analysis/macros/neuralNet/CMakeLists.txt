message(STATUS "building neuralNet ... ")


if    (ROOT_FOUND)
  ADD_EXECUTABLE(inputNN                     inputNN.cc)
  TARGET_LINK_LIBRARIES(inputNN              SuperFGD_Analysis ${ROOT_LIBRARIES})

  



  
endif (ROOT_FOUND)
