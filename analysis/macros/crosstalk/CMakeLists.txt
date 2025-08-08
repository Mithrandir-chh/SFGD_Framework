message(STATUS "building crosstalk ... ")


if    (ROOT_FOUND)
  ADD_EXECUTABLE(crosstalkDist                     crosstalkDist.cc)
  TARGET_LINK_LIBRARIES(crosstalkDist              SuperFGD_Analysis ${ROOT_LIBRARIES})

 
  
endif (ROOT_FOUND)
