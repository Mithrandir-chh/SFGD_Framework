message(STATUS "building template ... ")

if    (ROOT_FOUND)
  #ADD_EXECUTABLE(templateAna                     templateAna.cc)
  #TARGET_LINK_LIBRARIES(templateAna              SuperFGD_Analysis ${ROOT_LIBRARIES})  
endif (ROOT_FOUND)
