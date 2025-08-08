# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

# Set name of our project
project(SuperFGD_Analysis)
enable_testing()
include (CTest)
set(CTEST_MEMORYCHECK_COMMAND, "usr/bin/valgrind")
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS, "--leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11 -m64")

INCLUDE_DIRECTORIES(. ${CMAKE_CURRENT_SOURCE_DIR}/src)
LINK_DIRECTORIES(     ${CMAKE_CURRENT_SOURCE_DIR}/lib)


find_package(ROOT CONFIG REQUIRED Graf Gpad Hist RIO Tree)

if (ROOT_FOUND)
  message(STATUS "ROOT found")
  INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
  INCLUDE(${ROOT_USE_FILE}) 
  LINK_DIRECTORIES   (${ROOT_LIBRARY_DIR})
Endif (ROOT_FOUND)


message(STATUS "C flags      : " ${CMAKE_C_FLAGS})
message(STATUS "CXX flags    : " ${CMAKE_CXX_FLAGS})
message(STATUS "Linker flags : " ${CMAKE_EXE_LINKER_FLAGS}\n)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

MACRO(ADD_PROJECT PROJECT)
    #IF(${PROJECT} OR ALL)
    set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin/${PROJECT}")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/macros/${PROJECT})
    #ENDIF()
ENDMACRO()

#IF YOU ARE A NEW USER, COPY THE LINE BELOW AND REPLACE "template" 
#BY THE NAME OF YOUR PERSONAL DIRECTORY IN /macros/

#ADD_PROJECT(template)
ADD_PROJECT(wilfs)
ADD_PROJECT(che)
#ADD_PROJECT(neuralNet)
#ADD_PROJECT(algorithmsValidation)
#ADD_PROJECT(crosstalk)
#ADD_PROJECT(muonDecay)
#ADD_PROJECT(pionReco)


install(FILES ${CMAKE_BINARY_DIR}/src/libEventDict_rdict.pcm 
	        ${CMAKE_BINARY_DIR}/src/libEventDict.rootmap 
          ${CMAKE_BINARY_DIR}/src/EventDict.cxx 
	        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
