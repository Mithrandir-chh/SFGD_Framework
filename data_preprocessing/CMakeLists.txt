# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

# Set name of our project
project(unpacking-SuperFGD)


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")

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
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/macros)


install(FILES ${CMAKE_BINARY_DIR}/src/libEventDict_rdict.pcm 
	      ${CMAKE_BINARY_DIR}/src/libEventDict.rootmap 
              ${CMAKE_BINARY_DIR}/src/EventDict.cxx 
	DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)




