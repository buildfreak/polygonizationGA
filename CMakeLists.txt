# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
#set (CMAKE_CXX_FLAGS "-fopenmp -DTRACE")
#set (CMAKE_CXX_FLAGS "")
set (CMAKE_CXX_FLAGS "-fopenmp")
#set (CMAKE_CXX_FLAGS "-DNDEBUG -fopenmp")

project(genetic_main)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBS OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )


  file(GLOB SRC
    "*.h"
    "*.cpp"
  )
  #adding all .h and .cpp to target  
  add_executable(${PROJECT_NAME} ${SRC})
  
  #inclusions
  include_directories ("/usr/include")
  link_directories ("/usr/lib/")

#linking xercesc
  include(FindXercesC)

  find_package(XercesC REQUIRED)
  include_directories( ${XercesC_INCLUDE_DIR} )
  target_link_libraries (${PROJECT_NAME} ${XercesC_LIBRARY} )

#link boost
  find_package(Boost 1.58.0 COMPONENTS program_options timer filesystem REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries (${PROJECT_NAME} ${Boost_LIBRARIES})


target_compile_features(${PROJECT_NAME} PRIVATE cxx_range_for)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
