# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( triangulation_ )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # Including local libraries
  set(CMAKE_CXX_FLAGS "-fPIC -pthread -std=c++11 -Wall")
  add_library(Lock Lock.cpp Lock.h)
  add_library(Point Point.cpp Point.h)
  add_library(Vertex Vertex.cpp Vertex.h)
  add_library(Face Face.cpp Face.h)
  add_library(Face_circulator Face_circulator.cpp Face_circulator.h)
  add_library(Triangulation Triangulation.cpp Triangulation.h)
  add_library(Delaunay_triangulation Delaunay_triangulation.cpp Delaunay_triangulation.h)
  add_library(predicates predicates.c)
  set_target_properties(predicates PROPERTIES COMPILE_FLAGS "-fPIC")
  target_link_libraries(Vertex Lock)
  target_link_libraries(Vertex Point)
  target_link_libraries(Vertex Face)
  target_link_libraries(Face Vertex)
  target_link_libraries(Face_circulator Face)
  target_link_libraries(Triangulation Face)
  target_link_libraries(Triangulation Face_circulator)
  target_link_libraries(Triangulation predicates)
  target_link_libraries(Delaunay_triangulation Triangulation)

  create_single_source_cgal_program( "triangulate.cpp" )
  target_link_libraries(triangulate Delaunay_triangulation)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

