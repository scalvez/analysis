# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(Bayes VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-W -Wall -Wl,--no-as-needed -std=c++14")

# Will install module and support program
include(GNUInstallDirs)

find_package(ROOT 5.34.00 REQUIRED)
set(ROOT_COMPONENTS Minuit)
find_package(ROOT 5.34.00 REQUIRED ${ROOT_COMPONENTS})

# Ensure our code can see the ROOT etc headers
include_directories(${ROOT_INCLUDE_DIRS})

list(APPEND Bayes_HEADERS
  analysis_config.h
  stats_tools.h
  pseudo_generator.h
  )

list(APPEND Bayes_SOURCES
  main_bayes.cxx
  # main_bayes_example.cxx
  pseudo_generator.cc
  )

# Build a dynamic library from our sources
add_library(Bayes SHARED
  ${Bayes_HEADERS}
  ${Bayes_SOURCES}
  )

target_link_libraries(Bayes ${ROOT_LIBRARIES})

install(TARGETS Bayes DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(main_bayes main_bayes.cxx)
# add_executable(main_bayes main_bayes_example.cxx)
set_target_properties(main_bayes PROPERTIES COMPILE_DEFINITIONS ENABLE_BINRELOC)
target_link_libraries(main_bayes Bayes)
set_target_properties(main_bayes PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

install(TARGETS main_bayes RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
