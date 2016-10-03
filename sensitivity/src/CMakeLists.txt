# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(SensitivityMVA VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-W -Wall -Wl,--no-as-needed -std=c++14")

# Will install module and support program
include(GNUInstallDirs)

find_package(ROOT 5.34.00 REQUIRED)
set(ROOT_COMPONENTS Minuit)
find_package(ROOT 5.34.00 REQUIRED ${ROOT_COMPONENTS})

# Ensure our code can see the ROOT etc headers
include_directories(${ROOT_INCLUDE_DIRS})

list(APPEND SensitivityMVA_HEADERS
  sensitivity_constants.h
  sensitivity_measurements.h
  analysis_config.h
  channel_selection.h
  pseudo_generator.h
  multi_fit.h
  stats_tools.h
  )

list(APPEND SensitivityMVA_SOURCES
  analysis_config.cc
  channel_selection.cc
  pseudo_generator.cc
  multi_fit.cxx
  main_sensitivity.cxx
  )

# Build a dynamic library from our sources
add_library(SensitivityMVA SHARED
  ${SensitivityMVA_HEADERS}
  ${SensitivityMVA_SOURCES}
  )

target_link_libraries(SensitivityMVA ${ROOT_LIBRARIES})

install(TARGETS SensitivityMVA DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(main_sensitivity main_sensitivity.cxx)
set_target_properties(main_sensitivity PROPERTIES COMPILE_DEFINITIONS ENABLE_BINRELOC)
target_link_libraries(main_sensitivity SensitivityMVA)
set_target_properties(main_sensitivity PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

install(TARGETS main_sensitivity RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
