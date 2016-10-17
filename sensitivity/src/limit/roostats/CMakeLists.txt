# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(limit_roostats VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-W -Wall -Wl,--no-as-needed -std=c++14")

# Will install module and support program
include(GNUInstallDirs)

find_package(ROOT 6 REQUIRED)
set(ROOT_COMPONENTS Minuit)
find_package(ROOT 6 REQUIRED ${ROOT_COMPONENTS})

# Ensure our code can see the ROOT etc headers
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})

# list(APPEND limit_roostats_HEADERS)

list(APPEND limit_roostats_SOURCES
  rs101_limitexample.C
)

# Build a dynamic library from our sources
add_library(limit_roostats SHARED
  # ${limit_roostats_HEADERS}
  ${limit_roostats_SOURCES}
)

target_link_libraries(limit_roostats ${ROOT_LIBRARIES})

install(TARGETS limit_roostats DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(test_limit rs101_limitexample.C)
set_target_properties(test_limit PROPERTIES COMPILE_DEFINITIONS ENABLE_BINRELOC)
target_link_libraries(test_limit limit_roostats)
set_target_properties(test_limit PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

install(TARGETS test_limit RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
