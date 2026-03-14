cmake_minimum_required(VERSION 3.10)

# set the project name
project(gmsh2gambit)

# specify the C++ standard (11 is enough for now, according to the SCons file)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

cmake_policy(SET CMP0074 NEW)


add_executable(gmsh2gambit
    src/main.cpp
    src/GMSH.cpp
    src/GAMBITWriter.cpp)

#logging
set(LOG_LEVEL "warning" CACHE STRING "Log level for the code")
set(LOG_LEVEL_OPTIONS "debug" "info" "warning" "error")
set_property(CACHE LOG_LEVEL PROPERTY STRINGS ${LOG_LEVEL_OPTIONS})
if("${LOG_LEVEL}" STREQUAL "debug")
  target_compile_definitions(gmsh2gambit PUBLIC LOG_LEVEL=3)
elseif("${LOG_LEVEL}" STREQUAL "info")
  target_compile_definitions(gmsh2gambit PUBLIC LOG_LEVEL=2)
elseif("${LOG_LEVEL}" STREQUAL "warning")
  target_compile_definitions(gmsh2gambit PUBLIC LOG_LEVEL=1)
elseif("${LOG_LEVEL}" STREQUAL "error")
  target_compile_definitions(gmsh2gambit PUBLIC LOG_LEVEL=0)
endif()

#add some compiler specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_compile_options(gmsh2gambit PUBLIC -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas>)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(gmsh2gambit PUBLIC -pedantic $<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANGUAGE:C>>:-Wall -w3 -diag-disable:remark>)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

target_include_directories(gmsh2gambit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules)
