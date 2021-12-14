cmake_minimum_required(VERSION 3.9)

project(place_receivers LANGUAGES C CXX)

set(LOG_LEVEL "info" CACHE STRING "Log level for the code")
set_property(CACHE LOG_LEVEL PROPERTY STRINGS "debug" "info" "warning" "error")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()


add_executable(place_receivers
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader.cpp
  )

target_include_directories(place_receivers PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

set(HDF5_PREFER_PARALLEL False)
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(place_receivers PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(place_receivers PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
target_compile_definitions(place_receivers PUBLIC USE_HDF)

if(HDF5_IS_PARALLEL)
  find_package(MPI REQUIRED)
  target_link_libraries(place_receivers PUBLIC MPI::MPI_CXX)
  target_compile_definitions(place_receivers PUBLIC PARALLEL)
endif()

function(cast_log_level_to_int log_level_str log_level_int)
  if (${log_level_str} STREQUAL "debug")
    set(${log_level_int} 3 PARENT_SCOPE)
  elseif (${log_level_str} STREQUAL "info")
    set(${log_level_int} 2 PARENT_SCOPE)
  elseif (${log_level_str} STREQUAL "warning")
    set(${log_level_int} 1 PARENT_SCOPE)
  elseif (${log_level_str} STREQUAL "error")
    set(${log_level_int} 0 PARENT_SCOPE)
  endif()
endfunction()

cast_log_level_to_int(LOG_LEVEL LOG_LEVEL)
target_compile_definitions(place_receivers PUBLIC LOG_LEVEL=${LOG_LEVEL})

install(TARGETS place_receivers RUNTIME DESTINATION bin DESTINATION ${CMAKE_INSTALL_BINDIR})
