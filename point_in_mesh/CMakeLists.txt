cmake_minimum_required(VERSION 3.9)

project(point_in_mesh LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()

add_executable(point_in_mesh 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader.cpp
  )

target_include_directories(point_in_mesh PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

set(HDF5_PREFER_PARALLEL False)
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(point_in_mesh PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(point_in_mesh PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})

if(HDF5_IS_PARALLEL)
  find_package(MPI REQUIRED)
  target_link_libraries(point_in_mesh PUBLIC MPI::MPI_CXX)
  target_compile_definitions(point_in_mesh PUBLIC PARALLEL)
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(point_in_mesh PUBLIC OpenMP::OpenMP_CXX)

find_package(Eigen3 3.4 REQUIRED)
target_link_libraries(point_in_mesh PUBLIC Eigen3::Eigen)

function(cast_build_type_to_int build_type log_level_int)
  if (${build_type} STREQUAL "Debug")
    set(${log_level_int} 3 PARENT_SCOPE)
  elseif (${build_type} STREQUAL "RelWithDebInfo")
    set(${log_level_int} 2 PARENT_SCOPE)
  elseif (${build_type} STREQUAL "Release")
    set(${log_level_int} 2 PARENT_SCOPE)
  endif()
endfunction()

cast_build_type_to_int(CMAKE_BUILD_TYPE LOG_LEVEL)
target_compile_definitions(point_in_mesh PUBLIC LOG_LEVEL=${LOG_LEVEL})

install(TARGETS point_in_mesh RUNTIME DESTINATION bin DESTINATION ${CMAKE_INSTALL_BINDIR})
