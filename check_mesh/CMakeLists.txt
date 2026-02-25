cmake_minimum_required(VERSION 3.10)

# set the project name
project(check_mesh)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(check_mesh src/main.cpp src/mesh.cpp)

target_include_directories(check_mesh PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

find_package(MPI REQUIRED)
target_link_libraries(check_mesh PUBLIC MPI::MPI_CXX)
target_compile_definitions(check_mesh PUBLIC USE_MPI)


find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(check_mesh PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(check_mesh PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})

