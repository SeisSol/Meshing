cmake_minimum_required(VERSION 3.10)

# set the project name
project(check_correctness)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(check_correctness src/main.cpp src/mesh.cpp)

target_include_directories(check_correctness PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

find_package(MPI REQUIRED)
target_link_libraries(check_correctness PUBLIC MPI::MPI_CXX)
target_compile_definitions(check_correctness PUBLIC USE_MPI)


find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(check_correctness PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(check_correctness PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})

