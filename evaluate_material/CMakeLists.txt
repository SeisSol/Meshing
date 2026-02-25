cmake_minimum_required(VERSION 3.10)

# set the project name
project(evaluate_easi)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(evaluate_easi src/main.cpp src/mesh.cpp src/parameterDB.cpp src/writer.cpp)

target_include_directories(evaluate_easi PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
)

find_package(MPI REQUIRED)

target_include_directories(evaluate_easi SYSTEM PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(evaluate_easi PUBLIC MPI::MPI_C)

find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(evaluate_easi PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(evaluate_easi PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})

find_package(easi 1.0.0 REQUIRED)
target_link_libraries(evaluate_easi PUBLIC easi::easi)


