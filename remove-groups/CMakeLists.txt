cmake_minimum_required(VERSION 3.16)

project(remove-groups LANGUAGES C CXX)

option(sanitizer "Use address sanitizer when debugging" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo") # MinSizeRel is useless for us
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Set build type to Release as none was supplied.")
endif()


add_executable(remove-groups
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader.cpp
  )

target_include_directories(remove-groups PUBLIC src
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_include_directories(remove-groups PUBLIC SYSTEM src
        "${CMAKE_CURRENT_SOURCE_DIR}/submodules/"
        )

target_compile_options(remove-groups PUBLIC -Wall -Wextra)
target_link_libraries(remove-groups PUBLIC "-lstdc++fs")

set(HDF5_PREFER_PARALLEL True)
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(remove-groups PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(remove-groups PUBLIC ${HDF5_C_HL_LIBRARIES} ${HDF5_C_LIBRARIES})
target_compile_definitions(remove-groups PUBLIC USE_HDF)

find_package(MPI REQUIRED)
target_link_libraries(remove-groups PUBLIC MPI::MPI_CXX)
target_compile_definitions(remove-groups PUBLIC PARALLEL USE_MPI)

if (sanitizer)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
  set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
endif()

install(TARGETS remove-groups RUNTIME DESTINATION bin DESTINATION ${CMAKE_INSTALL_BINDIR})
