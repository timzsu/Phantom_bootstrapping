cmake_minimum_required(VERSION 3.26)

project(PhantomBootstrapping VERSION 1.0 LANGUAGES CXX CUDA)

# Apply the patch
execute_process(COMMAND (git apply "${PROJECT_SOURCE_DIR}/cmake/phantom-fhe.patch") || (git apply "${PROJECT_SOURCE_DIR}/cmake/phantom-fhe.patch" --reverse --check)
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extern/phantom-fhe")

# CMAKE config
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Debugging config
option(SHOW_WARNINGS "Enable warnings" OFF)
if (SHOW_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wno-comment")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -Wall,-Wextra,-Wpedantic,-Wno-comment")
endif()

# Language config
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# CUDA config
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES OR CMAKE_CUDA_ARCHITECTURES STREQUAL "")
  set(CMAKE_CUDA_ARCHITECTURES native) # auto-detect GPU
endif()
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 --default-stream per-thread")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Print CUDA info
message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CUDA Compiler: ${CMAKE_CUDA_COMPILER}")
message(STATUS "CUDA Compiler ID: ${CMAKE_CUDA_COMPILER_ID}")
message(STATUS "CUDA Compiler Version: ${CMAKE_CUDA_COMPILER_VERSION}")

# PhantomFHE config
option(PHANTOM_USE_CUDA_PTX "Use CUDA PTX Assembly" ON)
message(STATUS "Use CUDA PTX Assembly: ${PHANTOM_USE_CUDA_PTX}")
if (PHANTOM_USE_CUDA_PTX)
    add_compile_definitions(PHANTOM_USE_CUDA_PTX)
endif()

set(PHANTOM_INCLUDE_DIR 
${CMAKE_CURRENT_LIST_DIR}/extern/phantom-fhe/include ${CMAKE_SOURCE_DIR}/src
)

# PhantomFHE
file(GLOB_RECURSE PHANTOM_SRC CONFIGURE_DEPENDS "extern/phantom-fhe/src/*.cu")
add_library(Phantom SHARED ${PHANTOM_SRC})
target_include_directories(Phantom PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/extern/phantom-fhe/include>
  $<INSTALL_INTERFACE:include/phantom-fhe>
)
include_directories(${CMAKE_SOURCE_DIR}/extern/phantom-fhe/include)
target_compile_options(Phantom PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-src-in-ptx>")
target_compile_features(Phantom PUBLIC cxx_std_17 cuda_std_17)

# NexusCUDA Bootstrapping
file(GLOB_RECURSE BOOTSTRAP_SRC CONFIGURE_DEPENDS "src/bootstrapping/*.cu")
add_library(PhantomBootstrapping SHARED 
  ${BOOTSTRAP_SRC} 
  "src/ckks_evaluator.cu"
)
target_include_directories(PhantomBootstrapping PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/bootstrapping>
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/bootstrapping>
)
target_link_libraries(PhantomBootstrapping PUBLIC ntl gmp m pthread Phantom)

option(BUILD_EXAMPLES "Build Examples" On)
message(STATUS "Build Examples: ${BUILD_EXAMPLES}")
if (BUILD_EXAMPLES)
    add_subdirectory(example)
endif ()

# Install the libraries
install(
  TARGETS Phantom PhantomBootstrapping
  EXPORT PhantomBootstrappingTargets
  LIBRARY DESTINATION lib
)

# Install the headers
install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/extern/phantom-fhe/include/
  DESTINATION include/phantom-fhe
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/
  DESTINATION include
)

# Install the targets file
install(
  EXPORT PhantomBootstrappingTargets
  FILE PhantomBootstrappingTargets.cmake
  NAMESPACE PhantomBootstrapping::
  DESTINATION lib/cmake/PhantomBootstrapping
)

# Generate the package version file
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PhantomBootstrappingConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/PhantomBootstrappingConfig.cmake
  INSTALL_DESTINATION lib/cmake/PhantomBootstrapping
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/PhantomBootstrappingConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install the package configuration files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PhantomBootstrappingConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/PhantomBootstrappingConfigVersion.cmake
  DESTINATION lib/cmake/PhantomBootstrapping
)