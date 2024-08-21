# cmake-format: on

cmake_minimum_required(VERSION 3.24)

# 3.27: find_package() uses upper-case <PACKAGENAME>_ROOT variables
cmake_policy(SET CMP0144 NEW)

project(secondryVertex)

# CMake includes
include(CheckCXXCompilerFlag)

# Set default standard to C++20
set(CMAKE_CXX_STANDARD_MIN 20)
set(CMAKE_CXX_STANDARD
    ${CMAKE_CXX_STANDARD_MIN}
    CACHE STRING "C++ standard to be used")
if(CMAKE_CXX_STANDARD LESS CMAKE_CXX_STANDARD_MIN)
  message(
    FATAL_ERROR
      "Unsupported C++ standard: ${CMAKE_CXX_STANDARD} (at least ${CMAKE_CXX_STANDARD_MIN} required)"
  )
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Algorithms
find_package(algorithms 1.0.0 REQUIRED Core)

# PODIO, EDM4HEP, EDM4EIC event models
find_package(Eigen3 REQUIRED)
find_package(podio REQUIRED)
find_package(EDM4HEP 0.7.1 REQUIRED)
find_package(EDM4EIC 5.0 REQUIRED)

# spdlog
find_package(spdlog REQUIRED)

# Guidelines Support Library
find_package(Microsoft.GSL CONFIG)

# ACTS cmake-lint: disable=C0103
find_package(Acts REQUIRED COMPONENTS Core PluginDD4hep PluginJson)
set(Acts_VERSION_MIN "19.0.0")
set(Acts_VERSION
    "${Acts_VERSION_MAJOR}.${Acts_VERSION_MINOR}.${Acts_VERSION_PATCH}")
if(${Acts_VERSION} VERSION_LESS ${Acts_VERSION_MIN} AND NOT "${Acts_VERSION}"
                                                        STREQUAL "9.9.9")
  message(
    FATAL_ERROR
      "Acts version ${Acts_VERSION_MIN} or higher required, but ${Acts_VERSION} found"
  )
endif()
# ActsCore
get_target_property(ACTS_COMPILE_FEATURES ActsCore INTERFACE_COMPILE_FEATURES)
if(NOT "cxx_std_${CMAKE_CXX_STANDARD}" IN_LIST ACTS_COMPILE_FEATURES)
  message(
    WARNING
      "Acts has not been built with the requested C++ standard ${CMAKE_CXX_STANDARD}."
  )
endif()

# ROOT
find_package(ROOT REQUIRED)
# Check that ROOT is compiled with a modern enough c++ standard
get_target_property(ROOT_COMPILE_FEATURES ROOT::Core INTERFACE_COMPILE_FEATURES)
if(NOT "cxx_std_${CMAKE_CXX_STANDARD}" IN_LIST ROOT_COMPILE_FEATURES)
  message(
    WARNING
      "ROOT has not been built with the requested C++ standard ${CMAKE_CXX_STANDARD}."
  )
endif()
include(${ROOT_USE_FILE})

# Include my dirs
include_directories(${CMAKE_SOURCE_DIR}/src/include ${ROOT_INCLUDE_DIR})
file(GLOB sources ${CMAKE_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_SOURCE_DIR}/src/include/*.h)

#add binary 
add_executable(displacedVertex displacedvert.cc ${sources} ${headers})

#link executable 
#target_link_libraries(displacedVertex "-L${ROOT_LIBRARY_DIR}" ${ROOT_LIBRARIES})
target_link_libraries(displacedVertex PUBLIC ActsCore
      ActsExamplesFramework
      ActsExamplesIoJson
      ActsExamplesMagneticField
      "-L${ROOT_LIBRARY_DIR}" ${ROOT_LIBRARIES} spdlog::spdlog)
