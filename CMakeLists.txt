cmake_minimum_required(VERSION 3.20)


option(DOWNLOAD_BOOST "Use this to download boost into this directory, instead of searching the system." OFF)

# How the per-sequence expression scale factor ("beta") is treated during
# training.  This mirrors --with-beta-optimize from the retired autotools
# build, whose default was "together".
#   together : beta is an ordinary model parameter.  The objective function
#              uses the stored beta, so the value written to the parameter
#              file is the one the predictions were scored with.
#   separate : the objective solves for the best beta analytically on every
#              evaluation and discards it; the stored beta is never fitted.
#              (Predictions written by -fo are still multiplied by the
#              stored beta, so it does not match the reported objective.)
set(GEMSTAT_BETA_OPTIMIZE "together" CACHE STRING
    "How beta is optimized during training: together (default) or separate")
set_property(CACHE GEMSTAT_BETA_OPTIMIZE PROPERTY STRINGS together separate)


########## CMAKE SETUP ###########
#For building on much newer CMAKE
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif(POLICY CMP0079)

if(POLICY CMP0022)
    cmake_policy(SET CMP0022 NEW)
endif(POLICY CMP0022)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    )
include("PreventInTreeBuilds")



#
######### PROJECT ###########
# set the project name
project(GEMSTAT VERSION 1.9.0 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


add_subdirectory("src")

enable_testing()
add_subdirectory("tests")
