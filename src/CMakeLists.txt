cmake_minimum_required(VERSION 3.20)



#The configure file from the old autotools build...
configure_file(config.h.in config.h)
# The per-sequence scale beta is an ordinary model parameter that the objective
# uses ("together" in the old --with-beta-optimize terms). ObjFunc.cpp refuses
# to build without it.
add_compile_definitions(BETAOPTTOGETHER)

#required packages
find_package(GSL REQUIRED)
find_package(NLOPT REQUIRED)

# Warnings for our own code (nlopt was added above and keeps its own flags).
# The two -Werror= classes have each hidden a real bug in this code base:
# a switch written with plain labels instead of case labels, and a stray
# "std:" label.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
                        -Werror=unused-label -Werror=return-type)
endif()



#
############## LIBRARIES ############
#

#imported package in our own source-tree
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/snot")

add_subdirectory("utils")

add_subdirectory("tools")

add_subdirectory("annotation")

add_subdirectory("aspects")

add_library(gem_everything DataSet.cpp ExprFunc.cpp ExprFunc_rates.cpp ExprModel.cpp ExprPar.cpp ExprPredictor.cpp FactorIntFunc.cpp IO.cpp ObjFunc.cpp PredictorTrainer.cpp)

# Sequences are predicted independently of each other inside every objective
# evaluation; with OpenMP they are predicted in parallel (see --threads).
option(GEMSTAT_OPENMP "Predict sequences in parallel with OpenMP" ON)
if(GEMSTAT_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(gem_everything OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP found: sequences will be predicted in parallel")
    else()
        message(STATUS "OpenMP not found: single-threaded build")
    endif()
endif()
target_link_libraries(gem_everything GS_tools GS_annot)
target_link_libraries(gem_everything snot)
target_link_libraries(gem_everything "nlopt")#should get includes
target_include_directories(gem_everything PUBLIC ${NLOPT_INCLUDE_DIRS})

#Signalling extensions
add_library(GS_signalling DataSet_signal.cpp DataSet_signal.h)
target_link_libraries(GS_signalling GS_tools GS_annot)

#
############# BINARIES ####################
#
# add the executables
add_executable(seq2expr "seq2expr.cpp")
target_link_libraries(seq2expr main_aspects)
target_link_libraries(seq2expr gem_everything GS_annot GS_tools)
target_link_libraries(seq2expr GS_signalling)
target_link_libraries(seq2expr ${NLOPT_LIBRARIES} "nlopt")
target_link_libraries(seq2expr ${GSL_LIBRARIES})
