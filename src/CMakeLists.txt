cmake_minimum_required(VERSION 3.20)



#The configure file from the old autotools build...
configure_file(config.h.in config.h)

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


add_library(gem_everything DataSet.cpp ExprFunc.cpp ExprModel.cpp ExprPar.cpp ExprPredictor.cpp FactorIntFunc.cpp IO.cpp ObjFunc.cpp PredictorTrainer.cpp)

if(GEMSTAT_BETA_OPTIMIZE STREQUAL "together")
    target_compile_definitions(gem_everything PUBLIC BETAOPTTOGETHER)
elseif(GEMSTAT_BETA_OPTIMIZE STREQUAL "separate")
    # nothing defined: the objective functions solve for beta themselves
else()
    message(FATAL_ERROR "GEMSTAT_BETA_OPTIMIZE must be 'together' or 'separate', got '${GEMSTAT_BETA_OPTIMIZE}'")
endif()
message(STATUS "GEMSTAT_BETA_OPTIMIZE = ${GEMSTAT_BETA_OPTIMIZE}")
target_link_libraries(gem_everything GS_tools GS_annot)
target_link_libraries(gem_everything snot)
target_link_libraries(gem_everything "nlopt")#should get includes
target_include_directories(gem_everything PUBLIC ${NLOPT_INCLUDE_DIRS})


#
############# BINARIES ####################
#
# add the executables
add_executable(seq2expr seq2expr.cpp)
target_link_libraries(seq2expr gem_everything GS_annot GS_tools)
target_link_libraries(seq2expr ${NLOPT_LIBRARIES} "nlopt")
target_link_libraries(seq2expr ${GSL_LIBRARIES})
