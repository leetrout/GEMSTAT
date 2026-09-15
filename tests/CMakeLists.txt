# Golden-output regression tests, driven by ctest.
#
# Every directory under tests/cases/ that contains a case.json becomes one
# ctest test named golden.<dir>.  See tests/run_case.py for the case format.
#
# To regenerate the expected outputs after an intentional behaviour change:
#   make golden-update      (or: cmake --build . --target golden-update)
# then review the diff in tests/cases/*/expected/ before committing.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(WARNING "Python3 not found: golden regression tests are disabled")
    return()
endif()

set(GEMSTAT_TEST_RUNNER ${CMAKE_CURRENT_SOURCE_DIR}/run_case.py)
set(GEMSTAT_TEST_COMMON_ARGS
    --binary $<TARGET_FILE:seq2expr>
    --root ${PROJECT_SOURCE_DIR}
    --data ${PROJECT_SOURCE_DIR}/data
    --testdata ${CMAKE_CURRENT_SOURCE_DIR}/data)

file(GLOB _case_dirs LIST_DIRECTORIES true ${CMAKE_CURRENT_SOURCE_DIR}/cases/*)
foreach(_case_dir ${_case_dirs})
    if(NOT EXISTS ${_case_dir}/case.json)
        continue()
    endif()
    get_filename_component(_name ${_case_dir} NAME)
    add_test(NAME golden.${_name}
             COMMAND ${Python3_EXECUTABLE} ${GEMSTAT_TEST_RUNNER} ${GEMSTAT_TEST_COMMON_ARGS}
                     --case ${_case_dir}
                     --workdir ${CMAKE_CURRENT_BINARY_DIR}/work/${_name})
    file(READ ${_case_dir}/case.json _spec)
    string(JSON _nlabels ERROR_VARIABLE _err LENGTH ${_spec} labels)
    if(NOT _err AND _nlabels GREATER 0)
        math(EXPR _last "${_nlabels} - 1")
        set(_labels "")
        foreach(_i RANGE ${_last})
            string(JSON _l GET ${_spec} labels ${_i})
            list(APPEND _labels ${_l})
        endforeach()
        set_tests_properties(golden.${_name} PROPERTIES LABELS "${_labels}")
    endif()
    string(JSON _timeout ERROR_VARIABLE _err GET ${_spec} timeout)
    if(NOT _err)
        set_tests_properties(golden.${_name} PROPERTIES TIMEOUT ${_timeout})
    endif()
endforeach()

add_custom_target(golden-update
    COMMAND ${Python3_EXECUTABLE} ${GEMSTAT_TEST_RUNNER} ${GEMSTAT_TEST_COMMON_ARGS}
            --all-cases ${CMAKE_CURRENT_SOURCE_DIR}/cases --update
    DEPENDS seq2expr
    COMMENT "Regenerating golden expected outputs from the current seq2expr")
