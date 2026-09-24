cmake_minimum_required(VERSION 3.24)
project(WarpXCrossSections LANGUAGES CXX)
add_executable(export_rbeq export_rbeq.cpp)
target_compile_features(export_rbeq PRIVATE cxx_std_17)
target_include_directories(export_rbeq PRIVATE ../../Source)

if(TARGET lib_1d)
    find_package(Python 3.8 COMPONENTS Interpreter REQUIRED)
    add_executable(test_mcc_rotation test_mcc_rotation.cpp)
    target_link_libraries(test_mcc_rotation PRIVATE lib_1d)
    target_include_directories(test_mcc_rotation PRIVATE ../../Source)
    add_executable(test_secondary_angles test_secondary_angles.cpp)
    target_link_libraries(test_secondary_angles PRIVATE lib_1d)
    target_include_directories(test_secondary_angles PRIVATE ../../Source)
    if(WarpX_COMPUTE STREQUAL CUDA)
        setup_target_for_cuda_compilation(test_mcc_rotation)
        setup_target_for_cuda_compilation(test_secondary_angles)
    endif()
    add_test(NAME test_1d_secondary_angles COMMAND test_secondary_angles)
    add_test(NAME test_1d_rotation_reference
        COMMAND ${Python_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../../Examples/Tests/collision/analysis_rotation_reference.py
            ${CMAKE_CURRENT_BINARY_DIR}/rotation_reference)
    set_tests_properties(test_1d_rotation_reference PROPERTIES FIXTURES_SETUP rotation_reference)
    add_test(NAME test_1d_rotation_sampler
        COMMAND test_mcc_rotation
            reference=${CMAKE_CURRENT_BINARY_DIR}/rotation_reference/reference.txt
            output=${CMAKE_CURRENT_BINARY_DIR}/rotation_reference/sampler.txt samples=32768)
    set_tests_properties(test_1d_rotation_sampler PROPERTIES FIXTURES_REQUIRED rotation_reference)
    add_test(NAME test_1d_rotation_cumulative_reference
        COMMAND test_mcc_rotation cumulative=1
            reference=${CMAKE_CURRENT_BINARY_DIR}/rotation_reference/reference.txt
            output=${CMAKE_CURRENT_BINARY_DIR}/rotation_reference/cumulative.txt samples=32768)
    set_tests_properties(test_1d_rotation_cumulative_reference PROPERTIES FIXTURES_REQUIRED rotation_reference)
    add_test(NAME test_1d_rbeq_source_quadrature
        COMMAND ${Python_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../../Examples/Tests/collision/analysis_rbeq_sources.py
            $<TARGET_FILE:export_rbeq>)
endif()
