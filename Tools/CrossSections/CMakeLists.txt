cmake_minimum_required(VERSION 3.24)
project(WarpXCrossSections LANGUAGES CXX)
add_executable(export_rbeq export_rbeq.cpp)
target_compile_features(export_rbeq PRIVATE cxx_std_17)
target_include_directories(export_rbeq PRIVATE ../../Source)
set(WarpX_SPECTATOR_TEST_DATA "" CACHE PATH "Prepared offline IAA spectator MCC data")
set(WarpX_RECIPROCAL_TEST_DATA "" CACHE PATH "warpx-data MCC_cross_sections with prepared V6 bundles")

if(TARGET lib_1d)
    find_package(Python 3.8 COMPONENTS Interpreter REQUIRED)
    add_executable(test_reciprocal_rotation test_reciprocal_rotation.cpp)
    target_link_libraries(test_reciprocal_rotation PRIVATE lib_1d)
    target_include_directories(test_reciprocal_rotation PRIVATE ../../Source)
    add_executable(benchmark_reciprocal_mcc benchmark_reciprocal_mcc.cpp)
    target_link_libraries(benchmark_reciprocal_mcc PRIVATE lib_1d)
    target_include_directories(benchmark_reciprocal_mcc PRIVATE ../../Source)
    if(WarpX_COMPUTE STREQUAL CUDA)
        setup_target_for_cuda_compilation(test_reciprocal_rotation)
        setup_target_for_cuda_compilation(benchmark_reciprocal_mcc)
    endif()
    if(WarpX_RECIPROCAL_TEST_DATA)
        set(reciprocal_output "${CMAKE_CURRENT_BINARY_DIR}/reciprocal_tests")
        file(MAKE_DIRECTORY "${reciprocal_output}")
        foreach(gas IN ITEMS N2 O2)
            set(reciprocal_bundle "${WarpX_RECIPROCAL_TEST_DATA}/${gas}/IAA/reciprocal_hybrid_300K/thermal_rotation.rot")
            add_test(NAME test_1d_reciprocal_sampler_${gas}
                COMMAND test_reciprocal_rotation
                    file=${reciprocal_bundle}
                    cross_section=${WarpX_RECIPROCAL_TEST_DATA}/${gas}/IAA/elastic.txt
                    output=${reciprocal_output}/${gas}.txt samples=32768)
            set_tests_properties(test_1d_reciprocal_sampler_${gas} PROPERTIES
                FIXTURES_SETUP reciprocal_sampler_${gas})
            add_test(NAME test_1d_reciprocal_moments_${gas}
                COMMAND ${Python_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/reciprocal_rotation/sampling_reference.py
                    --bundle ${reciprocal_bundle} --sampler ${reciprocal_output}/${gas}.txt
                    --samples 32768)
            set_tests_properties(test_1d_reciprocal_moments_${gas} PROPERTIES
                FIXTURES_REQUIRED reciprocal_sampler_${gas}
                ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
        endforeach()
    endif()
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
    add_executable(test_rotation_thresholds test_rotation_thresholds.cpp)
    target_link_libraries(test_rotation_thresholds PRIVATE lib_1d)
    target_include_directories(test_rotation_thresholds PRIVATE ../../Source)
    if(WarpX_COMPUTE STREQUAL CUDA)
        setup_target_for_cuda_compilation(test_rotation_thresholds)
    endif()
    add_test(NAME test_1d_rotation_thresholds COMMAND test_rotation_thresholds)
    foreach(check IN ITEMS angles recoil)
        add_executable(test_proton_${check} ../Algorithms/ProtonImpactIonization/test_${check}.cpp)
        target_link_libraries(test_proton_${check} PRIVATE lib_1d)
        target_include_directories(test_proton_${check} PRIVATE ../.. ../../Source)
        if(WarpX_COMPUTE STREQUAL CUDA)
            setup_target_for_cuda_compilation(test_proton_${check})
        endif()
        add_test(NAME test_1d_proton_${check} COMMAND test_proton_${check})
    endforeach()
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
    if(WarpX_SPECTATOR_TEST_DATA)
        set(spectator_reference "${CMAKE_CURRENT_BINARY_DIR}/spectator_reference")
        add_test(NAME test_1d_spectator_source_quadrature
            COMMAND ${CMAKE_COMMAND} -E env OPENBLAS_NUM_THREADS=1
                ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_thermiaa_spectator.py
                --data-dir ${WarpX_SPECTATOR_TEST_DATA}
                --bundle ${WarpX_SPECTATOR_TEST_DATA}/N2/IAA/thermal_spectator.rot
                --output ${spectator_reference})
        set_tests_properties(test_1d_spectator_source_quadrature PROPERTIES
            FIXTURES_SETUP spectator_reference)
        foreach(mode IN ITEMS 0 1)
            add_test(NAME test_1d_spectator_sampler_${mode}
                COMMAND test_mcc_rotation model=iaa_spectator cumulative=${mode}
                    reference=${spectator_reference}/reference.txt
                    output=${spectator_reference}/sampler-${mode}.txt samples=32768
                    n2_dcs=${WarpX_SPECTATOR_TEST_DATA}/N2/IAA/elastic_dcs.txt)
            set_tests_properties(test_1d_spectator_sampler_${mode} PROPERTIES
                FIXTURES_REQUIRED spectator_reference)
        endforeach()
    endif()
    add_test(NAME test_1d_rbeq_source_quadrature
        COMMAND ${Python_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../../Examples/Tests/collision/analysis_rbeq_sources.py
            $<TARGET_FILE:export_rbeq>)
endif()
