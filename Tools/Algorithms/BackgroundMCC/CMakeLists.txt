cmake_minimum_required(VERSION 3.25)
project(BackgroundMCCPhysics LANGUAGES C CXX)

# Exercise the same host/device kernels with an existing AMReX package.
find_package(AMReX REQUIRED CONFIG)
if(AMReX_CUDA_FOUND)
    enable_language(CUDA)
    set_source_files_properties(test_kinematics.cpp test_interpolation.cpp
                                PROPERTIES LANGUAGE CUDA)
elseif(AMReX_HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(test_kinematics.cpp test_interpolation.cpp
                                PROPERTIES LANGUAGE HIP)
endif()

enable_testing()
foreach(test IN ITEMS kinematics interpolation)
    add_executable(test_mcc_${test} test_${test}.cpp)
    target_compile_features(test_mcc_${test} PRIVATE cxx_std_17)
    target_include_directories(test_mcc_${test} PRIVATE ../../../Source)
    target_link_libraries(test_mcc_${test} PRIVATE AMReX::amrex)
    add_test(NAME mcc_${test} COMMAND test_mcc_${test})
    set_tests_properties(mcc_${test} PROPERTIES TIMEOUT 30)
endforeach()
