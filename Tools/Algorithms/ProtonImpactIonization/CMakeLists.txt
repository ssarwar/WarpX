cmake_minimum_required(VERSION 3.25)
project(ProtonImpactIonizationPhysics LANGUAGES C CXX)

# Use an existing AMReX package; do not configure another WarpX build.
find_package(AMReX REQUIRED CONFIG)
if(AMReX_CUDA_FOUND)
    enable_language(CUDA)
    set_source_files_properties(test_relativistic.cpp PROPERTIES LANGUAGE CUDA)
elseif(AMReX_HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(test_relativistic.cpp PROPERTIES LANGUAGE HIP)
endif()

add_executable(test_proton_ionization_relativistic test_relativistic.cpp)
target_compile_features(test_proton_ionization_relativistic PRIVATE cxx_std_17)
target_include_directories(test_proton_ionization_relativistic PRIVATE ../../..)
target_link_libraries(test_proton_ionization_relativistic PRIVATE AMReX::amrex)

enable_testing()
add_test(NAME proton_ionization_relativistic COMMAND test_proton_ionization_relativistic)
set_tests_properties(proton_ionization_relativistic PROPERTIES TIMEOUT 30)
