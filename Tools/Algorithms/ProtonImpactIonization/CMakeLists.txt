cmake_minimum_required(VERSION 3.25)
project(ProtonImpactIonizationPhysics LANGUAGES C CXX)

# Use an existing AMReX package; do not configure another WarpX build.
find_package(AMReX REQUIRED CONFIG)
set(pjg_sources
    ../../../Source/Particles/Collision/ProtonImpactIonization/PJGModel.cpp
    ../../../Source/ablastr/utils/TextMsg.cpp
    ../../../Source/ablastr/utils/text/StringUtils.cpp)
if(AMReX_CUDA_FOUND)
    enable_language(CUDA)
    set_source_files_properties(test_relativistic.cpp test_sampling.cpp benchmark_sampling.cpp
                                test_pjg_model.cpp benchmark_pjg_model.cpp ${pjg_sources}
                                PROPERTIES LANGUAGE CUDA)
elseif(AMReX_HIP_FOUND)
    enable_language(HIP)
    set_source_files_properties(test_relativistic.cpp test_sampling.cpp benchmark_sampling.cpp
                                test_pjg_model.cpp benchmark_pjg_model.cpp ${pjg_sources}
                                PROPERTIES LANGUAGE HIP)
endif()

add_executable(test_proton_ionization_relativistic test_relativistic.cpp)
add_executable(test_proton_ionization_sampling test_sampling.cpp)
add_executable(benchmark_proton_ionization_sampling benchmark_sampling.cpp)
add_executable(test_proton_ionization_model test_pjg_model.cpp ${pjg_sources})
add_executable(benchmark_proton_ionization_model benchmark_pjg_model.cpp ${pjg_sources})
foreach(target IN ITEMS test_proton_ionization_relativistic test_proton_ionization_sampling
                        benchmark_proton_ionization_sampling test_proton_ionization_model
                        benchmark_proton_ionization_model)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_include_directories(${target} PRIVATE ../../.. ../../../Source)
    target_link_libraries(${target} PRIVATE AMReX::amrex)
endforeach()

enable_testing()
add_test(NAME proton_ionization_relativistic COMMAND test_proton_ionization_relativistic)
add_test(NAME proton_ionization_sampling COMMAND test_proton_ionization_sampling)
add_test(NAME proton_ionization_model COMMAND test_proton_ionization_model)
set_tests_properties(proton_ionization_model PROPERTIES TIMEOUT 60)
set_tests_properties(proton_ionization_relativistic proton_ionization_sampling PROPERTIES TIMEOUT 30)
