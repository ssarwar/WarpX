# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_rz_rigid_source_recoil_picmi RZ 1
    "inputs_test_rz_rigid_source_picmi.py --solver Yee --subcycles 2"
    OFF OFF OFF
)

foreach(solver IN ITEMS Yee semi_implicit_em semi_implicit_mm)
    foreach(subcycles IN ITEMS 1 4)
        add_warpx_test(
            test_1d_background_mcc_time_${solver}_${subcycles}_picmi 1 2
            "inputs_test_1d_background_mcc_time_picmi.py --solver ${solver} --subcycles ${subcycles}"
            OFF OFF OFF
        )
    endforeach()
    add_warpx_test(
        test_1d_background_mcc_time_${solver}_early_picmi 1 2
        "inputs_test_1d_background_mcc_time_picmi.py --solver ${solver} --subcycles 4 --early"
        OFF OFF OFF
    )
endforeach()

foreach(solver IN ITEMS Yee PSATD semi_implicit_em semi_implicit_mm)
    if(solver STREQUAL "PSATD" AND NOT WarpX_FFT)
        continue()
    endif()
    add_warpx_test(
        test_rz_rigid_source_${solver}_picmi RZ 2
        "inputs_test_rz_rigid_source_picmi.py --solver ${solver} --subcycles 2"
        "analysis_rz_prescribed_outputs.py --solver ${solver}" OFF OFF
    )
    add_warpx_test(
        test_rz_rigid_source_${solver}_restart_picmi RZ 1
        "inputs_test_rz_rigid_source_picmi.py --solver ${solver} --subcycles 2 --restart ../test_rz_rigid_source_${solver}_picmi/diags/chk000002"
        "analysis_rz_prescribed_outputs.py --solver ${solver}" OFF "test_rz_rigid_source_${solver}_picmi"
    )
    if(solver STREQUAL "Yee" OR solver STREQUAL "PSATD")
        add_warpx_test(
            test_rz_rigid_source_${solver}_load_balance_picmi RZ 2
            "inputs_test_rz_rigid_source_picmi.py --solver ${solver} --subcycles 2 --load-balance"
            "analysis_rz_prescribed_outputs.py --solver ${solver}" OFF "test_rz_rigid_source_${solver}_picmi"
        )
    endif()
endforeach()

if(WarpX_OPENPMD)
    add_warpx_test(
        test_rz_prescribed_openpmd_picmi RZ 2
        "inputs_test_rz_rigid_source_picmi.py --solver Yee --subcycles 2 --openpmd"
        "analysis_rz_prescribed_outputs.py --openpmd" OFF OFF
    )
    add_warpx_test(
        test_rz_prescribed_openpmd_restart_picmi RZ 1
        "inputs_test_rz_rigid_source_picmi.py --solver Yee --subcycles 2 --openpmd --restart ../test_rz_prescribed_openpmd_picmi/diags/chk000002"
        "analysis_rz_prescribed_outputs.py --openpmd" OFF "test_rz_prescribed_openpmd_picmi"
    )
endif()

add_warpx_test(
    test_rz_prescribed_restart_guards_picmi RZ 1
    inputs_test_rz_prescribed_restart_guards_picmi.py
    OFF OFF "test_rz_rigid_source_Yee_picmi"
)

foreach(kind IN ITEMS proton ionization attachment)
    foreach(shape RANGE 1 4)
        add_warpx_test(
            test_rz_immobile_${kind}_${shape}_picmi
            RZ 2
            "inputs_test_rz_immobile_products_picmi.py --kind ${kind} --shape ${shape} --walls"
            OFF OFF OFF
        )
    endforeach()
    if(WarpX_FFT)
        add_warpx_test(
            test_rz_immobile_${kind}_psatd_picmi
            RZ 2
            "inputs_test_rz_immobile_products_picmi.py --kind ${kind} --solver PSATD"
            OFF OFF OFF
        )
    endif()
    if(WarpX_COMPUTE STREQUAL "OMP")
        add_warpx_test(
            test_rz_immobile_${kind}_threaded_picmi RZ 2
            "inputs_test_rz_immobile_products_picmi.py --kind ${kind} --shape 4 --walls --small-tiles"
            OFF OFF OFF
        )
        set_tests_properties(test_rz_immobile_${kind}_threaded_picmi.run PROPERTIES
            ENVIRONMENT_MODIFICATION "OMP_NUM_THREADS=set:4"
        )
    endif()
endforeach()

add_warpx_test(
    test_1d_background_mcc_attachment_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_attachment_picmi.py  # inputs
    "analysis_background_mcc_attachment.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_bad_inputs_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_bad_inputs_picmi.py  # inputs
    OFF  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_iaa_elastic_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_iaa_elastic_picmi.py  # inputs
    "analysis_background_mcc_iaa_elastic.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_loose_majorant_picmi
    1
    1
    inputs_test_1d_background_mcc_loose_majorant_picmi.py
    "analysis_background_mcc_selector.py"
    OFF
    OFF
)

add_warpx_test(
    test_1d_background_mcc_many_attachment_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_many_attachment_picmi.py  # inputs
    "analysis_background_mcc_many_attachment.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_many_excitation_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_many_excitation_picmi.py  # inputs
    "analysis_background_mcc_many_excitation.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_nonrelativistic_majorant_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_nonrelativistic_majorant_picmi.py  # inputs
    "analysis_background_mcc_nonrelativistic_majorant.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_relativistic  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_relativistic  # inputs
    "analysis_background_mcc_relativistic.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_rbeq_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_rbeq_picmi.py  # inputs
    "analysis_background_mcc_rbeq.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_selector_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_selector_picmi.py  # inputs
    "analysis_background_mcc_selector.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_background_mcc_rbeq_snapshot_picmi
    1 1
    "inputs_test_1d_background_mcc_rbeq_picmi.py --rbeq-model elmolcs_b8643810"
    "analysis_background_mcc_rbeq.py"
    OFF OFF
)

add_warpx_test(
    test_1d_background_mcc_rotation_picmi
    1 1
    inputs_test_1d_background_mcc_rotation_picmi.py
    "analysis_background_mcc_rotation.py"
    OFF OFF
)

add_warpx_test(
    test_1d_background_mcc_rotation_cumulative_picmi
    1 1
    "inputs_test_1d_background_mcc_rotation_picmi.py --cumulative"
    "analysis_background_mcc_rotation.py"
    OFF OFF
)

add_warpx_test(
    test_1d_background_mcc_thermal_attachment_picmi  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_background_mcc_thermal_attachment_picmi.py  # inputs
    "analysis_background_mcc_thermal_attachment.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_collision_z  # name
    1  # dims
    2  # nprocs
    inputs_test_1d_collision_z  # inputs
    "analysis_collision_1d.py diags/diag1000600"  # analysis
    "analysis_default_regression.py --path diags/diag1000600"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_collision_z_modulus_shuffle  # name
    1  # dims
    2  # nprocs
    inputs_test_1d_collision_z_modulus_shuffle  # inputs
    "analysis_collision_1d.py diags/diag1000600"  # analysis
    "analysis_default_regression.py --path diags/diag1000600"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_collision_z_correct_conservation  # name
    1  # dims
    2  # nprocs
    inputs_test_1d_collision_z_correct_conservation  # inputs
    "analysis_collision_1d_correct_conservation.py diags/diag1000000 diags/diag1000010"  # analysis
    "analysis_default_regression.py --path diags/diag1000010"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_collision_xz  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_collision_xz  # inputs
    "analysis_collision_2d.py diags/diag1000150"  # analysis
    "analysis_default_regression.py --path diags/diag1000150"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_collision_xz_global_debye  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_collision_xz_global_debye  # inputs
    "analysis_collision_2d.py diags/diag1000150"  # analysis
    "analysis_default_regression.py --path diags/diag1000150"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_collision_xz_picmi  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_collision_xz_picmi.py  # inputs
    "analysis_collision_2d.py diags/diag1000150"  # analysis
    "analysis_default_regression.py --path diags/diag1000150"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_collision_pulsed_decay  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_collision_pulsed_decay  # inputs
    "analysis_collision_3d_pulsed_decay.py diags/"  # analysis
    "analysis_default_regression.py --path diags/diag1025000"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_collision_iso  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_collision_iso  # inputs
    "analysis_collision_3d_isotropization.py diags/diag1000100"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_collision_iso_modulus_shuffle  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_collision_iso_modulus_shuffle  # inputs
    "analysis_collision_3d_isotropization.py diags/diag1000100"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_collision_iso_subcycle  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_collision_iso_subcycle  # inputs
    "analysis_collision_3d_isotropization.py diags/diag1000010"  # analysis
    "analysis_default_regression.py --path diags/diag1000010"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_collision_xyz  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_collision_xyz  # inputs
    "analysis_collision_3d.py diags/diag1000100"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_collision_angles_picmi  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_collision_angles_picmi.py  # inputs
    OFF  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_proton_impact_ionization_picmi  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_proton_impact_ionization_picmi.py  # inputs
    "analysis_proton_impact_ionization.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_proton_impact_ionization_performance_picmi  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_proton_impact_ionization_performance_picmi.py  # inputs
    "analysis_proton_impact_ionization_performance.py"  # analysis
    OFF  # checksum
    OFF  # dependency
)

foreach(case IN ITEMS cold low_energy mixed)
    add_warpx_test(
        test_3d_proton_impact_ionization_${case}_picmi
        3
        1
        inputs_test_3d_proton_impact_ionization_${case}_picmi.py
        "analysis_proton_impact_ionization.py"
        OFF
        OFF
    )
endforeach()

foreach(case IN ITEMS bad_inputs budget)
    add_warpx_test(
        test_3d_proton_impact_ionization_${case}_picmi
        3
        1
        inputs_test_3d_proton_impact_ionization_${case}_picmi.py
        OFF
        OFF
        OFF
    )
endforeach()

add_warpx_test(
    test_rz_proton_impact_ionization_budget_picmi
    RZ
    1
    "inputs_test_3d_proton_impact_ionization_budget_picmi.py --rz"
    OFF
    OFF
    OFF
)

add_warpx_test(
    test_rz_rigid_source_bad_inputs_picmi
    RZ
    1
    "inputs_test_3d_proton_impact_ionization_bad_inputs_picmi.py --rz"
    OFF
    OFF
    OFF
)

add_warpx_test(
    test_3d_proton_impact_ionization_checkpoint_picmi
    3
    1
    "inputs_test_3d_proton_impact_ionization_budget_picmi.py --checkpoint"
    OFF
    OFF
    OFF
)

add_warpx_test(
    test_3d_proton_impact_ionization_restart_picmi
    3
    1
    "inputs_test_3d_proton_impact_ionization_budget_picmi.py --restart ../test_3d_proton_impact_ionization_checkpoint_picmi/diags/chk000002 --reference ../test_3d_proton_impact_ionization_checkpoint_picmi/proton_impact_ionization_budget_results.npz"
    OFF
    OFF
    test_3d_proton_impact_ionization_checkpoint_picmi
)

add_warpx_test(
    test_3d_proton_impact_ionization_o2_performance_picmi
    3
    1
    inputs_test_3d_proton_impact_ionization_o2_performance_picmi.py
    "analysis_proton_impact_ionization_performance.py"
    OFF
    OFF
)

add_warpx_test(
    test_rz_collision  # name
    RZ  # dims
    1  # nprocs
    inputs_test_rz_collision  # inputs
    "analysis_collision_rz.py diags/diag1000150"  # analysis
    "analysis_default_regression.py --path diags/diag1000150 --skip-particles"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_collision_z_Bremsstrahlung  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_collision_z_Bremsstrahlung  # inputs
    "analysis_collision_1d_Bremsstrahlung.py diags/diag1000100"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_1d_collision_z_inverse_Bremsstrahlung  # name
    1  # dims
    1  # nprocs
    inputs_test_1d_collision_z_inverse_Bremsstrahlung  # inputs
    "analysis_collision_1d_Bremsstrahlung.py diags/diag1000100"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_charge_exchange_dsmc  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_charge_exchange_dsmc  # inputs
    "analysis_test_2d_collisions_split_momentum_push.py --path diags/reducedfiles/"  # analysis
    "analysis_default_regression.py --path diags/diag1/"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_collisions_split_momentum_push_electrostatic  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_collisions_split_momentum_push_electrostatic  # inputs
    "analysis_test_2d_collisions_split_momentum_push.py --path diags/reducedfiles/"  # analysis
    "analysis_default_regression.py --path diags/diag1/"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_collisions_split_momentum_push_electromagnetic  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_collisions_split_momentum_push_electromagnetic  # inputs
    "analysis_test_2d_collisions_split_momentum_push.py --path diags/reducedfiles/"  # analysis
    "analysis_default_regression.py --path diags/diag1/"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_collisions_split_momentum_push_electromagnetic_vay  # name
    2  # dims
    1  # nprocs
    inputs_test_2d_collisions_split_momentum_push_electromagnetic_vay  # inputs
    "analysis_test_2d_collisions_split_momentum_push.py --path diags/reducedfiles/"  # analysis
    "analysis_default_regression.py --path diags/diag1/"  # checksum
    OFF  # dependency
)
foreach(solver IN ITEMS Yee PSATD semi_implicit_em semi_implicit_mm)
    if(solver STREQUAL "PSATD" AND NOT WarpX_FFT)
        continue()
    endif()
    add_warpx_test(
        test_rz_rigid_train_${solver}_picmi RZ 2
        "inputs_test_rz_rigid_train_picmi.py --solver ${solver}"
        OFF OFF OFF
    )
    add_warpx_test(
        test_rz_rigid_train_${solver}_restart_picmi RZ 1
        "inputs_test_rz_rigid_train_picmi.py --solver ${solver} --restart ../test_rz_rigid_train_${solver}_picmi/diags/chk000075"
        OFF OFF "test_rz_rigid_train_${solver}_picmi"
    )
endforeach()
foreach(solver IN ITEMS Yee PSATD semi_implicit_em semi_implicit_mm)
    if(solver STREQUAL "PSATD" AND NOT WarpX_FFT)
        continue()
    endif()
    add_warpx_test(
        test_rz_attachment_checkpoint_${solver}_picmi RZ 2
        "inputs_test_rz_immobile_products_picmi.py --kind attachment --solver ${solver}"
        OFF OFF OFF
    )
    add_warpx_test(
        test_rz_attachment_checkpoint_${solver}_restart_picmi RZ 1
        "inputs_test_rz_immobile_products_picmi.py --kind attachment --solver ${solver} --restart ../test_rz_attachment_checkpoint_${solver}_picmi/diags/chk000002"
        OFF OFF "test_rz_attachment_checkpoint_${solver}_picmi"
    )
endforeach()
add_warpx_test(
    test_rz_rigid_source_alpha_picmi RZ 2
    "inputs_test_rz_rigid_source_picmi.py --alpha --subcycles 2"
    OFF OFF OFF
)
