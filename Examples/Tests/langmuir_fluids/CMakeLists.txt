# Add tests (alphabetical order) ##############################################
#

foreach(solver IN ITEMS Yee PSATD semi_implicit_em semi_implicit_mm)
    if(solver STREQUAL "PSATD" AND NOT WarpX_FFT)
        continue()
    endif()
    foreach(shape RANGE 1 4)
        add_warpx_test(
            test_rz_rigid_beam_${solver}_${shape}_picmi RZ 2
            "inputs_test_rz_rigid_beam_picmi.py --solver ${solver} --shape ${shape}"
            OFF OFF OFF
        )
    endforeach()
    add_warpx_test(
        test_rz_rigid_beam_${solver}_self_fields_picmi RZ 2
        "inputs_test_rz_rigid_beam_picmi.py --solver ${solver} --negative --self-fields"
        OFF OFF OFF
    )
endforeach()

foreach(solver IN ITEMS Yee PSATD)
    if(solver STREQUAL "PSATD" AND NOT WarpX_FFT)
        continue()
    endif()
    add_warpx_test(
        test_rz_immobile_fluid_${solver}_picmi
        RZ
        2
        "inputs_test_rz_immobile_fluid_picmi.py --solver ${solver}"
        OFF OFF OFF
    )
    add_warpx_test(
        test_rz_immobile_fluid_${solver}_restart_picmi
        RZ
        1
        "inputs_test_rz_immobile_fluid_picmi.py --solver ${solver} --restart ../test_rz_immobile_fluid_${solver}_picmi/diags/chk000002"
        OFF OFF
        "test_rz_immobile_fluid_${solver}_picmi"
    )
endforeach()

add_warpx_test(
    test_1d_langmuir_fluid  # name
    1  # dims
    2  # nprocs
    inputs_test_1d_langmuir_fluid  # inputs
    "analysis_1d.py diags/diag1000080"  # analysis
    "analysis_default_regression.py --path diags/diag1000080"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_langmuir_fluid  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_langmuir_fluid  # inputs
    "analysis_2d.py diags/diag1000080"  # analysis
    "analysis_default_regression.py --path diags/diag1000080"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_langmuir_fluid  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_langmuir_fluid  # inputs
    "analysis_3d.py diags/diag1000040"  # analysis
    "analysis_default_regression.py --path diags/diag1000040"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_rz_langmuir_fluid  # name
    RZ  # dims
    2  # nprocs
    inputs_test_rz_langmuir_fluid  # inputs
    "analysis_rz.py diags/diag1000080"  # analysis
    "analysis_default_regression.py --path diags/diag1000080"  # checksum
    OFF  # dependency
)
