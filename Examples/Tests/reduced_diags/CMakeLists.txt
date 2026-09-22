# Add tests (alphabetical order) ##############################################
#

set(fluid_flux_solvers Yee semi_implicit_em semi_implicit_mm)
if(WarpX_FFT)
    list(APPEND fluid_flux_solvers PSATD)
endif()
foreach(solver IN LISTS fluid_flux_solvers)
    add_warpx_test(
        test_rz_fluid_flux_${solver}_picmi
        RZ
        2
        "inputs_test_rz_fluid_diagnostics_picmi.py --nonzero-flux --solver ${solver}"
        OFF
        OFF
        OFF
    )
    add_warpx_test(
        test_rz_fluid_flux_${solver}_restart_picmi
        RZ
        1
        "inputs_test_rz_fluid_diagnostics_picmi.py --nonzero-flux --solver ${solver} --continue-output --restart ../test_rz_fluid_flux_${solver}_picmi/diags/chk000002"
        "analysis_rz_fluid_diagnostics.py ../test_rz_fluid_flux_${solver}_picmi . --nonzero-flux"
        OFF
        test_rz_fluid_flux_${solver}_picmi
    )
endforeach()

add_warpx_test(
    test_rz_fluid_diagnostics_picmi
    RZ
    2
    inputs_test_rz_fluid_diagnostics_picmi.py
    OFF
    OFF
    OFF
)

add_warpx_test(
    test_rz_fluid_diagnostics_restart_picmi
    RZ
    1
    "inputs_test_rz_fluid_diagnostics_picmi.py --continue-output --restart ../test_rz_fluid_diagnostics_picmi/diags/chk000002"
    "analysis_rz_fluid_diagnostics.py ../test_rz_fluid_diagnostics_picmi ."
    OFF
    test_rz_fluid_diagnostics_picmi
)

add_warpx_test(
    test_rz_fluid_diagnostics_bad_restart_picmi
    RZ
    1
    "inputs_test_rz_fluid_diagnostics_bad_restart_picmi.py ../test_rz_fluid_diagnostics_picmi/diags/chk000002"
    OFF
    OFF
    test_rz_fluid_diagnostics_picmi
)

add_warpx_test(
    test_3d_reduced_diags  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_reduced_diags  # inputs
    "analysis_reduced_diags.py diags/diag1000200"  # analysis
    "analysis_default_regression.py --path diags/diag1000200"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_reduced_diags_load_balance_costs_heuristic  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_reduced_diags_load_balance_costs_heuristic  # inputs
    "analysis_reduced_diags_load_balance_costs.py diags/diag1000003"  # analysis
    "analysis_default_regression.py --path diags/diag1000003"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_reduced_diags_load_balance_costs_timers  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_reduced_diags_load_balance_costs_timers  # inputs
    "analysis_reduced_diags_load_balance_costs.py diags/diag1000003"  # analysis
    "analysis_default_regression.py --path diags/diag1000003"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_3d_reduced_diags_load_balance_costs_timers_picmi  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_reduced_diags_load_balance_costs_timers_picmi.py  # inputs
    "analysis_reduced_diags_load_balance_costs.py diags/diag1000003"  # analysis
    "analysis_default_regression.py --path diags/diag1000003"  # checksum
    OFF  # dependency
)

if(WarpX_FFT)
    add_warpx_test(
        test_3d_reduced_diags_load_balance_costs_timers_psatd  # name
        3  # dims
        2  # nprocs
        inputs_test_3d_reduced_diags_load_balance_costs_timers_psatd  # inputs
        "analysis_reduced_diags_load_balance_costs.py diags/diag1000003"  # analysis
        "analysis_default_regression.py --path diags/diag1000003"  # checksum
        OFF  # dependency
    )
endif()
