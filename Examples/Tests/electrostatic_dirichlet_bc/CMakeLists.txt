# Add tests (alphabetical order) ##############################################
#

if(WarpX_PRECISION STREQUAL "DOUBLE")
    foreach(geometry IN ITEMS axis annular)
        set(radial_options "")
        if(geometry STREQUAL "annular")
            set(radial_options "--annular")
        endif()
        add_warpx_test(
            test_rz_${geometry}_dirichlet_picmi
            RZ
            2
            "inputs_test_rz_radial_dirichlet_picmi.py ${radial_options}"
            OFF
            OFF
            OFF
        )
    endforeach()
endif()

add_warpx_test(
    test_2d_dirichlet_bc  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_dirichlet_bc  # inputs
    "analysis.py"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)

add_warpx_test(
    test_2d_dirichlet_bc_picmi  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_dirichlet_bc_picmi.py  # inputs
    "analysis.py"  # analysis
    "analysis_default_regression.py --path diags/diag1000100"  # checksum
    OFF  # dependency
)
