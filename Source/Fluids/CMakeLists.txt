foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        MultiFluidContainer.cpp
        PrescribedFluid.cpp
        RigidBeam.cpp
        FluidDiagnostics.cpp
        QdsmcParticleContainer.cpp
        WarpXFluidContainer.cpp
    )
endforeach()
