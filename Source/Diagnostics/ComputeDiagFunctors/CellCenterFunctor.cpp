#include "CellCenterFunctor.H"

#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

CellCenterFunctor::CellCenterFunctor(amrex::MultiFab const * mf_src, int lev,
                                     amrex::IntVect crse_ratio,
                                     bool convertRZmodes2cartesian, int ncomp, int scomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev),
      m_scomp(scomp),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{}

void
CellCenterFunctor::operator()(amrex::MultiFab& mf_dst, int dcomp, const int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_mf_src != nullptr, "m_mf_src can't be a nullptr.");
    AMREX_ASSUME(m_mf_src != nullptr);

    amrex::MultiFab source(*m_mf_src, amrex::make_alias, m_scomp,
        m_convertRZmodes2cartesian ? m_mf_src->nComp()-m_scomp : nComp());
    InterpolateMFForDiag(
        mf_dst, source, dcomp,
        warpx.DistributionMap(m_lev),m_convertRZmodes2cartesian);
}
