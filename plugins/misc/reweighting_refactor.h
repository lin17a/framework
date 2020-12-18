// -*- C++ -*-
// author: lina meyer
// short: for reweighting SM ttbar events

#include "misc/function_util.h"
#include "misc/constants.h"
#include "misc/numeric_vector.h"i

#include <complex>
#include "TLorentzVector.h"


template <typename Number = float>
const std::vector<std::pair<std::string, Number>>& 
compute_spin_correlation(Number pTop_pt, Number pTop_eta, Number pTop_phi, Number pTop_m,
                         Number aTop_pt, Number aTop_eta, Number aTop_phi, Number aTop_m,
                         Number pLep_pt, Number pLep_eta, Number pLep_phi, Number pLep_m,
                         Number aLep_pt, Number aLep_eta, Number aLep_phi, Number aLep_m)
{
  using namespace Framework;
  static std::vector<std::pair<std::string, Number>> m_spin_corr;
  static int initialize = 0;
  if (initialize == 0) {
    m_spin_corr.reserve(100); // exact size
    m_spin_corr.emplace_back("cLab", -9999.);
    m_spin_corr.emplace_back("dPhi", -9999.);
    m_spin_corr.emplace_back("dEta", -9999.);

    m_spin_corr.emplace_back("cpTP", -9999.);

    m_spin_corr.emplace_back("kdx", -9999.);
    m_spin_corr.emplace_back("kdy", -9999.);
    m_spin_corr.emplace_back("kdz", -9999.);

    m_spin_corr.emplace_back("rdx", -9999.);
    m_spin_corr.emplace_back("rdy", -9999.);
    m_spin_corr.emplace_back("rdz", -9999.);

    m_spin_corr.emplace_back("ndx", -9999.);
    m_spin_corr.emplace_back("ndy", -9999.);
    m_spin_corr.emplace_back("ndz", -9999.);

    m_spin_corr.emplace_back("b1k", -9999.);
    m_spin_corr.emplace_back("b2k", -9999.);

    m_spin_corr.emplace_back("b1j", -9999.);
    m_spin_corr.emplace_back("b2j", -9999.);

    m_spin_corr.emplace_back("b1r", -9999.);
    m_spin_corr.emplace_back("b2r", -9999.);

    m_spin_corr.emplace_back("b1q", -9999.);
    m_spin_corr.emplace_back("b2q", -9999.);

    m_spin_corr.emplace_back("b1n", -9999.);
    m_spin_corr.emplace_back("b2n", -9999.);

    m_spin_corr.emplace_back("b1x", -9999.);
    m_spin_corr.emplace_back("b2x", -9999.);

    m_spin_corr.emplace_back("b1y", -9999.);
    m_spin_corr.emplace_back("b2y", -9999.);

    m_spin_corr.emplace_back("b1z", -9999.);
    m_spin_corr.emplace_back("b2z", -9999.);

    m_spin_corr.emplace_back("bPkk", -9999.);
    m_spin_corr.emplace_back("bMkk", -9999.);

    m_spin_corr.emplace_back("bPjj", -9999.);
    m_spin_corr.emplace_back("bMjj", -9999.);

    m_spin_corr.emplace_back("bPrr", -9999.);
    m_spin_corr.emplace_back("bMrr", -9999.);

    m_spin_corr.emplace_back("bPqq", -9999.);
    m_spin_corr.emplace_back("bMqq", -9999.);

    m_spin_corr.emplace_back("bPnn", -9999.);
    m_spin_corr.emplace_back("bMnn", -9999.);

    m_spin_corr.emplace_back("bPxx", -9999.);
    m_spin_corr.emplace_back("bMxx", -9999.);

    m_spin_corr.emplace_back("bPyy", -9999.);
    m_spin_corr.emplace_back("bMyy", -9999.);

    m_spin_corr.emplace_back("bPzz", -9999.);
    m_spin_corr.emplace_back("bMzz", -9999.);

    m_spin_corr.emplace_back("ckk", -9999.);
    m_spin_corr.emplace_back("crr", -9999.);
    m_spin_corr.emplace_back("cnn", -9999.);

    m_spin_corr.emplace_back("crk", -9999.);
    m_spin_corr.emplace_back("ckr", -9999.);

    m_spin_corr.emplace_back("cnr", -9999.);
    m_spin_corr.emplace_back("crn", -9999.);

    m_spin_corr.emplace_back("cnk", -9999.);
    m_spin_corr.emplace_back("ckn", -9999.);

    m_spin_corr.emplace_back("cPrk", -9999.);
    m_spin_corr.emplace_back("cMrk", -9999.);

    m_spin_corr.emplace_back("cPnr", -9999.);
    m_spin_corr.emplace_back("cMnr", -9999.);

    m_spin_corr.emplace_back("cPnk", -9999.);
    m_spin_corr.emplace_back("cMnk", -9999.);

    m_spin_corr.emplace_back("cxx", -9999.);
    m_spin_corr.emplace_back("cyy", -9999.);
    m_spin_corr.emplace_back("czz", -9999.);

    m_spin_corr.emplace_back("cyx", -9999.);
    m_spin_corr.emplace_back("cxy", -9999.);

    m_spin_corr.emplace_back("czy", -9999.);
    m_spin_corr.emplace_back("cyz", -9999.);

    m_spin_corr.emplace_back("czx", -9999.);
    m_spin_corr.emplace_back("cxz", -9999.);

    m_spin_corr.emplace_back("cPyx", -9999.);
    m_spin_corr.emplace_back("cMyx", -9999.);

    m_spin_corr.emplace_back("cPzy", -9999.);
    m_spin_corr.emplace_back("cMzy", -9999.);

    m_spin_corr.emplace_back("cPzx", -9999.);
    m_spin_corr.emplace_back("cMzx", -9999.);

    m_spin_corr.emplace_back("cHel", -9999.);

    m_spin_corr.emplace_back("cHan", -9999.);
    m_spin_corr.emplace_back("cSca", -9999.);
    m_spin_corr.emplace_back("cTra", -9999.);

    m_spin_corr.emplace_back("crkP", -9999.);
    m_spin_corr.emplace_back("cnrP", -9999.);
    m_spin_corr.emplace_back("cnkP", -9999.);

    m_spin_corr.emplace_back("crkM", -9999.);
    m_spin_corr.emplace_back("cnrM", -9999.);
    m_spin_corr.emplace_back("cnkM", -9999.);

    m_spin_corr.emplace_back("cXxx", -9999.);
    m_spin_corr.emplace_back("cYyy", -9999.);
    m_spin_corr.emplace_back("cZzz", -9999.);

    m_spin_corr.emplace_back("cyxP", -9999.);
    m_spin_corr.emplace_back("czyP", -9999.);
    m_spin_corr.emplace_back("czxP", -9999.);

    m_spin_corr.emplace_back("cyxM", -9999.);
    m_spin_corr.emplace_back("czyM", -9999.);
    m_spin_corr.emplace_back("czxM", -9999.);

    m_spin_corr.emplace_back("kNorm", -9999.);
    m_spin_corr.emplace_back("rNorm", -9999.);
    m_spin_corr.emplace_back("nNorm", -9999.);

    m_spin_corr.emplace_back("phi0", -9999.);
    m_spin_corr.emplace_back("phi1", -9999.);

    m_spin_corr.emplace_back("cpTTT", -9999.);

    ++initialize;
  }

  static std::array<Number, 16> arg;
  if (arg[0]  == pTop_pt and arg[1]  == pTop_eta and arg[2]  == pTop_phi and arg[3]  == pTop_m and 
      arg[4]  == aTop_pt and arg[5]  == aTop_eta and arg[6]  == aTop_phi and arg[7]  == aTop_m and 
      arg[8]  == pLep_pt and arg[9]  == pLep_eta and arg[10] == pLep_phi and arg[11] == pLep_m and 
      arg[12] == aLep_pt and arg[13] == aLep_eta and arg[14] == aLep_phi and arg[15] == aLep_m)
    return m_spin_corr;

  arg[0]  = pTop_pt; arg[1]  = pTop_eta; arg[2]  = pTop_phi; arg[3]  = pTop_m;
  arg[4]  = aTop_pt; arg[5]  = aTop_eta; arg[6]  = aTop_phi; arg[7]  = aTop_m;
  arg[8]  = pLep_pt; arg[9]  = pLep_eta; arg[10] = pLep_phi; arg[11] = pLep_m;
  arg[12] = aLep_pt; arg[13] = aLep_eta; arg[14] = aLep_phi; arg[15] = aLep_m;

  static TLorentzVector p4lab_pTop, p4lab_aTop, p4lab_pLep, p4lab_aLep;
  p4lab_pTop.SetPtEtaPhiM(arg[0],  arg[1],  arg[2],  arg[3]);
  p4lab_aTop.SetPtEtaPhiM(arg[4],  arg[5],  arg[6],  arg[7]);
  p4lab_pLep.SetPtEtaPhiM(arg[8],  arg[9],  arg[10], arg[11]);
  p4lab_aLep.SetPtEtaPhiM(arg[12], arg[13], arg[14], arg[15]);

  // computes spin correlation variables
  static const int icLab = index_with_key(m_spin_corr, "cLab");
  m_spin_corr[icLab].second = p4lab_aLep.Vect().Unit().Dot( p4lab_pLep.Vect().Unit() );

  static const int idPhi = index_with_key(m_spin_corr, "dPhi");
  m_spin_corr[idPhi].second = dPhi(p4lab_aLep.Phi(), p4lab_pLep.Phi());

  static const int idEta = index_with_key(m_spin_corr, "dEta");
  m_spin_corr[idEta].second = absolute_difference(p4lab_aLep.Eta(), p4lab_pLep.Eta());

  // the various spin corr vars with boosting: Bernreuther 1508.05271
  // xyz coordinate system
  static const TVector3 xBase(1., 0., 0.);
  static const TVector3 yBase(0., 1., 0.);
  static const TVector3 zBase(0., 0., 1.);

  // the necessary boosting such that lepton ~ top spin vector
  const TLorentzVector p4lab_TT( p4lab_pTop + p4lab_aTop );

  const auto f_zmf_tt = [&p4lab_TT] (TLorentzVector p4) {
    p4.Boost( -1. * p4lab_TT.BoostVector() );
    return p4;
  };
  const TLorentzVector p4hel_pTop = f_zmf_tt(p4lab_pTop);
  const TLorentzVector p4hel_aTop = f_zmf_tt(p4lab_aTop);

  const auto f_zmf_top = [&f_zmf_tt] (const TLorentzVector &lep, const TLorentzVector &top) {
    auto p4 = f_zmf_tt(lep);
    p4.Boost( -1. * top.BoostVector() );
    return p4;
  };
  const TLorentzVector p4hel_aLep = f_zmf_top(p4lab_aLep, p4hel_pTop);
  const TLorentzVector p4hel_pLep = f_zmf_top(p4lab_pLep, p4hel_aTop);

  // calculating the top-beam angle for pTop only (defining pP as the +z beam proton)
  static const int icpTP = index_with_key(m_spin_corr, "cpTP");
  m_spin_corr[icpTP].second = p4hel_pTop.Vect().Unit().Dot(zBase);
  const Number spTP = std::sqrt(1. - (m_spin_corr[icpTP].second * m_spin_corr[icpTP].second));

  // the signs needed to account for Bose symmetry
  const Number sY = (m_spin_corr[icpTP].second >= 0.) ? 1. : -1.;
  const Number sD = ( std::abs(p4lab_pTop.Rapidity()) >= std::abs(p4lab_aTop.Rapidity()) ) ? 1. : -1.;

  // define the base vectors a
  // j and q base are the k* and r* respectively
  // b is always -a
  const TVector3 kBase = p4hel_pTop.Vect().Unit();
  const TVector3 jBase = sD * kBase;
  const TVector3 rBase = ( (sY / spTP) * (zBase - (m_spin_corr[icpTP].second * kBase)) ).Unit();
  const TVector3 qBase = sD * rBase;
  const TVector3 nBase = ( (sY / spTP) * zBase.Cross(kBase) ).Unit();

  // for resolution studies
  // naming eg ndx for n dot x -> x component of the n axis
  static const int ikdx = index_with_key(m_spin_corr, "kdx");
  static const int ikdy = index_with_key(m_spin_corr, "kdy");
  static const int ikdz = index_with_key(m_spin_corr, "kdz");

  static const int irdx = index_with_key(m_spin_corr, "rdx");
  static const int irdy = index_with_key(m_spin_corr, "rdy");
  static const int irdz = index_with_key(m_spin_corr, "rdz");

  static const int indx = index_with_key(m_spin_corr, "ndx");
  static const int indy = index_with_key(m_spin_corr, "ndy");
  static const int indz = index_with_key(m_spin_corr, "ndz");

  m_spin_corr[ikdx].second = kBase.x();
  m_spin_corr[ikdy].second = kBase.y();
  m_spin_corr[ikdz].second = kBase.z();

  m_spin_corr[irdx].second = rBase.x();
  m_spin_corr[irdy].second = rBase.y();
  m_spin_corr[irdz].second = rBase.z();

  m_spin_corr[indx].second = nBase.x();
  m_spin_corr[indy].second = nBase.y();
  m_spin_corr[indz].second = nBase.z();

  // find the polarization angles along bases
  static const int ib1k = index_with_key(m_spin_corr, "b1k");
  static const int ib2k = index_with_key(m_spin_corr, "b2k");

  static const int ib1j = index_with_key(m_spin_corr, "b1j");
  static const int ib2j = index_with_key(m_spin_corr, "b2j");

  static const int ib1r = index_with_key(m_spin_corr, "b1r");
  static const int ib2r = index_with_key(m_spin_corr, "b2r");

  static const int ib1q = index_with_key(m_spin_corr, "b1q");
  static const int ib2q = index_with_key(m_spin_corr, "b2q");

  static const int ib1n = index_with_key(m_spin_corr, "b1n");
  static const int ib2n = index_with_key(m_spin_corr, "b2n");

  static const int ib1x = index_with_key(m_spin_corr, "b1x");
  static const int ib2x = index_with_key(m_spin_corr, "b2x");

  static const int ib1y = index_with_key(m_spin_corr, "b1y");
  static const int ib2y = index_with_key(m_spin_corr, "b2y");

  static const int ib1z = index_with_key(m_spin_corr, "b1z");
  static const int ib2z = index_with_key(m_spin_corr, "b2z");

  m_spin_corr[ib1k].second = p4hel_aLep.Vect().Unit().Dot( kBase );
  m_spin_corr[ib2k].second = p4hel_pLep.Vect().Unit().Dot( -1. * kBase );

  m_spin_corr[ib1j].second = p4hel_aLep.Vect().Unit().Dot( jBase );
  m_spin_corr[ib2j].second = p4hel_pLep.Vect().Unit().Dot( -1. * jBase );

  m_spin_corr[ib1r].second = p4hel_aLep.Vect().Unit().Dot( rBase );
  m_spin_corr[ib2r].second = p4hel_pLep.Vect().Unit().Dot( -1. * rBase );

  m_spin_corr[ib1q].second = p4hel_aLep.Vect().Unit().Dot( qBase );
  m_spin_corr[ib2q].second = p4hel_pLep.Vect().Unit().Dot( -1. * qBase );

  m_spin_corr[ib1n].second = p4hel_aLep.Vect().Unit().Dot( nBase );
  m_spin_corr[ib2n].second = p4hel_pLep.Vect().Unit().Dot( -1. * nBase );

  m_spin_corr[ib1x].second = p4hel_aLep.Vect().Unit().Dot( xBase );
  m_spin_corr[ib2x].second = p4hel_pLep.Vect().Unit().Dot( -1. * xBase );

  m_spin_corr[ib1y].second = p4hel_aLep.Vect().Unit().Dot( yBase );
  m_spin_corr[ib2y].second = p4hel_pLep.Vect().Unit().Dot( -1. * yBase );

  m_spin_corr[ib1z].second = p4hel_aLep.Vect().Unit().Dot( zBase );
  m_spin_corr[ib2z].second = p4hel_pLep.Vect().Unit().Dot( -1. * zBase );

  // sums and differences of the polarization angles
  static const int ibPkk = index_with_key(m_spin_corr, "bPkk");
  static const int ibMkk = index_with_key(m_spin_corr, "bMkk");

  static const int ibPjj = index_with_key(m_spin_corr, "bPjj");
  static const int ibMjj = index_with_key(m_spin_corr, "bMjj");

  static const int ibPrr = index_with_key(m_spin_corr, "bPrr");
  static const int ibMrr = index_with_key(m_spin_corr, "bMrr");

  static const int ibPqq = index_with_key(m_spin_corr, "bPqq");
  static const int ibMqq = index_with_key(m_spin_corr, "bMqq");

  static const int ibPnn = index_with_key(m_spin_corr, "bPnn");
  static const int ibMnn = index_with_key(m_spin_corr, "bMnn");

  static const int ibPxx = index_with_key(m_spin_corr, "bPxx");
  static const int ibMxx = index_with_key(m_spin_corr, "bMxx");

  static const int ibPyy = index_with_key(m_spin_corr, "bPyy");
  static const int ibMyy = index_with_key(m_spin_corr, "bMyy");

  static const int ibPzz = index_with_key(m_spin_corr, "bPzz");
  static const int ibMzz = index_with_key(m_spin_corr, "bMzz");

  m_spin_corr[ibPkk].second = m_spin_corr[ib1k].second + m_spin_corr[ib2k].second;
  m_spin_corr[ibMkk].second = m_spin_corr[ib1k].second - m_spin_corr[ib2k].second;

  m_spin_corr[ibPjj].second = m_spin_corr[ib1j].second + m_spin_corr[ib2j].second;
  m_spin_corr[ibMjj].second = m_spin_corr[ib1j].second - m_spin_corr[ib2j].second;

  m_spin_corr[ibPrr].second = m_spin_corr[ib1r].second + m_spin_corr[ib2r].second;
  m_spin_corr[ibMrr].second = m_spin_corr[ib1r].second - m_spin_corr[ib2r].second;

  m_spin_corr[ibPqq].second = m_spin_corr[ib1q].second + m_spin_corr[ib2q].second;
  m_spin_corr[ibMqq].second = m_spin_corr[ib1q].second - m_spin_corr[ib2q].second;

  m_spin_corr[ibPnn].second = m_spin_corr[ib1n].second + m_spin_corr[ib2n].second;
  m_spin_corr[ibMnn].second = m_spin_corr[ib1n].second - m_spin_corr[ib2n].second;

  m_spin_corr[ibPxx].second = m_spin_corr[ib1x].second + m_spin_corr[ib2x].second;
  m_spin_corr[ibMxx].second = m_spin_corr[ib1x].second - m_spin_corr[ib2x].second;

  m_spin_corr[ibPyy].second = m_spin_corr[ib1y].second + m_spin_corr[ib2y].second;
  m_spin_corr[ibMyy].second = m_spin_corr[ib1y].second - m_spin_corr[ib2y].second;

  m_spin_corr[ibPzz].second = m_spin_corr[ib1z].second + m_spin_corr[ib2z].second;
  m_spin_corr[ibMzz].second = m_spin_corr[ib1z].second - m_spin_corr[ib2z].second;

  // Cab = -9<cab>
  static const int ickk = index_with_key(m_spin_corr, "ckk");
  static const int icrr = index_with_key(m_spin_corr, "crr");
  static const int icnn = index_with_key(m_spin_corr, "cnn");

  static const int icrk = index_with_key(m_spin_corr, "crk");
  static const int ickr = index_with_key(m_spin_corr, "ckr");

  static const int icnr = index_with_key(m_spin_corr, "cnr");
  static const int icrn = index_with_key(m_spin_corr, "crn");

  static const int icnk = index_with_key(m_spin_corr, "cnk");
  static const int ickn = index_with_key(m_spin_corr, "ckn");

  static const int icxx = index_with_key(m_spin_corr, "cxx");
  static const int icyy = index_with_key(m_spin_corr, "cyy");
  static const int iczz = index_with_key(m_spin_corr, "czz");

  static const int icyx = index_with_key(m_spin_corr, "cyx");
  static const int icxy = index_with_key(m_spin_corr, "cxy");

  static const int iczy = index_with_key(m_spin_corr, "czy");
  static const int icyz = index_with_key(m_spin_corr, "cyz");

  static const int iczx = index_with_key(m_spin_corr, "czx");
  static const int icxz = index_with_key(m_spin_corr, "cxz");

  m_spin_corr[ickk].second = m_spin_corr[ib1k].second * m_spin_corr[ib2k].second;
  m_spin_corr[icrr].second = m_spin_corr[ib1r].second * m_spin_corr[ib2r].second;
  m_spin_corr[icnn].second = m_spin_corr[ib1n].second * m_spin_corr[ib2n].second;

  m_spin_corr[icrk].second = m_spin_corr[ib1r].second * m_spin_corr[ib2k].second;
  m_spin_corr[ickr].second = m_spin_corr[ib1k].second * m_spin_corr[ib2r].second;

  m_spin_corr[icnr].second = m_spin_corr[ib1n].second * m_spin_corr[ib2r].second;
  m_spin_corr[icrn].second = m_spin_corr[ib1r].second * m_spin_corr[ib2n].second;

  m_spin_corr[icnk].second = m_spin_corr[ib1n].second * m_spin_corr[ib2k].second;
  m_spin_corr[ickn].second = m_spin_corr[ib1k].second * m_spin_corr[ib2n].second;

  m_spin_corr[icxx].second = m_spin_corr[ib1x].second * m_spin_corr[ib2x].second;
  m_spin_corr[icyy].second = m_spin_corr[ib1y].second * m_spin_corr[ib2y].second;
  m_spin_corr[iczz].second = m_spin_corr[ib1z].second * m_spin_corr[ib2z].second;

  m_spin_corr[icyx].second = m_spin_corr[ib1y].second * m_spin_corr[ib2x].second;
  m_spin_corr[icxy].second = m_spin_corr[ib1x].second * m_spin_corr[ib2y].second;

  m_spin_corr[iczy].second = m_spin_corr[ib1z].second * m_spin_corr[ib2y].second;
  m_spin_corr[icyz].second = m_spin_corr[ib1y].second * m_spin_corr[ib2z].second;

  m_spin_corr[iczx].second = m_spin_corr[ib1z].second * m_spin_corr[ib2x].second;
  m_spin_corr[icxz].second = m_spin_corr[ib1x].second * m_spin_corr[ib2z].second;

  // sums and differences of cij
  static const int icPrk = index_with_key(m_spin_corr, "cPrk");
  static const int icMrk = index_with_key(m_spin_corr, "cMrk");

  static const int icPnr = index_with_key(m_spin_corr, "cPnr");
  static const int icMnr = index_with_key(m_spin_corr, "cMnr");

  static const int icPnk = index_with_key(m_spin_corr, "cPnk");
  static const int icMnk = index_with_key(m_spin_corr, "cMnk");

  static const int icPyx = index_with_key(m_spin_corr, "cPyx");
  static const int icMyx = index_with_key(m_spin_corr, "cMyx");

  static const int icPzy = index_with_key(m_spin_corr, "cPzy");
  static const int icMzy = index_with_key(m_spin_corr, "cMzy");

  static const int icPzx = index_with_key(m_spin_corr, "cPzx");
  static const int icMzx = index_with_key(m_spin_corr, "cMzx");

  m_spin_corr[icPrk].second = m_spin_corr[icrk].second + m_spin_corr[ickr].second;
  m_spin_corr[icMrk].second = m_spin_corr[icrk].second - m_spin_corr[ickr].second;

  m_spin_corr[icPnr].second = m_spin_corr[icnr].second + m_spin_corr[icrn].second;
  m_spin_corr[icMnr].second = m_spin_corr[icnr].second - m_spin_corr[icrn].second;

  m_spin_corr[icPnk].second = m_spin_corr[icnk].second + m_spin_corr[ickn].second;
  m_spin_corr[icMnk].second = m_spin_corr[icnk].second - m_spin_corr[ickn].second;

  m_spin_corr[icPyx].second = m_spin_corr[icyx].second + m_spin_corr[icxy].second;
  m_spin_corr[icMyx].second = m_spin_corr[icyx].second - m_spin_corr[icxy].second;

  m_spin_corr[icPzy].second = m_spin_corr[iczy].second + m_spin_corr[icyz].second;
  m_spin_corr[icMzy].second = m_spin_corr[iczy].second - m_spin_corr[icyz].second;

  m_spin_corr[icPzx].second = m_spin_corr[iczx].second + m_spin_corr[icxz].second;
  m_spin_corr[icMzx].second = m_spin_corr[iczx].second - m_spin_corr[icxz].second;

  // opening angle between the spin vectors
  static const int icHel = index_with_key(m_spin_corr, "cHel");
  m_spin_corr[icHel].second = p4hel_aLep.Vect().Unit().Dot( p4hel_pLep.Vect().Unit() );

  // cHel with one spin vector flipped about a given axis
  // useful to measure Cii = 1.5(Di - D)
  // where Di is the coeff for these flipped cHel
  static const int icHan = index_with_key(m_spin_corr, "cHan");
  static const int icSca = index_with_key(m_spin_corr, "cSca");
  static const int icTra = index_with_key(m_spin_corr, "cTra");

  m_spin_corr[icHan].second = 0. + m_spin_corr[ickk].second - m_spin_corr[icrr].second - m_spin_corr[icnn].second;
  m_spin_corr[icSca].second = 0. - m_spin_corr[ickk].second + m_spin_corr[icrr].second - m_spin_corr[icnn].second;
  m_spin_corr[icTra].second = 0. - m_spin_corr[ickk].second - m_spin_corr[icrr].second + m_spin_corr[icnn].second;

  // cHel with one of the spin vectors having two of its components swapped
  // in other words, flipped about the bisector of the i and j axes
  // useful to measure Cij + Cji = -3Dij - Ckk for axes i, j and k
  static const int icrkP = index_with_key(m_spin_corr, "crkP");
  static const int icnrP = index_with_key(m_spin_corr, "cnrP");
  static const int icnkP = index_with_key(m_spin_corr, "cnkP");

  m_spin_corr[icnrP].second = 0. - m_spin_corr[ickk].second - m_spin_corr[icnr].second - m_spin_corr[icrn].second;
  m_spin_corr[icnkP].second = 0. - m_spin_corr[icnk].second - m_spin_corr[icrr].second - m_spin_corr[ickn].second;
  m_spin_corr[icrkP].second = 0. - m_spin_corr[icrk].second - m_spin_corr[ickr].second - m_spin_corr[icnn].second;

  // the above naturally lends itself to a measurement of Cij - Cji
  // with swapping and also sign flipping
  static const int icrkM = index_with_key(m_spin_corr, "crkM");
  static const int icnrM = index_with_key(m_spin_corr, "cnrM");
  static const int icnkM = index_with_key(m_spin_corr, "cnkM");

  m_spin_corr[icnrM].second = 0. - m_spin_corr[ickk].second - m_spin_corr[icnr].second + m_spin_corr[icrn].second;
  m_spin_corr[icnkM].second = 0. - m_spin_corr[icnk].second - m_spin_corr[icrr].second + m_spin_corr[ickn].second;
  m_spin_corr[icrkM].second = 0. - m_spin_corr[icrk].second + m_spin_corr[ickr].second - m_spin_corr[icnn].second;

  // the flipped cHel in the xyz system
  static const int icXxx = index_with_key(m_spin_corr, "cXxx");
  static const int icYyy = index_with_key(m_spin_corr, "cYyy");
  static const int icZzz = index_with_key(m_spin_corr, "cZzz");

  m_spin_corr[icXxx].second = 0. + m_spin_corr[icxx].second - m_spin_corr[icyy].second - m_spin_corr[iczz].second;
  m_spin_corr[icYyy].second = 0. - m_spin_corr[icxx].second + m_spin_corr[icyy].second - m_spin_corr[iczz].second;
  m_spin_corr[icZzz].second = 0. - m_spin_corr[icxx].second - m_spin_corr[icyy].second + m_spin_corr[iczz].second;

  static const int icyxP = index_with_key(m_spin_corr, "cyxP");
  static const int iczyP = index_with_key(m_spin_corr, "czyP");
  static const int iczxP = index_with_key(m_spin_corr, "czxP");

  m_spin_corr[iczyP].second = 0. - m_spin_corr[icxx].second - m_spin_corr[iczy].second - m_spin_corr[icyz].second;
  m_spin_corr[iczxP].second = 0. - m_spin_corr[iczx].second - m_spin_corr[icyy].second - m_spin_corr[icxz].second;
  m_spin_corr[icyxP].second = 0. - m_spin_corr[icyx].second - m_spin_corr[icxy].second - m_spin_corr[iczz].second;

  static const int icyxM = index_with_key(m_spin_corr, "cyxM");
  static const int iczyM = index_with_key(m_spin_corr, "czyM");
  static const int iczxM = index_with_key(m_spin_corr, "czxM");

  m_spin_corr[iczyM].second = 0. - m_spin_corr[icxx].second - m_spin_corr[iczy].second + m_spin_corr[icyz].second;
  m_spin_corr[iczxM].second = 0. - m_spin_corr[iczx].second - m_spin_corr[icyy].second + m_spin_corr[icxz].second;
  m_spin_corr[icyxM].second = 0. - m_spin_corr[icyx].second + m_spin_corr[icxy].second - m_spin_corr[iczz].second;

  // these are the O_CP1 and O_CP2 as in page 18 (why not O_CP3 with n base too)
  static const int ikNorm = index_with_key(m_spin_corr, "kNorm");
  static const int irNorm = index_with_key(m_spin_corr, "rNorm");
  static const int inNorm = index_with_key(m_spin_corr, "nNorm");

  const TVector3 llNorm = p4hel_aLep.Vect().Unit().Cross( p4hel_pLep.Vect().Unit() );
  m_spin_corr[ikNorm].second =  llNorm.Dot( kBase );
  m_spin_corr[irNorm].second =  llNorm.Dot( rBase );
  m_spin_corr[inNorm].second =  llNorm.Dot( nBase );

  // angles as in 1702.06063; phi0 = phi* and phi1 = phi*_CP
  // sensitive to CP mixtures (shows up as a phase in the distributions)
  static const int iphi0 = index_with_key(m_spin_corr, "phi0");
  static const int iphi1 = index_with_key(m_spin_corr, "phi1");

  const TVector3 t3_aLep = ( p4hel_aLep.Vect().Unit() - (m_spin_corr[ib1k].second * kBase) ).Unit();
  const TVector3 t3_pLep = ( p4hel_pLep.Vect().Unit() - (-1. * m_spin_corr[ib2k].second * kBase) ).Unit();

  m_spin_corr[iphi0].second = std::acos(t3_aLep.Dot(t3_pLep));
  m_spin_corr[iphi1].second = (m_spin_corr[ikNorm].second < 0.) ? (2. * constants::pi<Number>) - m_spin_corr[iphi0].second : m_spin_corr[iphi0].second;

  // cos theta top lep used in lj 2016 A/H
  // labeled by top rather than final state branch
  // other side not needed, pTop and aTop are back to back
  static const int icpTTT = index_with_key(m_spin_corr, "cpTTT");
  m_spin_corr[icpTTT].second = p4hel_pTop.Vect().Unit().Dot( p4lab_TT.Vect().Unit() );

  return m_spin_corr;
}



template <typename Number = float>
auto spin_correlation(const std::string &var)
{
  const auto &map = compute_spin_correlation(Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number());

  return [ivar = index_with_key(map, var)] (Number pTop_pt, Number pTop_eta, Number pTop_phi, Number pTop_m,
                                            Number aTop_pt, Number aTop_eta, Number aTop_phi, Number aTop_m,
                                            Number pLep_pt, Number pLep_eta, Number pLep_phi, Number pLep_m,
                                            Number aLep_pt, Number aLep_eta, Number aLep_phi, Number aLep_m)
  {
    const auto &map = compute_spin_correlation(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                                               aTop_pt, aTop_eta, aTop_phi, aTop_m,
                                               pLep_pt, pLep_eta, pLep_phi, pLep_m,
                                               aLep_pt, aLep_eta, aLep_phi, aLep_m);

    return map[ivar].second;
  };
}




double weight_ggHtt() {
      // Check that the chosen options are valid
      if (httInput.HIGGS_OPTION>2 || httInput.WIDTH_OPTION>3) {
            printf("What? HIGGS_OPTION=%d, WIDTH_OPTION=%d; returning weight=1 !!\n", httInput.HIGGS_OPTION, httInput.WIDTH_OPTION);
            return 1.;
      }

      // setting up variables with input parameters
      double MH2_REF = pow(httInput.MH,2);
      double MA2_REF = pow(httInput.MA,2);
      double ytop2 = pow(httInput.YTOP,2);

      TLorentzVector P4gen_t1 = httInput.P4gen_t[0]; 
      TLorentzVector P4gen_t2 = httInput.P4gen_t[1]; 
      TLorentzVector P4gen_d1 = httInput.P4gen_d[0]; 
      TLorentzVector P4gen_d2 = httInput.P4gen_d[1]; 

      // setting constants
      const double top_mass_ref = 172.5;
      const double top_mass_ref_sq = MT_REF*MT_REF;
      //Fermi-Konstante
      const double fermi_const = 1.16637876e-5;
      const double pi = TMath::Pi();
      const double pi_sq = pi*pi;
      const double 2_sqrt = sqrt(2.);

      
      TLorentzVector P4gen_tt = P4gen_t1+P4gen_t2;
      double s_sqrt = P4gen_tt.M();
      double s_abs = s_sqrt*s_sqrt;
      double top_mass = httInput.P4gen_t[0].M();
      double antitop_mass = httInput.P4gen_t[1].M();
      double beta_sq = (1-(top_mass-antitop_mass)*(top_mass-antitop_mass)/sqrts2)*(1-(top_mass+antitop_mass)*(top_mass+antitop_mass)/sqrts2);
      double beta_abs = sqrt(beta2);
      double beta_3 = beta*beta2;
      double mtrefOverE = 2*MT_REF/s_sqrt;
      double beta2Ref = 1-mtrefOverE*mtrefOverE;

      // zboost vector in -z direction, with z-part of the ttbar impulse
      TVector3 zBoost(0.,0.,-P4gen_tt.Pz()/P4gen_tt.E()); 
      P4gen_tt.Boost(zBoost);
      // point reflections ttbar vector
      TVector3 transverseBoost(-P4gen_tt.Px()/P4gen_tt.E(),-P4gen_tt.Py()/P4gen_tt.E(),0.); 
      P4gen_tt.Boost(transverseBoost);
      // do -z-boost and transverse boost for all other lorentz vectors
      P4gen_t1.Boost(zBoost); 
      P4gen_t1.Boost(transverseBoost);
      P4gen_d1.Boost(zBoost); 
      P4gen_d1.Boost(transverseBoost);
      P4gen_t2.Boost(zBoost); 
      P4gen_t2.Boost(transverseBoost);
      P4gen_d2.Boost(zBoost); 
      P4gen_d2.Boost(transverseBoost);
      // z1 = top CosTheta
      double z1 = P4gen_t1.CosTheta();
      double z2 = z1*z1;

      // rotate things
      double minusPhi1Plane = -P4gen_t1.Phi();
      TVector3 zVect(0.,0.,1.);
      P4gen_t1.Rotate(minusPhi1Plane,zVect);
      P4gen_d1.Rotate(minusPhi1Plane,zVect);
      P4gen_t2.Rotate(minusPhi1Plane,zVect);
      P4gen_d2.Rotate(minusPhi1Plane,zVect);
      double t1Yangle = -atan2(P4gen_t1.Px(),P4gen_t1.Pz());
      TVector3 yVect(0.,1.,0.);
      P4gen_t1.Rotate(t1Yangle,yVect);
      P4gen_d1.Rotate(t1Yangle,yVect);
      P4gen_t2.Rotate(t1Yangle,yVect);
      P4gen_d2.Rotate(t1Yangle,yVect);

      TVector3 t1Boost(-P4gen_t1.Px()/P4gen_t1.E(),-P4gen_t1.Py()/P4gen_t1.E(),-P4gen_t1.Pz()/P4gen_t1.E()); 
      P4gen_d1.Boost(t1Boost);
      TLorentzVector P4gen_d1_trest = P4gen_d1;
      TVector3 t2Boost(-P4gen_t2.Px()/P4gen_t2.E(),-P4gen_t2.Py()/P4gen_t2.E(),-P4gen_t2.Pz()/P4gen_t2.E()); 
      P4gen_d2.Boost(t2Boost);
      TLorentzVector P4gen_d2_trest = P4gen_d2;

      double common_factor_for_M2_gg_QCD = pi2/12*(7+9*beta2*z2)/pow(1-beta2*z2,2);
      double common_factor_for_M2_gg_ss_nointerf = pi2/12*(5+9*beta2*z2)/pow(1-beta2*z2,2);
      double common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta2*z2);

      const std::complex<double> i_cmplx = 1i;

      std::complex<double> common_factor_for_M_H = 0.;
      double mh_gh_partial = 0.;
      double mh_gh = 0.;
      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
            if (httInput.MH>2*MT_REF) {
                  double beta3_tdecay = pow(1-4*MT2_REF/MH2_REF,1.5);
                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
                        mh_gh_partial += 3*GF*MT2_REF*sqrts2/(4*pi*sqrt2)*beta3_tdecay*ytop2;
                  } else {
                        mh_gh_partial += 3*GF*MT2_REF*MH2_REF/(4*pi*sqrt2)*beta3_tdecay*ytop2;
                  }
            }

            if (httInput.WIDTH_OPTION==1) {
                  mh_gh = httInput.GH*sqrts2/httInput.MH;
            } else if (httInput.WIDTH_OPTION==0) {
                  mh_gh = httInput.GH*httInput.MH;
            } else {
                  mh_gh = mh_gh_partial;
            }

            double NB_real, NB_imag;
            double auxh = 1.5*(1-beta2Ref)*ytop2;
            if (beta2Ref>0) {
                  double betaRef = sqrt(beta2Ref);
                  NB_real = auxh*(1-beta2Ref/4*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2));
                  NB_imag = auxh*(beta2Ref/2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  NB_real = auxh*(1+beta2Ref*pow(asin(s_sqrt/2/MT_REF),2));
                  NB_imag = 0.;
            }
            std::complex<double> NB = NB_real + NB_imag*1i;
            std::complex<double> denomH = sqrts2 - MH2_REF + mh_gh*1i;

            common_factor_for_M_H = NB / denomH;
            common_factor_for_M_H *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
      }

      std::complex<double> common_factor_for_M_A = 0.;
      double ma_ga_partial = 0.;
      double ma_ga = 0.;
      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
            if (httInput.MA>2*MT_REF) {
                  double beta_tdecay = pow(1-4*MT2_REF/MA2_REF,0.5);
                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
                        ma_ga_partial += 3*GF*MT2_REF*sqrts2/(4*pi*sqrt2)*beta_tdecay*ytop2;
                  } else {
                        ma_ga_partial += 3*GF*MT2_REF*MA2_REF/(4*pi*sqrt2)*beta_tdecay*ytop2;
                  }
            }
            if (httInput.WIDTH_OPTION==1) {
                        ma_ga = httInput.GA*sqrts2/httInput.MA;
            } else if (httInput.WIDTH_OPTION==0) {
                        ma_ga = httInput.GA*httInput.MA;
            } else {
                        ma_ga = ma_ga_partial;
            }

            double PB_real, PB_imag;
            double auxa = -1.5*(1-beta2Ref)/4.*ytop2;
            if (beta2Ref>0) {
                  double betaRef = sqrt(beta2Ref);
                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2);
                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  PB_real = -4*auxa*pow(asin(s_sqrt/2/MT_REF),2);
                  PB_imag = 0.;
            }
            std::complex<double> PB = PB_real + PB_imag*1i;
            std::complex<double> denomA = sqrts2 - MA2_REF + ma_ga*1i;

            common_factor_for_M_A = PB / denomA;
            common_factor_for_M_A *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
      }

      // Printing global input information (once)
      static bool printOnce = true;
      if (printOnce) {
            printOnce = false;
            printf("\nHTT_INPUT >>>> (printed only once)\n");
            printf("\tYTOP=%.3f times the SM Higgs-ttbar coupling\n", httInput.YTOP);
            if (httInput.HIGGS_OPTION==0) {
                  printf("\tOnly SCALAR CASE, MH=%.3f GeV\n", httInput.MH);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed width: GH=%.3f GeV\n", httInput.GH);
                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning width: GH(MH)=%.3f GeV\n", httInput.GH);
                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/sqrts2*httInput.MH);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed width: GH(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", mh_gh/sqrts2*httInput.MH);
                  }
            } else if (httInput.HIGGS_OPTION==1) {
                  printf("\tOnly PSEUDO-SCALAR CASE, MA=%.3f GeV\n", httInput.MA);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed width: GA=%.3f GeV\n", httInput.GA);
                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning width: GA(MA)=%.3f GeV\n", httInput.GA);
                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/sqrts2*httInput.MA);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed width: GA(top+bottom)=%.3f GeV\n", ma_ga/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", ma_ga/sqrts2*httInput.MA);
                  }
            } else {
                  printf("\tSCALAR + PSEUDO-SCALAR CASE, MH=%.3f GeV, MA=%.3f GeV\n", httInput.MH, httInput.MA);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed widths: GH=%.3f GeV, GA=%.3f GeV\n", httInput.GH, httInput.GA);
                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH, ma_ga_partial/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning widths: GH(MH)=%.3f GeV, GA(MA)=%.3f GeV\n", httInput.GH, httInput.GA);
                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/sqrts2*httInput.MH, ma_ga_partial/sqrts2*httInput.MA);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH, ma_ga/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/sqrts2*httInput.MH, ma_ga/sqrts2*httInput.MA);
                  }
            }
            printf("\n");
      }

      double c1 = P4gen_d1_trest.CosTheta();
      double c2 = P4gen_d2_trest.CosTheta();
      double s1 = sqrt(1-c1*c1);
      double s2 = sqrt(1-c2*c2);
      double phi1 = P4gen_d1_trest.Phi();
      double phi2 = P4gen_d2_trest.Phi();
      double os_factor = (1+z2+(1-z2)*(1-beta2Ref))*(1.-c1*c2) + 2*(1-beta2Ref)*(1-z2)*c1*c2
            - beta2Ref*(1-z2)*s1*s2*cos(phi1+phi2)
            -2*(1-beta2Ref)*(1-z2)*s1*s2*cos(phi1)*cos(phi2)
            +2*sqrt(1-beta2Ref)*z1*sqrt(1-z2)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
      double ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
      double ss_factor_beta = beta2 * (1+c1*c2 - s1*s2*cos(phi1-phi2));

      double M2_QCD = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
                    + common_factor_for_M2_gg_QCD*(1-beta2Ref)*(ss_factor_1+ss_factor_beta);

      double M2 = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
                + common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
                + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
                + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;

      // Final weight
      double weight = M2 / M2_QCD;
      if (isnan(weight)) {
            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
            return 0.;
      }

      return weight;
};

using attr_func_type_ttbar = const std::function<float(float,float,float,float,float,float,float,float)>;
using attr_func_type_spin = const std::function<float(float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float)>;

attr_func_type_ttbar get_attr_function(const std::function<float(TLorentzVector)> extract_attr_func) {
  return [extract_attr_func](float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    TLorentzVector sum = p1 + p2;
    return extract_attr_func(sum);
  };
}

const auto add_attribute = [&](const std::string &funcname, attr_func_type_ttbar func) {
    gen_ttbar.add_attribute(funcname, func, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                                            "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  };




auto reweight(const std::string &var)
{
  const auto &map = compute_weight(float, float, float, float,
                                   float, float, float, float,
                                   float, float, float, float,
                                   float, float, float, float,
                                   float, float, float, float,
                                   float, float, float, float,
                                   int, int);
  return float
  {
    const auto &map = compute_weight(                           aTop_pt, aTop_eta, aTop_phi, aTop_m,
                                               pLep_pt, pLep_eta, pLep_phi, pLep_m,
                                               aLep_pt, aLep_eta, aLep_phi, aLep_m);

    return weight;
  };
}



bool is_gg_initial_state(const int& pdgId1, const int& pdgId2, const double& pz_hard_radiation) {
/*
 * pdgId1 is the the PDG id of first parton, along the +Z direction
 * pdgId2 is the the PDG id of other parton, along the -Z direction
 * pz_hard_radiation is the pz of the system made of the hard radiated partons recoiling 
 *    against the ttbar system, measured in the LAB system
*/
      
      // Ultra-simple assignment ito gg or qqbar for qg initial states
      if (pdgId1==21 && pdgId2==21) {
            // gg
            return true;
      } else if (abs(pdgId1)<6 && pdgId2==-pdgId1) {
            // qqbar
            return false;
      } else if (pdgId1==21 || pdgId2==21) {
            // qg
            if ((pdgId1!=21 && pz_hard_radiation>0) || (pdgId2!=21 && pz_hard_radiation<0)) {
                  // radiated quark more collinear with initial quark: assume it is gg
                  return true;
            } else {
                  // radiated quark more collinear with initial gluon: assume it is qqbar
                  return false;
            }
      } 

      // If we got here then it is a rare initial state with multiple hard radiated partons: assume gg
      return true;
}
