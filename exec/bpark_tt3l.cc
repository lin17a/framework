// execution macro for testing the fwk devs
// compile:
// g++ $(root-config --cflags --evelibs) -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I ../plugins/ -I ../src/ -o bpark_tt3l bpark_tt3l.cc

#include "../src/Dataset.h"
#include "../src/Collection.h"
#include "../src/Aggregate.h"
#include "../src/Histogram.h"
#include "../src/Tree.h"

#include "misc/string_io.h"
#include "misc/function_util.h"
#include "misc/numeric_vector.h"
#include "misc/constants.h"
#include "misc/input_dataset.h"

// http://tclap.sourceforge.net/
#include "tclap/CmdLine.h"

int main(int argc, char** argv) {
  TCLAP::CmdLine cmdbase("b parking tt -> 3l differential analysis, 2018 dataset", ' ', "0.01");
  //TCLAP::ValueArg<std::string> cmdchannel("", "channel","Analysis channel - ee, emu, mumu", false, "emu", "string", cmdbase);
  TCLAP::ValueArg<std::string> cmddataset("", "dataset","Dataset to be used -  consult the input_dataset.h plugin", false, 
                                          "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", "string", cmdbase);
  //TCLAP::SwitchArg cmdtest("", "test", "Dummy boolean argument", cmdbase, false);
  cmdbase.parse( argc, argv );

  // get the value parsed by each arg 
  //std::string channel = cmdchannel.getValue();
  std::string dataset = cmddataset.getValue();
  //bool test = cmdtest.getValue();

  using namespace Framework;

  Dataset<TChain> dat(dataset, "Events");
  dat.set_files(file_list(dataset), 2);

  std::vector<std::string> hlt;
  Collection<boolean, uint, unsigned long long> meta("meta", 43);
  meta.add_attribute("run", "run", 1U);
  meta.add_attribute("lumi", "luminosityBlock", 1U);
  meta.add_attribute("event", "event", 1ULL);
  for (uint ip = 0; ip < 5; ++ip) {
    hlt.emplace_back("HLT_Mu7_IP4_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu8_IP3_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu8_IP5_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu8_IP6_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu9_IP4_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu9_IP5_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu9_IP6_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);

    hlt.emplace_back("HLT_Mu12_IP6_part" + to_str(ip));
    meta.add_attribute(hlt.back(), hlt.back(), true);
  }
  //hlt.emplace_back("HLT_Ele32_WPTight_Gsf");
  //meta.add_attribute(hlt.back(), hlt.back(), true);

  Collection<int, float> gen_particle("gen_particle", "nGenPart", 10, 256);
  gen_particle.add_attribute("gpmass", "GenPart_mass", 1.f);
  gen_particle.add_attribute("pt", "GenPart_pt", 1.f);
  gen_particle.add_attribute("eta", "GenPart_eta", 1.f);
  gen_particle.add_attribute("phi", "GenPart_phi", 1.f);
  gen_particle.add_attribute("pdg", "GenPart_pdgId", 1);
  gen_particle.add_attribute("status", "GenPart_status", 1);
  gen_particle.add_attribute("flag", "GenPart_statusFlags", 1);
  gen_particle.add_attribute("mother", "GenPart_genPartIdxMother", 1);
  gen_particle.transform_attribute("dileptonic_ttbar", 
                                   [&pdgs = gen_particle.get<int>("pdg"), &idxs = gen_particle.get<int>("mother"), 
                                    &flags = gen_particle.get<int>("flag")] 
                                   (int pdg, int flag, int idx) -> int {
                                     // integer flag for particles which are part of a generator-level dileptonic ttbar system
                                     // 1 top
                                     // 2 antitop
                                     // 3 W+
                                     // 4 W-
                                     // 5 bottom (parton)
                                     // 6 antibottom (parton)
                                     // 7 W+ antilepton - e, mu
                                     // 8 W- lepton - e, mu
                                     // 9 W+ neutrino - e, mu 
                                     // 10 W- antineutrino - e, mu
                                     // 11 W+ antilepton - tau
                                     // 12 W- lepton - tau
                                     // 13 W+ neutrino - tau 
                                     // 14 W- antineutrino - tau

                                     // technically not dileptonic anymore but ok
                                     // 15 bottom muon
                                     // 16 antibottom muon

                                     // final top quarks
                                     if (std::abs(pdg) == 6 and flag & 8192)
                                       return (pdg > 0) ? 1 : 2;

                                     // W boson block
                                     // choose the final copy prior to decay like the top
                                     if (std::abs(pdg) == 24 and flag & 8192) {
                                       while (idx > -1) {
                                         if (std::abs(pdgs[idx]) == 6 and flags[idx] & 8192)
                                           return (pdg > 0) ? 3 : 4;

                                         // mother pdg is not top, but also not the same as self
                                         // probably the result of some other decay instead of radiation
                                         if (pdgs[idx] != pdg)
                                           return 0;

                                         // otherwise we update the mother index to mother of mother and recheck
                                         idx = idxs[idx];
                                       }
                                     }

                                     // bottom quark block
                                     // parton level so simply check that immediate mother is a final copy top
                                     if (std::abs(pdg) == 5) {
                                       if (std::abs(pdgs[idx]) == 6 and flags[idx] & 8192)
                                         return (pdg > 0) ? 5 : 6;
                                     }

                                     // leptonic W daughter block
                                     // parton level so the immediate mother needs to be a W that comes from top
                                     if (std::abs(pdg) > 10 and std::abs(pdg) < 17) {
                                       if (std::abs(pdgs[idx]) == 24 and flags[idx] & 8192) {
                                         idx = idxs[idx];

                                         while (idx > -1) {
                                           if (std::abs(pdgs[idx]) == 6 and flags[idx] & 8192) {
                                             const bool thirdgen = std::abs(pdg) > 14;
                                             const bool charged = std::abs(pdg) % 2 == 1;
                                             const bool particle = pdg > 0;

                                             if (!thirdgen and charged and particle)
                                               return 8;
                                             if (!thirdgen and charged and !particle)
                                               return 7;

                                             if (!thirdgen and !charged and particle)
                                               return 9;
                                             if (!thirdgen and !charged and !particle)
                                               return 10;

                                             if (thirdgen and charged and particle)
                                               return 12;
                                             if (thirdgen and charged and !particle)
                                               return 11;

                                             if (thirdgen and !charged and particle)
                                               return 13;
                                             if (thirdgen and !charged and !particle)
                                               return 14;
                                           }

                                           idx = idxs[idx];
                                         }
                                       }
                                     }

                                     // the b muons
                                     // FIXME what about b -> gb -> b'b' b(...) -> u'b'(...) b(...)?
                                     // FIXME multi-matches b -> uc -> uus -> (...)?
                                     if (std::abs(pdg) == 13) {
                                       while (idx > -1) {
                                         if (std::abs(pdgs[idx]) == 5) {
                                           if (std::abs(pdgs[idxs[idx]]) == 6 and flags[idxs[idx]] & 8192)
                                             return (pdgs[idx] > 0) ? 15 : 16;
                                         }

                                         idx = idxs[idx];
                                       }
                                     }

                                     // everything else not relevant for us
                                     return 0;
                                   }, "pdg", "flag", "mother");

  // to get around the issue of naod saving only some mass and not the other
  gen_particle.transform_attribute("mass", [] (int tt2l, float mass, int pdg) -> float {
      if (tt2l < 5)
        return mass;

      if (tt2l == 5 or tt2l == 6)
        return constants::m_bottom<>;

      if (tt2l == 7 or tt2l == 8) {
        if (std::abs(pdg) == 11)
          return constants::m_electron<>;
        if (std::abs(pdg) == 13)
          return constants::m_muon<>;
      }

      if (tt2l == 11 or tt2l == 12)
        return constants::m_tau<>;

      if (tt2l == 15 or tt2l == 16)
        return constants::m_muon<>;

      return mass;
    }, "dileptonic_ttbar", "gpmass", "pdg");

  Collection<unsigned char, boolean, int, float> electron("electron", "nElectron", 16, 16);
  electron.add_attribute("pt", "Electron_pt", 1.f);
  electron.add_attribute("eta", "Electron_eta", 1.f);
  electron.add_attribute("phi", "Electron_phi", 1.f);
  electron.add_attribute("mass", "Electron_mass", 1.f);
  electron.add_attribute("charge", "Electron_charge", 1);
  electron.add_attribute("deta_sc", "Electron_deltaEtaSC", 1.f);
  electron.transform_attribute("sc_eta", std::plus<float>(), "deta_sc", "eta");
  electron.add_attribute("id_cutbased", "Electron_cutBased", 1);

  Collection<unsigned char, boolean, int, float> muon("muon", "nMuon", 16, 16);
  muon.add_attribute("pt", "Muon_pt", 1.f);
  muon.add_attribute("eta", "Muon_eta", 1.f);
  muon.add_attribute("phi", "Muon_phi", 1.f);
  muon.add_attribute("mass", "Muon_mass", 1.f);
  muon.add_attribute("charge", "Muon_charge", 1);
  muon.add_attribute("sip3d", "Muon_sip3d", 1.f);
  muon.add_attribute("ip3d", "Muon_ip3d", 1.f);
  muon.add_attribute("dxy", "Muon_dxy", 1.f);
  muon.add_attribute("dxyErr", "Muon_dxyErr", 1.f);
  muon.transform_attribute("sip2d", [] (float dxy, float dxyErr) { return std::abs(dxy / dxyErr); }, "dxy", "dxyErr");
  muon.add_attribute("dxybs", "Muon_dxybs", 1.f);
  muon.add_attribute("id_cutloose", "Muon_looseId", true);
  muon.add_attribute("gpf", "Muon_genPartFlav", std::numeric_limits<unsigned char>::max());

  dat.associate(meta, gen_particle, electron, muon);

  Aggregate euu("euu", 18, 1, electron, muon, muon);
  euu.set_indexer([] (auto &ge, auto &gu, auto &gb) -> std::vector<std::array<int, 3>> {
      if (ge.n_elements() < 1 or gu.n_elements() < 2)
        return {};

      auto idxe = ge.sort_descending("pt");
      auto idxu = gu.indices();
      idxu.clear();

      auto ipdb = gb.sort_absolute_descending("sip3d");
      for (int iu = 1; iu < ipdb.size(); ++iu)
        idxu.emplace_back(ipdb[iu]);

      gu.update_indices(idxu);
      idxu = gu.sort_descending("pt");
      idxu.emplace_back(ipdb[0]);
      gu.update_indices(idxu);

      return {{idxe[0], idxu[0], ipdb[0]}};
    });
  euu.add_attribute("e_pt", identity<>, "electron::pt");
  euu.add_attribute("e_eta", identity<>, "electron::eta");
  euu.add_attribute("e_phi", identity<>, "electron::phi");

  euu.add_attribute("u_pt", identity<>, "muon::pt");
  euu.add_attribute("u_eta", identity<>, "muon::eta");
  euu.add_attribute("u_phi", identity<>, "muon::phi");

  euu.add_attribute("b_pt", [] (float, float f2) {return f2;}, "muon::pt", "muon::pt");
  euu.add_attribute("b_eta", [] (float, float f2) {return f2;}, "muon::eta", "muon::eta");
  euu.add_attribute("b_phi", [] (float, float f2) {return f2;}, "muon::phi", "muon::phi");

  euu.add_attribute("eu_deta", absolute_difference<>, "electron::eta", "muon::eta");
  euu.add_attribute("eu_dphi", dphi<>, "electron::phi", "muon::phi");
  euu.transform_attribute("eu_dR", quadratic_sum<>, "eu_deta", "eu_dphi");
  euu.add_attribute("eu_charge", std::plus<int>(), "electron::charge", "muon::charge");
  euu.add_attribute("eu_mass", invariant_mass<2>(),
                    "electron::pt", "electron::eta", "electron::phi", "electron::mass", 
                    "muon::pt", "muon::eta", "muon::phi", "muon::mass");
  euu.add_attribute("lb_deta", [] (float eeta, int eq, float ueta, int uq, float beta, int bq) {
      if (eq * bq < 0)
        return absolute_difference(eeta, beta);
      else if (uq * bq < 0)
        return absolute_difference(ueta, beta);

      return -9999.f;
    }, "electron::eta", "electron::charge", "muon::eta", "muon::charge", "muon::eta", "muon::charge");
  euu.add_attribute("lb_dphi", [] (float ephi, int eq, float uphi, int uq, float bphi, int bq) {
      if (eq * bq < 0)
        return dphi(ephi, bphi);
      else if (uq * bq < 0)
        return dphi(uphi, bphi);

      return -9999.f;
    }, "electron::phi", "electron::charge", "muon::phi", "muon::charge", "muon::phi", "muon::charge");
  euu.add_attribute("lb_mass", [] (float ept, float eeta, float ephi, float em, int eq,
                                   float upt, float ueta, float uphi, float um, int uq,
                                   float bpt, float beta, float bphi, float bm, int bq) -> float {
                      static TLorentzVector p1, p2;
                      p1.SetPtEtaPhiM(bpt, beta, bphi, bm);

                      if (eq * bq < 0)
                        p2.SetPtEtaPhiM(ept, eeta, ephi, em);
                      else if (uq * bq < 0)
                        p2.SetPtEtaPhiM(upt, ueta, uphi, um);
                      else
                        return -9999.f;

                      return (p1 + p2).M();
                    }, 
                    "electron::pt", "electron::eta", "electron::phi", "electron::mass", "electron::charge", 
                    "muon::pt", "muon::eta", "muon::phi", "muon::mass", "muon::charge",
                    "muon::pt", "muon::eta", "muon::phi", "muon::mass", "muon::charge");

  euu.add_attribute("eub_mass", invariant_mass<3>(),
                    "electron::pt", "electron::eta", "electron::phi", "electron::mass", 
                    "muon::pt", "muon::eta", "muon::phi", "muon::mass",
                    "muon::pt", "muon::eta", "muon::phi", "muon::mass");

  Histogram hist_euu;
  hist_euu.make_histogram<TH1I>(filler_count(electron), "ele_nele", "", 10, 0, 10);
  hist_euu.make_histogram<TH1I>(filler_count(muon), "muo_nmuo", "", 10, 0, 10);
  hist_euu.make_histogram<TH2I>(filler_count(electron, muon), "euu_nele_nmuo", "", 10, 0, 10, 10, 0, 10);

  hist_euu.make_histogram<TH1F>(filler_first_of(electron, "pt"), "ele_pt_1", "", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_all_of(electron, "pt"), "ele_pt_a", "", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(electron, "eta"), "ele_eta_1", "", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH1F>(filler_all_of(electron, "eta"), "ele_eta_a", "", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH2F>(filler_first_of(electron, "pt", "eta"), "ele_pt_eta_1", "", 60, 0.f, 300.f, 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH2F>(filler_all_of(electron, "pt", "eta"), "ele_pt_eta_a", "", 60, 0.f, 300.f, 50, -2.5f, 2.5f);

  hist_euu.make_histogram<TH1F>(filler_first_of(muon, "pt"), "muo_pt_1", "", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_all_of(muon, "pt"), "muo_pt_a", "", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(muon, "eta"), "muo_eta_1", "", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH1F>(filler_all_of(muon, "eta"), "muo_eta_a", "", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH2F>(filler_first_of(muon, "pt", "eta"), "muo_pt_eta_1", "", 60, 0.f, 300.f, 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH2F>(filler_all_of(muon, "pt", "eta"), "muo_pt_eta_a", "", 60, 0.f, 300.f, 50, -2.5f, 2.5f);

  static constexpr float pif = constants::pi<>;
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "e_pt"), "euu_e_pt_euu", ";e pt", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "e_eta"), "euu_e_eta_euu", ";e eta", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "e_phi"), "euu_e_phi_euu", ";e phi", 48, -pif, pif);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "u_pt"), "euu_u_pt_euu", ";u pt", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "u_eta"), "euu_u_eta_euu", ";u eta", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "u_phi"), "euu_u_phi_euu", ";u phi", 48, -pif, pif);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "b_pt"), "euu_b_pt_euu", ";b pt", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "b_eta"), "euu_b_eta_euu", ";b eta", 50, -2.5f, 2.5f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "b_phi"), "euu_b_phi_euu", ";b phi", 48, -pif, pif);

  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "eu_dphi"), "euu_eu_dphi", ";eu dphi", 24, 0.f, pif);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "eu_deta"), "euu_eu_deta", ";eu deta", 24, 0.f, 6.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "eu_mass"), "euu_eu_mass", ";eu mass", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH2F>(filler_first_of(euu, "eu_deta", "eu_dphi"), "euu_eu_deta_dphi", ";eu deta;eu dphi", 24, 0.f, 6.f, 24, 0.f, pif);

  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "lb_dphi"), "euu_lb_dphi", ";lb dphi", 24, 0.f, pif);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "lb_deta"), "euu_lb_deta", ";lb deta", 24, 0.f, 6.f);
  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "lb_mass"), "euu_lb_mass", ";lb mass", 60, 0.f, 300.f);
  hist_euu.make_histogram<TH2F>(filler_first_of(euu, "lb_deta", "lb_dphi"), "euu_lb_deta_dphi", ";lb deta;lb dphi", 24, 0.f, 6.f, 24, 0.f, pif);

  hist_euu.make_histogram<TH1F>(filler_first_of(euu, "eub_mass"), "euu_eub_mass", ";eub mass", 100, 0.f, 1000.f);

  Aggregate gen_llu("gen_llu", 22, 1, gen_particle, gen_particle, gen_particle);
  gen_llu.set_indexer([] (auto &g1, auto &, auto &)
                      -> std::vector<std::array<int, 3>> {
                        auto lepton = g1.filter_equal("dileptonic_ttbar", 8);
                        auto antilepton = g1.filter_equal("dileptonic_ttbar", 7);
                        auto bottom = g1.filter_equal("dileptonic_ttbar", 15);
                        auto antibottom = g1.filter_equal("dileptonic_ttbar", 16);

                        if (lepton.size() != 1 or antilepton.size() != 1)
                          return {};
                        if (bottom.size() == 0 and antibottom.size() == 0)
                          return {};
                        // this agg is for single-sided matches only
                        if (bottom.size() != 0 and antibottom.size() != 0)
                          return {};

                        // grab the b muon with the highest pt
                        // for now this proxies the muon "closest in decay chain" to the b
                        auto idxs = g1.indices();
                        g1.update_indices((bottom.size() != 0) ? bottom : antibottom);
                        bottom = g1.sort_descending("pt");
                        g1.update_indices(idxs);

                        // can be either bottom or antibottom due to the above
                        return {{lepton[0], antilepton[0], bottom[0]}};
                      });
  gen_llu.add_attribute("plep_pdg", identity<int>, "gen_particle::pdg");
  gen_llu.add_attribute("alep_pdg", [] (int, int pdg) { return pdg; }, "gen_particle::pdg", "gen_particle::pdg");

  gen_llu.add_attribute("lb_correct_charge", [] (int , int , int id, int , int , int pdg) {
      int charge = (pdg > 0) ? -1 : 1;
      return (id == 15) ? charge + 1 : charge - 1;
    }, 
    "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
    "gen_particle::pdg", "gen_particle::pdg", "gen_particle::pdg");
  gen_llu.add_attribute("lb_correct_mass", [] (int , int , int id,
                                               float pt1, float eta1, float phi1, float m1,
                                               float pt2, float eta2, float phi2, float m2,
                                               float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &p4 = (id == 15) ? 
                            multivector_system(pt3, eta3, phi3, m3, pt2, eta2, phi2, m2) :
                            multivector_system(pt3, eta3, phi3, m3, pt1, eta1, phi1, m1);

                          return p4.M();
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("lb_correct_pt", [] (int , int , int id,
                                             float pt1, float eta1, float phi1, float m1,
                                             float pt2, float eta2, float phi2, float m2,
                                             float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &p4 = (id == 15) ? 
                            multivector_system(pt3, eta3, phi3, m3, pt2, eta2, phi2, m2) :
                            multivector_system(pt3, eta3, phi3, m3, pt1, eta1, phi1, m1);

                          return p4.Pt();
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("lb_correct_rapidity", [] (int , int , int id,
                                                   float pt1, float eta1, float phi1, float m1,
                                                   float pt2, float eta2, float phi2, float m2,
                                                   float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &p4 = (id == 15) ? 
                            multivector_system(pt3, eta3, phi3, m3, pt2, eta2, phi2, m2) :
                            multivector_system(pt3, eta3, phi3, m3, pt1, eta1, phi1, m1);

                          return p4.Rapidity();
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("lb_correct_deta", [] (int , int , int id,
                                               float eta1, float eta2, float eta3) {
                          return (id == 15) ? absolute_difference(eta3, eta2) : absolute_difference(eta3, eta1);
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");
  gen_llu.add_attribute("lb_correct_dphi", [] (int , int , int id,
                                               float phi1, float phi2, float phi3) {
                          return (id == 15) ? dphi(phi3, phi2) : dphi(phi3, phi1);
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::phi", "gen_particle::phi", "gen_particle::phi");
  gen_llu.transform_attribute("lb_correct_dR", quadratic_sum<>, "lb_correct_deta", "lb_correct_dphi");
  gen_llu.add_attribute("lb_wrong_charge", [] (int , int , int id, int , int , int pdg) {
      int charge = (pdg > 0) ? -1 : 1;
      return (id == 16) ? charge + 1 : charge - 1;
    }, 
    "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
    "gen_particle::pdg", "gen_particle::pdg", "gen_particle::pdg");
  gen_llu.add_attribute("lb_wrong_mass", [] (int , int , int id,
                                             float pt1, float eta1, float phi1, float m1,
                                             float pt2, float eta2, float phi2, float m2,
                                             float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &p4 = (id == 16) ? 
                            multivector_system(pt3, eta3, phi3, m3, pt2, eta2, phi2, m2) :
                            multivector_system(pt3, eta3, phi3, m3, pt1, eta1, phi1, m1);

                          return p4.M();
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("lb_wrong_pt", [] (int , int , int id,
                                           float pt1, float eta1, float phi1, float m1,
                                           float pt2, float eta2, float phi2, float m2,
                                           float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &p4 = (id == 16) ? 
                            multivector_system(pt3, eta3, phi3, m3, pt2, eta2, phi2, m2) :
                            multivector_system(pt3, eta3, phi3, m3, pt1, eta1, phi1, m1);

                          return p4.Pt();
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("lb_wrong_rapidity", [] (int , int , int id,
                                                 float pt1, float eta1, float phi1, float m1,
                                                 float pt2, float eta2, float phi2, float m2,
                                                 float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &p4 = (id == 16) ? 
                            multivector_system(pt3, eta3, phi3, m3, pt2, eta2, phi2, m2) :
                            multivector_system(pt3, eta3, phi3, m3, pt1, eta1, phi1, m1);

                          return p4.Rapidity();
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("lb_wrong_deta", [] (int , int , int id,
                                             float eta1, float eta2, float eta3) {
                          return (id == 16) ? absolute_difference(eta3, eta2) : absolute_difference(eta3, eta1);
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");
  gen_llu.add_attribute("lb_wrong_dphi", [] (int , int , int id,
                                             float phi1, float phi2, float phi3) {
                          return (id == 16) ? dphi(phi3, phi2) : dphi(phi3, phi1);
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::phi", "gen_particle::phi", "gen_particle::phi");
  gen_llu.transform_attribute("lb_wrong_dR", quadratic_sum<>, "lb_wrong_deta", "lb_wrong_dphi");
  gen_llu.add_attribute("ll_mass", invariant_mass<2>(), 
                         "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                         "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("ll_pt", system_pt<2>(), 
                         "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                         "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("ll_rapidity", system_rapidity<2>(), 
                         "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                         "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("ll_deta", absolute_difference<>, 
                         "gen_particle::eta", "gen_particle::eta");
  gen_llu.add_attribute("ll_dphi", dphi<>, 
                         "gen_particle::phi", "gen_particle::phi");
  gen_llu.add_attribute("llb_mass", invariant_mass<3>(), 
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("llb_pt", system_pt<3>(), 
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  gen_llu.add_attribute("llb_rapidity", system_rapidity<3>(), 
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  // does this have a good discrimination?
  auto f_lb_scatter = [] (TLorentzVector &pl, const TLorentzVector &plb) -> float {
    pl.Boost( -1. * plb.BoostVector() );
    return pl.Vect().Unit().Dot( plb.Vect().Unit() );
  };

  gen_llu.add_attribute("lb_correct_scatter", [&f_lb_scatter] (int , int , int id,
                                                  float pt1, float eta1, float phi1, float m1,
                                                  float pt2, float eta2, float phi2, float m2,
                                                  float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &pb = multivector_system(pt3, eta3, phi3, m3);
                          static TLorentzVector pl, plb;
                          if (id == 15)
                            pl.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                          else
                            pl.SetPtEtaPhiM(pt1, eta1, phi1, m1);

                          plb = pl + pb;
                          return f_lb_scatter(pl, plb);
                        }, 
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_llu.add_attribute("lb_wrong_scatter", [&f_lb_scatter] (int , int , int id,
                                                             float pt1, float eta1, float phi1, float m1,
                                                             float pt2, float eta2, float phi2, float m2,
                                                             float pt3, float eta3, float phi3, float m3) -> float {
                          const auto &pb = multivector_system(pt3, eta3, phi3, m3);
                          static TLorentzVector pl, plb;
                          if (id == 16)
                            pl.SetPtEtaPhiM(pt2, eta2, phi2, m2);
                          else
                            pl.SetPtEtaPhiM(pt1, eta1, phi1, m1);

                          plb = pl + pb;
                          return f_lb_scatter(pl, plb);
                        },
                        "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar", "gen_particle::dileptonic_ttbar",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");


  Histogram hist_llu;
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_charge"), "gen_llu_lb_correct_charge", ";charge", 5, -2.5f, 2.5f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_mass"), "gen_llu_lb_correct_mass", ";mass", 60, 0.f, 300.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_pt"), "gen_llu_lb_correct_pt", ";pt", 60, 0.f, 300.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_rapidity"), "gen_llu_lb_correct_rapidity", ";rapidity", 80, -4.f, 4.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_deta"), "gen_llu_lb_correct_deta", ";deta", 24, 0.f, 6.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_dphi"), "gen_llu_lb_correct_dphi", ";dphi", 24, 0.f, pif);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_correct_dR"), "gen_llu_lb_correct_dR", ";dR", 24, 0.f, 6.8f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_deta", "lb_correct_dphi"), 
                                "gen_llu_lb_correct_deta_dphi", ";deta;dphi", 24, 0.f, 6.f, 24, 0.f, pif);

  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_charge"), "gen_llu_lb_wrong_charge", ";charge", 5, -2.5f, 2.5f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_mass"), "gen_llu_lb_wrong_mass", ";mass", 60, 0.f, 300.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_pt"), "gen_llu_lb_wrong_pt", ";pt", 60, 0.f, 300.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_rapidity"), "gen_llu_lb_wrong_rapidity", ";rapidity", 80, -4.f, 4.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_deta"), "gen_llu_lb_wrong_deta", ";deta", 24, 0.f, 6.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_dphi"), "gen_llu_lb_wrong_dphi", ";dphi", 24, 0.f, pif);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "lb_wrong_dR"), "gen_llu_lb_wrong_dR", ";dR", 24, 0.f, 6.8f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_wrong_deta", "lb_wrong_dphi"), 
                                "gen_llu_lb_wrong_deta_dphi", ";deta;dphi", 24, 0.f, 6.f, 24, 0.f, pif);

  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_charge", "lb_wrong_charge"), 
                                "gen_llu_lb_charge_correct_wrong", ";correct charge;wrong charge", 5, -2.5f, 2.5f, 5, -2.5f, 2.5f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_mass", "lb_wrong_mass"), 
                                "gen_llu_lb_mass_correct_wrong", ";correct mass;wrong mass", 60, 0.f, 300.f, 60, 0.f, 300.f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_pt", "lb_wrong_pt"), 
                                "gen_llu_lb_pt_correct_wrong", ";correct pt;wrong pt", 60, 0.f, 300.f, 60, 0.f, 300.f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_rapidity", "lb_wrong_rapidity"), 
                                "gen_llu_lb_rapidity_correct_wrong", ";correct rapidity;wrong rapidity", 80, -4.f, 4.f, 80, -4.f, 4.f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_deta", "lb_wrong_deta"), 
                                "gen_llu_lb_deta_correct_wrong", ";correct deta;wrong deta", 24, 0.f, 6.f, 24, 0.f, 6.f);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_dphi", "lb_wrong_dphi"), 
                                "gen_llu_lb_dphi_correct_wrong", ";correct dphi;wrong dphi", 24, 0.f, pif, 24, 0.f, pif);
  hist_llu.make_histogram<TH2F>(filler_first_of(gen_llu, "lb_correct_dR", "lb_wrong_dR"), 
                                "gen_llu_lb_dR_correct_wrong", ";correct dR;wrong dR", 24, 0.f, 6.8f, 24, 0.f, 6.8f);

  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "ll_mass"), "gen_llu_ll_mass", "", 120, 0.f, 600.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "ll_pt"), "gen_llu_ll_pt", "", 80, 0.f, 400.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "ll_rapidity"), "gen_llu_ll_rapidity", "", 80, -4.f, 4.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "ll_deta"), "gen_llu_ll_deta", "", 24, 0.f, 6.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "ll_dphi"), "gen_llu_ll_dphi", "", 24, 0.f, pif);

  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "llb_mass"), "gen_llu_llb_mass", "", 160, 0.f, 800.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "llb_pt"), "gen_llu_llb_pt", "", 60, 0.f, 300.f);
  hist_llu.make_histogram<TH1F>(filler_first_of(gen_llu, "llb_rapidity"), "gen_llu_llb_rapidity", "", 80, -4.f, 4.f);

  Tree tree_llu("tree_smtt_bpark_llu.root", "tree");
  tree_llu.make_single_branches(gen_llu, 
                                "lb_correct_charge", "lb_correct_mass", "lb_correct_pt", "lb_correct_rapidity",
                                "lb_correct_deta", "lb_correct_dphi", "lb_correct_dR",
                                "lb_wrong_charge", "lb_wrong_mass", "lb_wrong_pt", "lb_wrong_rapidity",
                                "lb_wrong_deta", "lb_wrong_dphi", "lb_wrong_dR",
                                "lb_correct_scatter", "lb_wrong_scatter",
                                "ll_mass", "ll_pt", "ll_rapidity",
                                "ll_deta", "ll_dphi",
                                "llb_mass", "llb_pt", "llb_rapidity");
  //tree_llu.make_array_branches(gen_particle, "mass", "pt", "eta", "phi", "pdg", "dileptonic_ttbar");
  tree_llu.make_array_branches(muon, "pt", "gpf");

  /*/ FIXME finish this
  Aggregate gr_llu("gr_llu", 22, 1, gen_llu, euu, gen_llu, euu, gen_llu, euu);
  gr_llu.set_indexer([] (auto &gen, auto &reco, auto &, auto &) -> std::vector<std::array<int, 6>> {
      if (gen.n_elements() != 1 or reco.n_elements() != 0)
        return {};

      auto &ig = gen.ref_to_indices();
      auto &ir = reco.ref_to_indices();

      if (gen.count_equal("plep_pdg", 11) == 1 and gen.count_equal("alep_pdg", -13) == 1)
        return {{ig[0], ir[0], ig[1], ir[1], ig[2], ir[2]}};

      if (gen.count_equal("plep_pdg", 13) == 1 and gen.count_equal("alep_pdg", -11) == 1)
        return {{ig[1], ir[0], ig[0], ir[1], ig[2], ir[2]}};
    });
  gr_llu.add_attribute("e_deta", absolute_difference<>, "gen_llu::e_eta");
  */
  //auto printer = [] (auto &p) {std::cout << int(p) << "\n";};

  std::array<int, 2> pass_trigger_3l = {0, 0};
  auto f_analyze = [&pass_trigger_3l, &hlt, &meta, &gen_particle, &electron, &muon, &euu, &hist_euu, &gen_llu, &hist_llu, &tree_llu] 
    (long long entry) {
    meta.populate(entry);

    gen_particle.populate(entry);
    gen_llu.populate(entry);
    if (gen_llu.n_elements()) {
      hist_llu.fill();
      tree_llu.fill();
    }

    bool pass_trigger = false;
    for (auto &path : hlt) {
      auto bit = meta.filter_equal(path, true);

      if (bit.size())
        pass_trigger = true;
    }

    if (pass_trigger)
      pass_trigger_3l[0] += 1;
    else
      return;

    muon.populate(entry);
    muon.update_indices( muon.filter_greater("pt", 5.f) );
    muon.update_indices( muon.filter_in("eta", -2.4f, 2.4f) );
    muon.update_indices( muon.filter_equal("id_cutloose", true) );
    if (muon.n_elements() < 2)
      return;

    electron.populate(entry);
    electron.update_indices( electron.filter_greater("pt", 5.f) );
    electron.update_indices( electron.filter_in("sc_eta", -2.4f, 2.4f) );
    electron.update_indices( electron.filter_greater("id_cutbased", 0) );
    electron.update_indices( electron.sort_descending("pt") );
    if (electron.n_elements() < 1)
      return;

    pass_trigger_3l[1] += 1;

    euu.populate(entry);
    if (euu.n_elements())
      hist_euu.fill();

    //muon.iterate(printer, -1, -1, "gpf");
  };

  dat.set_analyzer(f_analyze);
  dat.analyze();
  std::cout << pass_trigger_3l[0] << " pass trigger\n";
  std::cout << pass_trigger_3l[1] << " also contain euu triplet" << std::endl;

  tree_llu.save();
  save_all_as("hist_smtt_bpark.root", hist_llu, hist_euu);

  return 0;
}

