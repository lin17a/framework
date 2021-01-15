// -*- C++ -*-
// author: afiq anuar
// short: functions pertaining to locating and making lists of input datasets for the analysis
// note: currently the assumption is that all datasets are located locally at DESY

#ifndef FWK_INPUT_DATASET_H
#define FWK_INPUT_DATASET_H

#include "TSystem.h"
#include "TString.h"

// poor man's version of file finder by extension
std::vector<std::string> file_by_ext(const std::string &dir, const std::string &ext)
{
  // which really relies on ROOT's ability to run shell commands aha
  TString allfile = gSystem->GetFromPipe(("find " + dir + " -type f -name '*" + ext + "'").c_str());
  TString file;
  Ssiz_t index = 0;

  std::vector<std::string> v_file;
  while (allfile.Tokenize(file, index, "\n"))
    v_file.emplace_back(file.Data());

  return v_file;
}



std::vector<std::string> file_list(const std::string &dataset)
{
  // a list of known datasets and their location
  static const std::string prefix = "/pnfs/desy.de/cms/tier2/store/";

  static const std::vector<std::array<std::string, 2>> v_dataset = {
    {"DoubleMuon-2018A", 
     "data/Run2018A/DoubleMuon/NANOAOD/02Apr2020-v1/"},

    {"DoubleMuon-2018B", 
     "data/Run2018B/DoubleMuon/NANOAOD/02Apr2020-v1/"},

    {"DoubleMuon-2018C", 
     "data/Run2018C/DoubleMuon/NANOAOD/02Apr2020-v1/"},

    {"DoubleMuon-2018D", 
     "data/Run2018D/DoubleMuon/NANOAOD/02Apr2020-v1/"},

    {"EGamma-2018A", 
     "data/Run2018A/EGamma/NANOAOD/02Apr2020-v1/"},

    {"EGamma-2018B", 
     "data/Run2018B/EGamma/NANOAOD/02Apr2020-v1/"},

    {"EGamma-2018C", 
     "data/Run2018C/EGamma/NANOAOD/02Apr2020-v1/"},

    {"EGamma-2018D", 
     "data/Run2018D/EGamma/NANOAOD/02Apr2020-v1/"},

    {"MuonEG-2018A", 
     "/data/Run2018A/MuonEG/NANOAOD/02Apr2020-v1/"},

    {"MuonEG-2018B", 
     "/data/Run2018B/MuonEG/NANOAOD/02Apr2020-v1/"},

    {"MuonEG-2018C", 
     "/data/Run2018C/MuonEG/NANOAOD/02Apr2020-v1/"},

    {"MuonEG-2018D", 
     "/data/Run2018D/MuonEG/NANOAOD/02Apr2020-v1/"},

    {"SingleMuon-2018A", 
     "/data/Run2018A/SingleMuon/NANOAOD/02Apr2020-v1/"},

    {"SingleMuon-2018B", 
     "/data/Run2018B/SingleMuon/NANOAOD/02Apr2020-v1/"},

    {"SingleMuon-2018C", 
     "/data/Run2018C/SingleMuon/NANOAOD/02Apr2020-v1/"},

    {"SingleMuon-2018D", 
     "/data/Run2018D/SingleMuon/NANOAOD/02Apr2020-v1/"},

    {"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},


    {"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext2-v1/"},

    {"ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/"},

    {"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/"},

    {"TTTo2L2Nu_TuneCP5CR1_QCDbased_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5CR1_QCDbased_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_TuneCP5CR2_GluonMove_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5CR2_GluonMove_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_TuneCP5_erdON_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_erdON_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_TuneCP5down_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5down_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_TuneCP5up_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5up_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_hdampDOWN_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_mtop169p5_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_mtop169p5_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTTo2L2Nu_mtop175p5_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_mtop175p5_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_TuneCP5CR1_QCDbased_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_TuneCP5CR1_QCDbased_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_TuneCP5CR2_GluonMove_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_TuneCP5CR2_GluonMove_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/"},

    {"TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/"},

    {"TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"WW_TuneCP5_13TeV-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"WZ_TuneCP5_13TeV-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/WZ_TuneCP5_13TeV-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"},

    {"ZZ_TuneCP5_13TeV-pythia8", 
     "mc/RunIIAutumn18NanoAODv7/ZZ_TuneCP5_13TeV-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/"}
  };

  for (const auto &d : v_dataset) {
    if (dataset == d[0])
      return file_by_ext(prefix + d[1], ".root");
  }

  return {};
}

#endif
