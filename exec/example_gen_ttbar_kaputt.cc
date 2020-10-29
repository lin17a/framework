// example execution macro showing a generator-level ttbar analysis assuming the CMS nanoAOD format
// compile:
// filename=example_gen_ttbar; g++ $(root-config --cflags --evelibs) -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -I ../plugins/ -I ../src/ -o ${filename} ${filename}.cc

// core framework headers
#include "Dataset.h"
#include "Collection.h"
#include "Aggregate.h"
#include "Histogram.h"
#include "../src/Tree.h"

// additional headers that aid in defining analysis-dependent functions
#include "TLorentzVector.h"
#include "misc/function_util.h"
#include "misc/numeric_vector.h"
#include "misc/spin_correlation.h"

// things I included
#include <dirent.h>
#include <string>


using attr_func_type_ttbar = const std::function<float(float,float,float,float,float,float,float,float)>;

attr_func_type_ttbar get_attr_function(const std::function<float(TLorentzVector)> extract_attr_func) {
  return [extract_attr_func](float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2) {
    TLorentzVector p1, p2;
    p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    TLorentzVector sum = p1 + p2;
    return extract_attr_func(sum);
  };
  }

// declare in advance a few functions we will need in the analysis
float invariant_mass(float pt1, float eta1, float phi1, float mass1,
                     float pt2, float eta2, float phi2, float mass2)
{
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);

  return (p1 + p2).M();
}

float pt_vector_sum(float pt1, float phi1, float pt2, float phi2)
{
  float px1 = pt1 * std::cos(phi1), px2 = pt2 * std::cos(phi2);
  float py1 = pt1 * std::sin(phi1), py2 = pt2 * std::sin(phi2);

  return std::sqrt(((px1 + px2) * (px1 + px2)) + ((py1 + py2) * (py1 + py2)));
}


float angle(float pt1, float eta1, float phi1, float mass1,
            float pt2, float eta2, float phi2, float mass2)
{ 
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  TVector3 p1_T3 = p1.Vect();
  TVector3 p2_T3 = p2.Vect();
  return p1_T3.Angle(p2_T3);
}

int main() {
  // the core part of the framework are all within this namespace
  using namespace Framework;

  // first and foremost, we specify the input files we will be looking at
  // this is done by constructing a dataset object
  // in this example we will analyze flat trees, so we specify that our dataset is of type TChain
  // the first argument to the constructor is the dataset name, and second is the name of the TTree we will be analyzing
  // note: Dataset<TChain> is used also when analyzing one ROOT file, since Dataset<TTree> is dedicated to text file analysis
  Dataset<TChain> dat("mc", "Events");
  // add the files to be analyzed
  // dat.add_file("/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/022107FA-F567-1B44-B139-A18ADC996FCF.root");
  // dat.add_file("/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/0EF179F9-428D-B944-8DB3-63E04ED9AE8E.root");
  // dat.add_file("/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/0EF179F9-428D-B944-8DB3-63E04ED9AE8E.root");

  std::string dir_string("/nfs/dust/cms/user/meyerlin/ba/runs/ttbarlo_normal/root_files/");
  struct dirent *entry = nullptr;
  DIR *dp = opendir(dir_string.c_str());
  if (dp == nullptr) {
     perror("opendir: Path does not exist or could not be read.");
     return EXIT_FAILURE;
  }
  while ((entry = readdir(dp))) {
       puts(entry->d_name);
       std::string file_string = dir_string + std::string(entry->d_name);
       dat.add_file(file_string);
  }
  closedir(dp);

  // next step is to specify the attributes to be included in the analysis
  // for this we will make use of two data structures, collections and aggregates
  // both inherit from the group data structure, which at this stage can be taken as a set of attributes that belong together
  // a notion that will be made clearer as we go through this example 

  // a collection is used when the attributes are directly read off flat ROOT files
  // here we initialize a non-array type collection
  // which is used when the read branches are single numbers per event
  // C++ is a strongly-typed language, so we need to specify the types of attributes we want the collection to cover within the <> bracket
  // the arguments of the constructor are the name of the collection
  // and an integer specifying the number of attributes the collection is expected to contain
  // which serves as a hint for the memory layout to be allocated for the analysis
  // the framework is able to adapt this layout as the need arises, but it is good practice to give an accurate estimate
  Collection<uint, unsigned long long, float> metadata("metadata", 5);

  // here we add an attribute to the collection
  // the arguments are the name of the attribute, the name of the branch associated to this attribute, and a hint number
  // the exact value of the hint number does not matter, but it needs to be of the same type as the data contained in the branch
  // the type of the branch can be obtained by calling Print() on the tree
  // in this case the type is an unsigned int i.e. uint 
  metadata.add_attribute("run", "run", 1U);

  // we add second and third attributes
  // attribute names (the first argument) must be unique within each Collection
  metadata.add_attribute("lumi", "luminosityBlock", 1U);
  metadata.add_attribute("event", "event", 1ULL);

  // it's worth remarking here that the default type for literal floating point number is double in C++
  // to obtain a literal float one needs to append an 'f' suffix at the end
  metadata.add_attribute("weight", "genWeight", 1.f);
  metadata.add_attribute("lhe_orixwgtup", "LHEWeight_originalXWGTUP", 1.f);

  // next we initialize an array-type collection
  // the constructor arguments in this case are: 
  // 1- the collection name
  // 2- the name of a non-array branch that holds the number of elements each array contains in each branch
  // 3- the number of attributes the collection is expected to contain
  // 4- the number of elements each array is expected to contain 
  // as with 3- the framework adapts the memory layout as the need arises, but it is good practice to provide a helpful hint
  // note the type list within <>, for technical reasons we can not use the type bool for boolean branches
  // but use the custom boolean type instead, which functions the same for us 
  Collection<boolean, int, float> gen_particle("gen_particle", "nGenPart", 11, 256);
  gen_particle.add_attribute("mass", "GenPart_mass", 1.f);
  gen_particle.add_attribute("pt", "GenPart_pt", 1.f);
  gen_particle.add_attribute("eta", "GenPart_eta", 1.f);
  gen_particle.add_attribute("phi", "GenPart_phi", 1.f);
  gen_particle.add_attribute("pdg", "GenPart_pdgId", 1);
  gen_particle.add_attribute("status", "GenPart_status", 1);
  gen_particle.add_attribute("flag", "GenPart_statusFlags", 1);
  gen_particle.add_attribute("mother", "GenPart_genPartIdxMother", 1);

  // on top of adding attributes directly read from the branches
  // we can also add attributes which are transformed from existing attributes
  // for example, in generator-level analysis we are typically interested in final copy particles of a particular pdg id
  // in this case we provide a boolean flag for whether the particle is a final copy top (anti-)quark
  // the arguments for transformed attributes are:
  // 1- the name of the new attribute
  // 2- a function, which can be either a lambda or a regular function, to transform the input arguments into a final value
  // 3- the list of attributes that serve as input arguments to 2-
  // note that transformed attribute has to be well-definable for every element in the collection
  gen_particle.transform_attribute("final_top", 
                                   // explicitly specify a boolean return type
                                   // otherwise true and false default to bool
                                   [] (int pdg, int flag) -> boolean {
                                     // 8192 is 2^13, which means isLastCopy in nanoAOD flag bitset
                                     if (std::abs(pdg) == 6 and flag & 8192)
                                       return true;
                                     else return false;
                                   }, "pdg", "flag");

  // often we are interested not just in the current particle attributes, but also in its provenance
  // for example if we want to tag only final state W bosons coming from top quark decays 
  // in nanoAOD the branch GenPart_genPartIdxMother contains the index of the mother of the current particle
  // the usage of which requires us to view the attribute arrays in their entirety
  // however the collection deals with its elements one by one, and similarly the functions take as arguments only the current attributes
  // so it may not be obvious how this can be achieved given these restrictions
  // this example shows how, exploiting the lambda captures to create an impure function
  gen_particle.transform_attribute("final_w_top_daughter", 
                                   // we capture references to the full array of mother indices and the final top tag
                                   // ordering matters; we can capture final_top only after defining it
                                   [&idxs = gen_particle.get<int>("mother"), &tops = gen_particle.get<boolean>("final_top")] 
                                   (int pdg, int flag, int idx) -> boolean {
                                     // first the finality check similar to the top case
                                     if (std::abs(pdg) == 24 and flag & 8192) {
                                       // the particle history log may contain radiation
                                       // which is written as the mother particle having the same pdg id as the particle itself
                                       // as such it is not sufficient to check the immediate mother of the particle
                                       // mother index == -1 is nanoAOD for mother is not saved in the array 
                                       while (idx > -1) {
                                         // we are done if the mother is a top
                                         if (tops[idx])
                                           return true;

                                         // otherwise we update the mother index to mother of mother and recheck
                                         idx = idxs[idx];
                                       }
                                     }
                                     return false;
                                   }, "pdg", "flag", "mother");

  // of course, the transformation can be as complex as desired
  // in this case we tag the entire set of particles of interest in a dileptonic ttbar decay tt -> WbWb -> lvblvb
  // restricting leptons to only electron or muon as is commonly done in experimental analyses
  gen_particle.transform_attribute("dileptonic_ttbar", 
                                   [&pdgs = gen_particle.get<int>("pdg"), &idxs = gen_particle.get<int>("mother"), 
                                    &flags = gen_particle.get<int>("flag")] 
                                   (int pdg, int flag, int idx) -> int {
                                     // integer flag for particles which are part of a generator-level dileptonic ttbar system
                                     // 1 top
                                     // 2 W+
                                     // 3 bottom (parton)
                                     // 4 antilepton - e, mu only
                                     // 5 neutrino - e, mu only 
                                     // 6 antitop
                                     // 7 W-
                                     // 8 antibottom (parton)
                                     // 9 lepton - e, mu only
                                     // 10 antineutrino - e, mu only

                                     // final top quarks
                                     if (pdg == 6 and flag & 8192)
                                       return 1;
                                     if (pdg == -6 and flag & 8192)
                                       return 6;

                                     // W boson block
                                     if (std::abs(pdg) == 24 and flag & 8192) {
                                       while (idx > -1) {
                                         if (std::abs(pdgs[idx]) == 6 and flags[idx] & 8192) {
                                           if (pdg > 0)
                                             return 2;
                                           else
                                             return 7;
                                         }

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
                                       if (std::abs(pdgs[idx]) == 6 and flags[idx] & 8192) {
                                         if (pdg > 0)
                                           return 3;
                                         else
                                           return 8;
                                       }
                                     }

                                     // leptonic W daughter block
                                     if (std::abs(pdg) > 10 and std::abs(pdg) < 15) {
                                       if (std::abs(pdgs[idx]) == 24 and flags[idx] & 8192) {
                                         idx = idxs[idx];

                                         while (idx > -1) {
                                           if (std::abs(pdgs[idx]) == 6 and flags[idx] & 8192) {
                                             if (pdg % 2) {
                                               if (pdg > 0)
                                                 return 9;
                                               else 
                                                 return 4;
                                             }
                                             else {
                                               if (pdg > 0)
                                                 return 5;
                                               else 
                                                 return 10;
                                             }
                                           }

                                           idx = idxs[idx];
                                         }
                                       }
                                     }

                                     // everything else not relevant for us
                                     return 0;
                                   }, "pdg", "flag", "mother");

  // having specified all the branches we are interested in, we associate the collections with the dataset
  // this is done by the call below, where the arguments are simply all the collections we are considering
  // this call is equivalent to SetBranchAddress(...) etc steps in a more traditional flat tree analyses
  // be sure to include all the collections in the call, as step-wise association is currently not supported
  dat.associate(metadata, gen_particle);

  // now we move to the case of attributes that are well-defined only for some selection of elements from the collections
  // for example, the invariant mass of the system of final top quark pair is relevant only for gen_particle with attribute dileptonic_ttbar == 1 or 6
  // for this we have aggregate, which is used to combine multiple not necessarily distinct groups according to some arbitrary indexing rule
  // the arguments for the constructor are:
  // 1- the name of the aggregate
  // 2- the number of attributes it is expected to contain
  // 3- the number of elements each array is expected to contain
  // 4- the collections that contribute an index to the aggregate (specified once per index)
  // here we want the 2-particle system made of the final top quark pair, so we provide the gen_particle twice, once for top and once for antitop
  // it is worth noting that currently an aggregate can only be made from groups of the same type
  Aggregate gen_ttbar("gen_ttbar", 7, 1, gen_particle, gen_particle);

  // when using aggregates, one must specify the indexing rule i.e. how the elements from the underlying groups are to be combined
  // this is done by providing a function, whose arguments are references to the groups
  // the return type of the function is a vector of array of indices; the array size corresponds to the number of underlying index
  // the first argument is identified with the first group given to the aggregate constructor and so on
  gen_ttbar.set_indexer([] (const auto &g1, const auto &g2)
                       -> std::vector<std::array<int, 2>> {
                          // we use the tags defined above to find the top and antitop
                          // filter_XXX returns a vector of indices of elements fullfilling the criteria
                          // list of currently supported filter operations are in the group header file
                          auto top = g1.filter_equal("dileptonic_ttbar", 1);
                          auto antitop = g2.filter_equal("dileptonic_ttbar", 6);

                          // check that the collection contains exactly one top and one antitop
                          // if not the case, return an empty index list
                          if (top.size() != 1 or antitop.size() != 1)
                            return {};

                          return {{top[0], antitop[0]}};
                       });

  // having specified the index, we can now specify the attributes
  // for aggregates, the argument to add_attributes are:
  // 1- the name of the attribute
  // 2- the function to calculate the attribute from the underlying group attributes
  // we now provide a regular function pointer to 2-, unlike lambda function as we did before
  // both are supported in either case
  // 3- the list of underlying group attributes in the same order as used in 2-
  // the syntax is underlying_group::attribute
  // when the same underlying attribute is used multiple times, the aggregate takes care of the indexing
  // here the first gen_particle::pt is read off the top index, and the second time from the antitop index
  gen_ttbar.add_attribute("ttbar_mass", invariant_mass<2, float>, 
                          "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                          "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_ttbar.add_attribute("ttbar_pt", pt_vector_sum, 
                          "gen_particle::pt", "gen_particle::phi", "gen_particle::pt", "gen_particle::phi");
  
  // to make ttbar add_attributes shorter
  // const auto add_attribute = [&](const std::string &funcname, attr_func_type_ttbar func) {
  //  gen_ttbar.add_attribute(funcname, func, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
  // 					    "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");  
  // };

  
//  add_attribute("ttbar_mass", get_attr_function([](TLorentzVector sum){return sum.M();}));
//   add_attribute("ttbar_rapidity", get_attr_function([](TLorentzVector sum){return sum.Rapidity();}));
//   add_attribute("ttbar_beta", get_attr_function([](TLorentzVector sum){return sum.Beta();}));
//   add_attribute("ttbar_pseudo_rapidity", get_attr_function([](TLorentzVector sum){return sum.PseudoRapidity();}));
//   add_attribute("ttbar_gamma", get_attr_function([](TLorentzVector sum){return sum.Gamma();}));
//   add_attribute("ttbar_e", get_attr_function([](TLorentzVector sum){return sum.E();}));
//   add_attribute("ttbar_energy", get_attr_function([](TLorentzVector sum){return sum.Energy();}));
//   add_attribute("ttbar_et", get_attr_function([](TLorentzVector sum){return sum.Et();}));
//   add_attribute("ttbar_et2", get_attr_function([](TLorentzVector sum){return sum.Et2();}));
//   add_attribute("ttbar_m2", get_attr_function([](TLorentzVector sum){return sum.M2();}));
//   add_attribute("ttbar_mag", get_attr_function([](TLorentzVector sum){return sum.Mag();}));
//   add_attribute("ttbar_mag2", get_attr_function([](TLorentzVector sum){return sum.Mag2();}));
//   add_attribute("ttbar_minus", get_attr_function([](TLorentzVector sum){return sum.Minus();}));
//   add_attribute("ttbar_mt", get_attr_function([](TLorentzVector sum){return sum.Mt();}));
//   add_attribute("ttbar_mt2", get_attr_function([](TLorentzVector sum){return sum.Mt2();}));
//   add_attribute("ttbar_p", get_attr_function([](TLorentzVector sum){return sum.P();}));
//   add_attribute("ttbar_perp", get_attr_function([](TLorentzVector sum){return sum.Perp();}));
//   add_attribute("ttbar_perp2", get_attr_function([](TLorentzVector sum){return sum.Perp2();}));
//   add_attribute("ttbar_phi", get_attr_function([](TLorentzVector sum){return sum.Phi();}));
//   add_attribute("ttbar_plus", get_attr_function([](TLorentzVector sum){return sum.Plus();}));
//   add_attribute("ttbar_px", get_attr_function([](TLorentzVector sum){return sum.Px();}));
//   add_attribute("ttbar_py", get_attr_function([](TLorentzVector sum){return sum.Py();}));
//   add_attribute("ttbar_pz", get_attr_function([](TLorentzVector sum){return sum.Pz();}));
//   add_attribute("ttbar_rho", get_attr_function([](TLorentzVector sum){return sum.Rho();}));
//   add_attribute("ttbar_t", get_attr_function([](TLorentzVector sum){return sum.T();}));
//   add_attribute("ttbar_theta", get_attr_function([](TLorentzVector sum){return sum.Theta();}));
//   add_attribute("ttbar_x", get_attr_function([](TLorentzVector sum){return sum.X();}));
//   add_attribute("ttbar_y", get_attr_function([](TLorentzVector sum){return sum.Y();}));
//   add_attribute("ttbar_z", get_attr_function([](TLorentzVector sum){return sum.Z();}));
  // possibly redundant, but simple copies of the underlying attribute at the index of interest is also possible
  // by masking the index through unnamed arguments
  // we will do this twice in this example, so let's declare the lambdas in advance
  auto return_first = [] (float f1, float ) {return f1;};
  auto return_second = [] (float , float f2) {return f2;};

  // of course for the first index no masking is necessary, we can simply call the relevant attributes once
  gen_ttbar.add_attribute("top_pt", return_first, "gen_particle::pt", "gen_particle::pt");
  // but for second index onwards masking is needed
  gen_ttbar.add_attribute("antitop_pt", return_second, "gen_particle::pt", "gen_particle::pt");

  // do the same thing for phi
  gen_ttbar.add_attribute("top_phi", return_first, "gen_particle::phi", "gen_particle::phi");
  gen_ttbar.add_attribute("antitop_phi", return_second, "gen_particle::phi", "gen_particle::phi");

  // and redundantly obtain the ttbar pt again
  gen_ttbar.transform_attribute("ttbar_pt_transform", pt_vector_sum, 
                                "top_pt", "top_phi", "antitop_pt", "antitop_phi");

  // now clearly the index masking examples above are rather artificial, serving only to highlight features
  // rather than something one might actually be interested in doing in an actual analysis
  // so now let us consider a more typical example of a dileptonic ttbar system
  // where we are interested in six particles: tops, charged leptons and bottoms
  Aggregate gen_tt_ll_bb("gen_tt_ll_bb", 15, 1, gen_particle, gen_particle, gen_particle, gen_particle, gen_particle, gen_particle);

  // set the indices similarly as above
  gen_tt_ll_bb.set_indexer([] (const auto &g1, const auto &g2, const auto &g3, const auto &g4, const auto &g5, const auto &g6)
                       -> std::vector<std::array<int, 6>> {
                          auto top = g1.filter_equal("dileptonic_ttbar", 1);
                          auto antitop = g2.filter_equal("dileptonic_ttbar", 6);
                          auto lepton = g3.filter_equal("dileptonic_ttbar", 9);
                          auto antilepton = g4.filter_equal("dileptonic_ttbar", 4);
                          auto bottom = g5.filter_equal("dileptonic_ttbar", 3);
                          auto antibottom = g6.filter_equal("dileptonic_ttbar", 8);

                          // check that the collection contains exactly one particle of interest
                          // if not, return an empty index list
                          if (top.size() != 1 or antitop.size() != 1 or lepton.size() != 1 or antilepton.size() != 1 or 
                              bottom.size() != 1 or antibottom.size() != 1)
                            return {};

                          return {{top[0], antitop[0], lepton[0], antilepton[0], bottom[0], antibottom[0]}};
                       });

  // for transferring the attributes of some daughter particles
  auto return_third = [] (float , float , float f3) {return f3;};
  auto return_fourth = [] (float , float , float , float f4) {return f4;};
  auto return_fifth = [] (float , float , float , float , float f5) {return f5;};
  auto return_sixth = [] (float , float , float , float , float , float f6) {return f6;};

  gen_tt_ll_bb.add_attribute("lepton_pt", return_third, 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("lepton_eta", return_third, 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  gen_tt_ll_bb.add_attribute("antilepton_pt", return_fourth, 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("antilepton_eta", return_fourth, 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  gen_tt_ll_bb.add_attribute("bottom_pt", return_fifth, 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("bottom_eta", return_fifth, 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  gen_tt_ll_bb.add_attribute("antibottom_pt", return_sixth, 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("antibottom_eta", return_sixth, 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  // let's also consider some invariant masses
  gen_tt_ll_bb.add_attribute("ttbar_mass", invariant_mass<2, float>, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("llbar_mass", [] (float , float , float , float , 
                                               float , float , float , float ,
                                               float pt1, float eta1, float phi1, float m1,
                                               float pt2, float eta2, float phi2, float m2)
                             { return invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("bbbar_mass", [] (float , float , float , float , 
                                               float , float , float , float ,
                                               float , float , float , float , 
                                               float , float , float , float ,
                                               float pt1, float eta1, float phi1, float m1,
                                               float pt2, float eta2, float phi2, float m2)
                             { return invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("lbbar_mass", [] (float , float , float , float , 
                                               float , float , float , float ,
                                               float pt1, float eta1, float phi1, float m1,
                                               float , float , float , float , 
                                               float , float , float , float ,
                                               float pt2, float eta2, float phi2, float m2)
                             { return invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("lbarb_mass", [] (float , float , float , float , 
                                               float , float , float , float ,
                                               float , float , float , float ,
                                               float pt1, float eta1, float phi1, float m1,
                                               float pt2, float eta2, float phi2, float m2, 
                                               float , float , float , float)
                             { return invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("lb_mass", [] (float , float , float , float , 
                                            float , float , float , float ,
                                            float pt1, float eta1, float phi1, float m1,
                                            float , float , float , float ,
                                            float pt2, float eta2, float phi2, float m2, 
                                            float , float , float , float)
                             { return invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("lbarbbar_mass", [] (float , float , float , float , 
                                                  float , float , float , float ,
                                                  float , float , float , float ,
                                                  float pt1, float eta1, float phi1, float m1,
                                                  float , float , float , float ,
                                                  float pt2, float eta2, float phi2, float m2)
                             { return invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  // let's histogram the attributes we defined above
  // this is done through the histogram class, which handles a group of histograms sharing the same weights and to be filled at the same time
  // in this example we will have two instances, before and after some acceptance cuts
  Histogram hist_no_cut;

  // we define an argument-less impure function that computes the weight for each histogram entry
  // here we only take the per-event weight from the metadata collection
  // we can see here that internally a non-array collection is in fact an array collection of size 1
  // if no weighter is defined, histograms are filled with weight 1
  hist_no_cut.set_weighter([&weight = metadata.get<float>("weight")] () { return weight[0]; });

  // next we define the histograms, where the histogram type are given inside the <> bracket
  // all histogram types supported by ROOT are supported
  // the first argument is an impure function instructing how the histograms should be filled
  // the function takes two arguments, a histogram pointer and a weight
  // the remaining arguments are those expected by ROOT histogram constructor of the type being used
  hist_no_cut.make_histogram<TH1F>([&gen_ttbar] (TH1F *hist, double weight) {
      // do not fill if the aggregate has no elements
      if (gen_ttbar.n_elements() != 1)
        return;

      auto &mass = gen_ttbar.get<float>("ttbar_mass");
      hist->Fill(mass[0], weight);
    }, "ttbar_mass_no_cut", "", 120, 300.f, 1500.f);

  // in many cases we will be filling the histograms in similar ways
  // e.g. check for presence, and if yes, fill the first/all elements
  // and having to write out the filling function every time can be cumbersome
  // so in the plugins some utility functions are provided for these commonly used functions
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt"), "ttbar_pt_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt_transform"), "ttbar_pt_transform_no_cut", "", 120, 0.f, 1200.f);

  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_pt"), "lepton_pt_no_cut", "", 100, 0.f, 400.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_eta"), "lepton_eta_no_cut", "", 100, -5.f, 5.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_pt"), "antilepton_pt_no_cut", "", 100, 0.f, 400.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_eta"), "antilepton_eta_no_cut", "", 100, -5.f, 5.f);

  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_pt"), "bottom_pt_no_cut", "", 100, 0.f, 400.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_eta"), "bottom_eta_no_cut", "", 100, -5.f, 5.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_pt"), "antibottom_pt_no_cut", "", 100, 0.f, 400.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_eta"), "antibottom_eta_no_cut", "", 100, -5.f, 5.f);

  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ttbar_mass"), "ttbar_mass_2_no_cut", "", 120, 300.f, 1500.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_mass"), "llbar_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bbbar_mass"), "bbbar_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lb_mass"), "lb_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarbbar_mass"), "lbarbbar_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbbar_mass"), "lbbar_mass_no_cut", "", 100, 0.f, 200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarb_mass"), "lbarb_mass_no_cut", "", 100, 0.f, 200.f);

  // let's define another histogram instance but now with acceptance cuts
  // we can, but don't need to, define the cuts in the filling function themselves
  // so the histogram instance is defined identically as above except the histogram names
  Histogram hist_cut;
  hist_cut.set_weighter([&weight = metadata.get<float>("weight")] () { return weight[0]; });
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_mass"), "ttbar_mass_cut", "", 120, 300.f, 1500.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt"), "ttbar_pt_cut", "", 120, 0.f, 1200.f);

  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_pt"), "lepton_pt_cut", "", 100, 0.f, 400.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_eta"), "lepton_eta_cut", "", 100, -5.f, 5.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_pt"), "antilepton_pt_cut", "", 100, 0.f, 400.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_eta"), "antilepton_eta_cut", "", 100, -5.f, 5.f);

  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_pt"), "bottom_pt_cut", "", 100, 0.f, 400.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_eta"), "bottom_eta_cut", "", 100, -5.f, 5.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_pt"), "antibottom_pt_cut", "", 100, 0.f, 400.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_eta"), "antibottom_eta_cut", "", 100, -5.f, 5.f);

  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ttbar_mass"), "ttbar_mass_2_cut", "", 120, 300.f, 1500.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_mass"), "llbar_mass_cut", "", 120, 0.f, 1200.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bbbar_mass"), "bbbar_mass_cut", "", 120, 0.f, 1200.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lb_mass"), "lb_mass_cut", "", 120, 0.f, 1200.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarbbar_mass"), "lbarbbar_mass_cut", "", 120, 0.f, 1200.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbbar_mass"), "lbbar_mass_cut", "", 100, 0.f, 200.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarb_mass"), "lbarb_mass_cut", "", 100, 0.f, 200.f);

  // if the unbinned values are needed we can save them as flat trees
  // just like the histogram object we start by instantiating the object
  // args are the file and tree names we want to save out
  // optionally also the compression setting
  Tree tree_gen("gen_smtt_kinematic_cut.root", "tree");

  // two types of branches are supported - single and array
  // by calling the respective methods as shown below
  // the args are the group contributing and its attributes that one would like to be saved
  // gen_ttbar aggregate always have only one element, so it's suitably saved as single branches
  // the branches are saved with name groupname_attribute
  tree_gen.make_single_branches(gen_ttbar, "ttbar_mass", "ttbar_pt");

  // use array branches when the group can have more than one element and we are interested in them all
  // in this case, in addition to the attribute array branches, one also gets the n_groupname branch
  tree_gen.make_array_branches(gen_particle, "mass", "pt", "eta", "phi", "pdg", "dileptonic_ttbar");

  // to add more branches simply call the make_*_branches again
  // note: each group can contribute branches to a tree exactly once
  tree_gen.make_single_branches(gen_tt_ll_bb, "llbar_mass", "bbbar_mass", "lb_mass", "lbarbbar_mass", "lbbar_mass", "lbarb_mass");

  // so far we have defined the inputs we would like to read, how to transform them and the form of our analysis output 
  // there is one last piece of preparation we need to do
  // and that is to define how we would like our analysis to be performed
  // by this point you can probably guess that this is achieved by defining an impure function
  // that captures the references to all the collections, aggregates and histograms we defined above
  // the only argument to this function is the entry number
  // one way to think about this function is that it contains the instructions on how to analyze a single event
  auto f_analyze = [&metadata, &gen_particle, &gen_ttbar, &gen_tt_ll_bb, &hist_no_cut, &hist_cut, &tree_gen] (long long entry) {
    // first we start by populating the collections
    // this is essentially equivalent of the tree->GetEntry(entry)
    // with the (compulsory) freedom of timing the call separately for each group
    metadata.populate(entry);
    gen_particle.populate(entry);

    // since the collections serve as input to the aggregates, they need to be populated first
    gen_ttbar.populate(entry);
    gen_tt_ll_bb.populate(entry);

    // we make an oversimplification here, considering only the events where gen_tt_ll_bb contain an element
    // this is because in the above, we have grouped the gen_ttbar and gen_tt_ll_bb histograms together
    // despite the fact that the requirements of gen_tt_ll_bb is strictly tighter than gen_ttbar
    // without this restriction, when filling the hist_no_cut histograms we get spurious entries in the gen_tt_ll_bb histograms
    // when this aggregate is empty e.g. when we have taus in the event
    if (!gen_tt_ll_bb.n_elements())
      return;

    // fill the no (acceptance) cut histograms
    hist_no_cut.fill();

    // as promised above we would like some acceptance cuts 
    // which we impose on the charged leptons and bottom quarks
    // we could have used the daughters that we transferred to the gen_tt_ll_bb aggregate
    // but to highlight some additional features we will instead do it through the gen_particle collections instead
    // begin by selecting the daughters among all the gen particles using a generic filter method
    // which needs a function that evaluates to true or false based on a list of attributes
    auto lepton_bottom_passing_pt_eta_cut = gen_particle.filter([] (int tag, float pt, float eta) {
        // not lepton or bottom, reject
        if (tag != 9 and tag != 4 and tag != 3 and tag != 8)
          return false;

        if (pt > 20.f and std::abs(eta) < 2.4f)
          return true;
        else
          return false;
      }, "dileptonic_ttbar", "pt", "eta");

    // recall that filter methods return a list of indices
    // to overwrite the indices list of the group, we use the update_indices method
    gen_particle.update_indices(lepton_bottom_passing_pt_eta_cut);

    // if all four objects pass the cut, then gen_particle will have 4 elements left
    // fill also our tree at this point
    if (gen_particle.n_elements() == 4) {
      hist_cut.fill();
      tree_gen.fill();
    }

    /*/ here is the way to perform equivalent filtering using the gen_tt_ll_bb aggregate
    // by stacking multiple filter_XXX calls
    // do check that it gives equivalent results as the above!
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_greater("lepton_pt", 20.f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_in("lepton_eta", -2.4f, 2.4f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_greater("antilepton_pt", 20.f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_in("antilepton_eta", -2.4f, 2.4f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_greater("bottom_pt", 20.f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_in("bottom_eta", -2.4f, 2.4f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_greater("antibottom_pt", 20.f) );
    gen_tt_ll_bb.update_indices( gen_tt_ll_bb.filter_in("antibottom_eta", -2.4f, 2.4f) );

    if (gen_tt_ll_bb.n_elements() == 1)
      hist_cut.fill();
    */
  };

  // tell the dataset instance about our event analyzer function
  dat.set_analyzer(f_analyze);

  // and run it!
  // for analyzing only a subset, provide as argument the desired number of events
  dat.analyze();

  // when all is said and done, we collect the output
  // which we can plot, or perform statistical tests etc
  hist_no_cut.save_as("hist_no_cut.root");
  hist_cut.save_as("hist_cut.root");
  tree_gen.save();

  return 0;
}

