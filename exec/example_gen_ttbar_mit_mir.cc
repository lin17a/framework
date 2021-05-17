// example execution macro showing a generator-level ttbar analysis assuming the CMS nanoAOD format
// compile from exec/ directory:
// make example_gen_ttbar -f ../Makefile
// or from anywhere else:
// make example_gen_ttbar FWK_BASE_DIR="<FWK_DIR>" -f <FWK_DIR>/Makefile

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
#include "misc/reweighting_new.h"

// things I included
#include <dirent.h>
#include <string>
#include <sstream>
#include <iomanip>
// this one is to add light weight paticle masses, if it would work
// #include "./misc/constants.h"


// enum to string things



std::string calc_weight_version_to_str(calc_weight_version version){
using cwv = calc_weight_version;
  switch ( version )
      {
         case cwv::juan_paper:
            return "juan_paper";
         case cwv::juan_code:
            return "juan_code";
         case cwv::juan_paper_different_M2:
            return "juan_paper_different_M2";
         default:
            return "undefined_version";
      }
}

std::string higgs_type_to_str(higgs_type_t type){
using ht = higgs_type_t;
  switch ( type )
      {
         case ht::scalar:
            return "scalar";
         case ht::pseudo_scalar:
            return "pseudo_scalar";
         default:
            return "undefined_higgs_type";
      }
}

std::string res_int_to_str(res_int_t res_int){
//using ri = res_int;
  switch ( res_int )
      {
         case res_int_t::resonance:
            return "resonance";
         case res_int_t::interference:
            return "interference";
         case res_int_t::both:
            return "res_int";
         default:
            return "undefined_res_int_type";
      }
}

std::string create_filename(std::string custom_prefix, higgs_type_t higgs_type, float mass, float width, calc_weight_version calc_weight_version, res_int_t res_int, std::string cut){
 
    //std::string folder = 

    std::string calc_weight_version_str = calc_weight_version_to_str(calc_weight_version);
    std::string higgs_type_str = higgs_type_to_str(higgs_type);
    std::string res_int_str = res_int_to_str(res_int);
  
    std::stringstream filename;
    filename << custom_prefix << "_" << higgs_type_str << "_m" << std::fixed << std::setprecision(0) << mass << "_w" << width << "_" << calc_weight_version_str << "_" << res_int_str << "_" << cut << ".root";
    std::cout << "filename: " << filename.str() << std::endl;
    return filename.str();
}

// 
// cool things from Ruben to make the add_attribute shorter
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

// declare in advance a few functions we will need in the analysis
float system_invariant_mass(float pt1, float eta1, float phi1, float mass1,
                            float pt2, float eta2, float phi2, float mass2)
{
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);

  return (p1 + p2).M();
}

// declare in advance a few functions we will need in the analysis
float system_invariant_mass_four_particles(float pt1, float eta1, float phi1, float mass1,
		                           float pt2, float eta2, float phi2, float mass2,
		                           float pt3, float eta3, float phi3, float mass3,
		                           float pt4, float eta4, float phi4, float mass4)
{
  TLorentzVector p1, p2, p3, p4;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  p3.SetPtEtaPhiM(pt3, eta3, phi3, mass3);
  p4.SetPtEtaPhiM(pt4, eta4, phi4, mass4);

  return (p1 + p2 + p3 + p4).M();
}

// calculate invariant mass of just one particle 
float invariant_mass_one_particle(float pt1, float eta1, float phi1, float mass1)
{
    TLorentzVector p1;
    p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    // std::cout << p1.M() << "\t"  << pt1 << "\t" << eta1 << "\t" << phi1 << "\t" << mass1 <<  "\n";    return p1.M();
    return p1.M();
}

float pt_vector_sum(float pt1, float phi1, float pt2, float phi2)
{
  float px1 = pt1 * std::cos(phi1), px2 = pt2 * std::cos(phi2);
  float py1 = pt1 * std::sin(phi1), py2 = pt2 * std::sin(phi2);

  return std::sqrt(((px1 + px2) * (px1 + px2)) + ((py1 + py2) * (py1 + py2)));
}

float cos_theta_old(float pt1, float eta1, float phi1, float mass1,
                float pt2, float eta2, float phi2, float mass2)
{
    TLorentzVector p1, p2, sum_p1_p2;
    p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    sum_p1_p2 = p1 + p2; 
    TVector3 zBoost(0.,0.,-sum_p1_p2.Pz()/sum_p1_p2.E()); sum_p1_p2.Boost(zBoost);
    TVector3 transverseBoost(-sum_p1_p2.Px()/sum_p1_p2.E(),-sum_p1_p2.Py()/sum_p1_p2.E(),0.); sum_p1_p2.Boost(transverseBoost);
    p1.Boost(zBoost); 
    p1.Boost(transverseBoost);
    return p1.CosTheta();
}

//TLorentzVector f_zmf_tt (TLorentzVector p4, TLorentzVector p4lab_TT) {
//    p4.Boost( -1. * p4lab_TT.BoostVector() );
//    return p4;
//}

float cpTTT(float pt1, float eta1, float phi1, float mass1,
                float pt2, float eta2, float phi2, float mass2)
{
    TLorentzVector p4lab_pTop, p4lab_aTop, p4lab_TT, p4hel_pTop; 
    p4lab_pTop.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    p4lab_aTop.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    p4lab_TT = p4lab_pTop + p4lab_aTop; 
    p4hel_pTop = f_zmf_tt(p4lab_pTop, p4lab_TT);
    float cpTTT = p4hel_pTop.Vect().Unit().Dot( p4lab_TT.Vect().Unit() ); 
    return cpTTT;
}


float diff_cos_theta_cpTP(float pt1, float eta1, float phi1, float mass1,
                          float pt2, float eta2, float phi2, float mass2)
{
    TLorentzVector p1, p2, sum_p1_p2;
    p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
    p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
    sum_p1_p2 = p1 + p2; 
    static const TVector3 zBase(0., 0., 1.);
    const TLorentzVector p4hel_pTop = f_zmf_tt(p1, sum_p1_p2);
    float cpTP = p4hel_pTop.Vect().Unit().Dot(zBase);
    
    TVector3 zBoost(0.,0.,-sum_p1_p2.Pz()/sum_p1_p2.E()); sum_p1_p2.Boost(zBoost);
    TVector3 transverseBoost(-sum_p1_p2.Px()/sum_p1_p2.E(),-sum_p1_p2.Py()/sum_p1_p2.E(),0.); sum_p1_p2.Boost(transverseBoost);
    p1.Boost(zBoost); 
    p1.Boost(transverseBoost);
    float cos_theta = p1.CosTheta();

    return cpTP - cos_theta;
}

float my_angle(float pt1, float eta1, float phi1, float mass1,
            float pt2, float eta2, float phi2, float mass2)
{ 
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  TVector3 p1_T3 = p1.Vect();
  TVector3 p2_T3 = p2.Vect();
  return p1_T3.Angle(p2_T3);
}

float llbar_phi(float pt1, float eta1, float phi1, float mass1,
            float pt2, float eta2, float phi2, float mass2,
            float pt3, float eta3, float phi3, float mass3,
            float pt4, float eta4, float phi4, float mass4)
{ 
  TLorentzVector p1, p2, p3, p4;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  p3.SetPtEtaPhiM(pt3, eta3, phi3, mass3);
  p4.SetPtEtaPhiM(pt4, eta4, phi4, mass4);
  return (p3 + p4).Phi();
}

float llbar_dphi(float pt1, float eta1, float phi1, float mass1,
            float pt2, float eta2, float phi2, float mass2,
            float pt3, float eta3, float phi3, float mass3,
            float pt4, float eta4, float phi4, float mass4)
{ 
  TLorentzVector p1, p2, p3, p4;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  p3.SetPtEtaPhiM(pt3, eta3, phi3, mass3);
  p4.SetPtEtaPhiM(pt4, eta4, phi4, mass4);
  return dPhi(p3.Phi(),p4.Phi());
}

// printing out things for all events like b mass or so
auto printer_normal = [] (auto &p) {std::cout << p << "\n";};
auto printer_reweighted = [] (auto &p) {std::cout << "reweighted: " << p << "\n";};

int main(int argc, char **argv) {
//int main() {
  // the core part of the framework are all within this namespace
  // note: an example of CL arguments are provided in bparking/bpark_tt3l.cc file
  // which is not yet integrated into the example
  using namespace Framework;
  int index;
  int c;

  opterr = 0;
  float mass;
  float width;
  higgs_type_t higgs_type = higgs_type_t::undefined;
  res_int_t res_int = res_int_t::undefined;
  std::istringstream ss;
  int calc_resonance = 0;
  int calc_interference = 0;
  int pseudo_scalar = 0;
  int scalar = 0;
  calc_weight_version calc_weight_variant = calc_weight_version::undefined;

  while ((c = getopt (argc, argv, "m:w:ripsv:")) != -1){
    std::string arg = "";
    if (optarg) arg = std::string(optarg);
    switch (c)
      {
      case 'm':{
        std::string mass_str;
        mass_str = arg;
        ss.clear();
        ss.str(mass_str);
        ss >> mass;
        if (ss && ss.eof()){
            std::cout << "mass: " << mass << std::endl;
        } 
        else {
            std::cerr << "mass must be a number." << std::endl;
            return 1;
        }
        break;}
      case 'w':{
        std::string width_str;
        width_str = arg;
        ss.clear();
        ss.str(width_str);
        ss >> width;
        if (ss && ss.eof()){
           std::cout << "width: " << width << std::endl;
        }
        else {
           std::cerr << "width must be a number." << std::endl;
           return 1;}
        break;}
      case 'r':
        calc_resonance = 1;
        res_int = res_int_t::resonance;
        std::cout << "calculating resonance" << std::endl;
        break;
      case 'i':
        calc_interference = 1;
        res_int = res_int_t::interference;
        std::cout << "calculating interference" << std::endl;
        break;
      case 's':
        scalar = 1;
        higgs_type = higgs_type_t::scalar;
        std::cout << "higgs type: scalar" << std::endl; 
        break;
      case 'p':
        pseudo_scalar = 1;
        higgs_type = higgs_type_t::pseudo_scalar;
        std::cout << "higgs type: pseudo scalar" << std::endl; 
        break;
      case 'v':{
        std::string variant_str = arg;
        if (variant_str ==  "juan_code"){
            calc_weight_variant = calc_weight_version::juan_code;
            std::cout << "Using calculation variant juan_code" << std::endl;
        }
        else if (variant_str == "juan_paper") {
            calc_weight_variant = calc_weight_version::juan_paper;
            std::cout << "Using calculation variant juan_paper" << std::endl;
        }
        else {
            std::cerr << "The calculation variant must be one of the following: juan_code or  juan_paper." << std::endl;
            return 1;
        }
        break;}
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
    }
  if (!mass){
    std::cerr << "You have to specify the the mass of the higgs (-m <mass in GeV>)." << std::endl;
  }
  if (!width){
    std::cerr << "You have to specify the the width of the higgs (-w <width>)." << std::endl;
  }
  
  if (calc_resonance + calc_interference == 0){
    std::cerr << "You must choose either resonance (-r) and/or interference (-i)." << std::endl;
    return 0;
  } 
  else if (calc_resonance + calc_interference == 2){
    res_int = res_int_t::both;
  }
  if (scalar + pseudo_scalar != 1){
    std::cerr << "You must choose either scalar (-s) or pseudo scalar (-p)." << std::endl;
    return 0;
  }
  if (calc_weight_variant == calc_weight_version::undefined){
    std::cerr << "You have to specify a calculation method for reweighting." << std::endl;
    return 0;
  } 

  
  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);





  // first and foremost, we specify the input files we will be looking at
  // this is done by constructing a dataset object
  // in this example we will analyze flat trees, so we specify that our dataset is of type TChain
  // the first argument to the constructor is the dataset name, and second is the name of the TTree we will be analyzing
  // note: Dataset<TChain> is used also when analyzing one ROOT file, since Dataset<TTree> is dedicated to text file analysis
  Dataset<TChain> dat("mc", "Events");
  // add the files to be analyzed
  //dat.add_file("/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/022107FA-F567-1B44-B139-A18ADC996FCF.root");
  //dat.add_file("/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/0EF179F9-428D-B944-8DB3-63E04ED9AE8E.root");
  //dat.add_file("/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/0EF179F9-428D-B944-8DB3-63E04ED9AE8E.root");
  //dat.add_file("/nfs/dust/cms/user/meyerlin/ba/runs/ttbarlo_new_nanoaod/root_files/ttbarlo_new_nanoaod__00000.root");
  //dat.add_file("/nfs/dust/cms/user/meyerlin/ba/runs/heavyhiggs_m600_w15_000_RES_PSEUDO_lo_cfg/root_files/heavyhiggs_m600_w15_000_RES_PSEUDO_lo_cfg__00000.root");
  //dat.add_file("/nfs/dust/cms/user/meyerlin/ba/runs/heavyhiggs_m600_w15_000_INT_PSEUDO_lo_cfg/root_files/heavyhiggs_m600_w15_000_INT_PSEUDO_lo_cfg__00000.root");

  std::string dir_string("/nfs/dust/cms/user/meyerlin/ba/runs/ttbarlo_new_nanoaod/root_files/");
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
  
  // 1- the collection name
  // 2- the name of a non-array branch that holds the number of elements each array contains in each branch
  // 3- the number of attributes the collection is expected to contain
  // 4- the number of elements each array is expected to contain 
  // as with 3- the framework adapts the memory layout as the need arises, but it is good practice to provide a helpful hint
  // note the type list within <>, for technical reasons we can not use the type bool for boolean branches
  // but use the custom boolean type instead, which functions the same for us 
  Collection<boolean, int, float, double> gen_particle("gen_particle", "nGenPart", 12, 8192);
  gen_particle.add_attribute("default_mass", "GenPart_mass", 1.f);
  gen_particle.add_attribute("pt", "GenPart_pt", 1.f);
  gen_particle.add_attribute("eta", "GenPart_eta", 1.f);
  gen_particle.add_attribute("phi", "GenPart_phi", 1.f);
  gen_particle.add_attribute("pdg", "GenPart_pdgId", 1);
  gen_particle.add_attribute("status", "GenPart_status", 1);
  gen_particle.add_attribute("flag", "GenPart_statusFlags", 1);
  gen_particle.add_attribute("mother", "GenPart_genPartIdxMother", 1);

  std::cout << "after add attributes to gen-particle" << std::endl << std::endl;;

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
  
  gen_particle.transform_attribute("mass", [] (float mass, int id, int pdg_id) { 
    if (mass != 0) 
      return mass; 
    else if (id == 4 or id == 9){ 
      if (std::abs(pdg_id) == 11)  
    	return 0.00051099895f;   //electron_mass_const
      else if (std::abs(pdg_id) == 13) 
        return 0.1056583745f; } //muon_mass
    return mass;
    }, "default_mass", "pdg", "dileptonic_ttbar" ); 
 
//  std::cout << "after gen_particle transforms" << std::endl;
  
  // add LHE particle Collection for reweighting
  Collection<boolean, int, float, double> lhe_particle("lhe_particle", "nLHEPart", 12, 16);
//  std::cout << "before adding sth" << std::endl;
  lhe_particle.add_attribute("pt", "LHEPart_pt", 1.f);
  lhe_particle.add_attribute("eta", "LHEPart_eta", 1.f);
  lhe_particle.add_attribute("phi", "LHEPart_phi", 1.f);
  lhe_particle.add_attribute("default_mass", "LHEPart_mass", 1.f);
  lhe_particle.add_attribute("pdg", "LHEPart_pdgId", 1);
  lhe_particle.add_attribute("ipz", "LHEPart_incomingpz", 1.f);

  Aggregate lhe_event("lhe_event", 7, 1, gen_particle, gen_particle, lhe_particle, lhe_particle);

//  std::cout << "before lhe_indexer" << std::endl;

  //  when using aggregates, one must specify the indexing rule i.e. how the elements from the underlying groups are to be combined
  //  this is done by providing a function, whose arguments are references to the groups
  //  the return type of the function is a vector of array of indices; the array size corresponds to the number of underlying index
  //  the first argument is identified with the first group given to the aggregate constructor and so on
  lhe_event.set_indexer([&g = lhe_particle] (const auto &g1, const auto &g2, const Group<boolean, int, float, double> &g3, const Group<boolean, int, float, double> &g4)
                       -> std::vector<std::array<int, 4>> {
                          // we use the tags defined above to find the top and antitop
                          // filter_XXX returns a vector of indices of elements fullfilling the criteria
                          // list of currently supported filter operations are in the group header file
                          // g.update_indices( g.filter_not("ipz", 0.f) );
                         
//                          std::cout << "lhe_event indexer beginning" << std::endl;
 
                          auto initial = g.filter_not("ipz", 0.f);
                          auto ngluon = g.filter_equal("pdg", 21, initial).size();
//                          std::cout << "lhe_event indexer after ngluon" << std::endl;
                            
                            
                          auto lepton = merge(g3.filter_equal("pdg", 11), g3.filter_equal("pdg", 13), g3.filter_equal("pdg", 15));
                          auto antilepton = merge(g4.filter_equal("pdg", -11), g4.filter_equal("pdg", -13), g4.filter_equal("pdg", -15));
                          auto bottoms = g.sort_descending("pdg", merge(g.filter_equal("pdg", 5), g.filter_equal("pdg", -5)));
                          auto neutrino = merge(g3.filter_equal("pdg", 12), g3.filter_equal("pdg", 14), g3.filter_equal("pdg", 16));
                          auto antineutrino = merge(g4.filter_equal("pdg", -12), g4.filter_equal("pdg", -14), g4.filter_equal("pdg", -16));
   
//                          float eta_lepton = g3.template get<float>("eta")[lepton[0]]; 
//                          float eta_antilepton = g4.template get<float>("eta")[antilepton[0]];
//                          float phi_lepton = g3.template get<float>("phi")[lepton[0]];
//                          float phi_antilepton = g4.template get<float>("phi")[antilepton[0]];
//                          float pt_lepton = g3.template get<float>("pt")[lepton[0]];
//                          float pt_antilepton = g4.template get<float>("pt")[antilepton[0]];
//                          float m_lepton = g3.template get<float>("default_mass")[lepton[0]];
//                          float m_antilepton = g4.template get<float>("default_mass")[antilepton[0]];
//                          
//                          float eta_bottom = g3.template get<float>("eta")[bottoms[0]]; 
//                          float eta_antibottom = g3.template get<float>("eta")[bottoms[1]];
//                          float phi_bottom = g3.template get<float>("phi")[bottoms[0]];
//                          float phi_antibottom = g3.template get<float>("phi")[bottoms[1]];
//                          float pt_bottom = g3.template get<float>("pt")[bottoms[0]];
//                          float pt_antibottom = g3.template get<float>("pt")[bottoms[1]];
//                          float m_bottom = g3.template get<float>("default_mass")[bottoms[0]];
//                          float m_antibottom = g3.template get<float>("default_mass")[bottoms[1]];
//                          
//                          float eta_neutrino = g3.template get<float>("eta")[neutrino[0]]; 
//                          float eta_antineutrino = g4.template get<float>("eta")[antineutrino[0]];
//                          float phi_neutrino = g3.template get<float>("phi")[neutrino[0]];
//                          float phi_antineutrino = g4.template get<float>("phi")[antineutrino[0]];
//                          float pt_neutrino = g3.template get<float>("pt")[neutrino[0]];
//                          float pt_antineutrino = g4.template get<float>("pt")[antineutrino[0]];
//                          float m_neutrino = g3.template get<float>("default_mass")[neutrino[0]];
//                          float m_antineutrino = g4.template get<float>("default_mass")[antineutrino[0]];
//                            
//                          TLorentzVector lepton_vec, antilepton_vec, bottom_vec, antibottom_vec, neutrino_vec, antineutrino_vec;
//                          lepton_vec.SetPtEtaPhiM(pt_lepton, eta_lepton, phi_lepton, m_lepton);
//                          antilepton_vec.SetPtEtaPhiM(pt_antilepton, eta_antilepton, phi_antilepton, m_antilepton);
//                          bottom_vec.SetPtEtaPhiM(pt_bottom, eta_bottom, phi_bottom, m_bottom);
//                          antibottom_vec.SetPtEtaPhiM(pt_antibottom, eta_antibottom, phi_antibottom, m_antibottom);
//                          neutrino_vec.SetPtEtaPhiM(pt_neutrino, eta_neutrino, phi_neutrino, m_neutrino);
//                          antineutrino_vec.SetPtEtaPhiM(pt_antineutrino, eta_antineutrino, phi_antineutrino, m_antineutrino);
//                         
//                          TLorentzVector sum_lbar_b_n = antilepton_vec + bottom_vec + neutrino_vec;
//                          double_t sum_lbar_b_n_X = sum_lbar_b_n.X();
//                          double_t sum_lbar_b_n_Y = sum_lbar_b_n.Y();
//                          double_t sum_lbar_b_n_Z = sum_lbar_b_n.Z();
//                          double_t sum_lbar_b_n_T = sum_lbar_b_n.T();
//                          std::cout << "sum_lbar_b_n: ( " << sum_lbar_b_n_X << ", " << sum_lbar_b_n_Y << ", " << sum_lbar_b_n_Z << ", " << sum_lbar_b_n_T << " )" << std::endl;
//                          
//                          TLorentzVector sum_l_bbar_nbar = lepton_vec + antibottom_vec + antineutrino_vec;
//                          double_t sum_l_bbar_nbar_X = sum_l_bbar_nbar.X();
//                          double_t sum_l_bbar_nbar_Y = sum_l_bbar_nbar.Y();
//                          double_t sum_l_bbar_nbar_Z = sum_l_bbar_nbar.Z();
//                          double_t sum_l_bbar_nbar_T = sum_l_bbar_nbar.T();
//                          std::cout << "sum_l_bbar_nbar: ( " << sum_l_bbar_nbar_X << ", " << sum_l_bbar_nbar_Y << ", " << sum_l_bbar_nbar_Z << ", " << sum_l_bbar_nbar_T << " )" << std::endl;
  
                          if (ngluon == 1)
                          {
                            std::cout << "lhe_event indexer ngluon == 1" << std::endl;
                            int sign_pz = ((g.get<float>("ipz"))[initial[0]] > 0) ? 1 : -1;
                            
                          
                            //auto tops = g.sort_descending("pdg", merge(g.filter_equal("pdg", 6), g.filter_equal("pdg", -6)));
                            auto top = g1.filter_equal("pdg", 6).filter_equal("status", 22);
                            auto antitop = g2.filter_equal("pdg", -6).filter_equal("status", 22);
                            float eta_top = g1.template get<float>("eta")[top[0]]; 
                            float eta_antitop = g2.template get<float>("eta")[antitop[0]];
                            float phi_top = g1.template get<float>("phi")[top[0]];
                            float phi_antitop = g2.template get<float>("phi")[antitop[0]];
                            float pt_top = g1.template get<float>("pt")[top[0]];
                            float pt_antitop = g2.template get<float>("pt")[antitop[0]];
                            float m_top = g1.template get<float>("default_mass")[top[0]];
                            float m_antitop = g2.template get<float>("default_mass")[antitop[0]];
                            
                            TLorentzVector top_vec, antitop_vec;
                            top_vec.SetPtEtaPhiM(pt_top, eta_top, phi_top, m_top);
                            antitop_vec.SetPtEtaPhiM(pt_antitop, eta_antitop, phi_antitop, m_antitop);
                            
                            
                            TLorentzVector sum = top_vec + antitop_vec;
                            int sign_eta = sum.Eta() > 0 ? 1 : -1;
 
                            if (g.get<int>("pdg")[initial[0]] == 21)
                            {
                              if (sign_eta != sign_pz) {
                                //std::cout << "ngluon = 1 was made to ngluon = 2, first particle gluon" << std::endl;
                                ngluon = 2;}
                              else {
                                //std::cout << "ngluon = 1 was made to ngluon = 0, first particle gluon" << std::endl;
                                ngluon = 0; }
                            } 
                            else
                            {
                              if (sign_eta != sign_pz) {
                                //std::cout << "ngluon = 1 was made to ngluon = 0, second particle gluon" << std::endl;
                                ngluon = 0; }
                              else {
                                //std::cout << "ngluon = 1 was made to ngluon = 2, second particle gluon" << std::endl;
                                ngluon = 2; }
                            }
                          }
                          
                          if (ngluon == 0) { 
                            std::cout << "lhe_event indexer ngluon == 0" << std::endl;
                            return {}; }
                        
                          if (ngluon == 2)
                          {
                            std::cout << "lhe_event indexer ngluon == 2" << std::endl;
                             auto top = g1.filter_equal("pdg", 6).filter_equal("status", 22);
 //                            std::cout << "lhe_event indexer ngluon == 2, after top" << std::endl;
                             auto antitop = g2.filter_equal("pdg", -6).filter_equal("status", 22);
 //                            std::cout << "lhe_event indexer ngluon == 2, after antitop" << std::endl;
                          
                             //auto lepton = g3.filter_equal("pdg", 11, g3.filter_equal("pdg", 13, g3.filter_equal("pdg", 15)));
                             auto lepton = merge(g3.filter_equal("pdg", 11), g3.filter_equal("pdg", 13), g3.filter_equal("pdg", 15));
  //                           std::cout << "lhe_event indexer ngluon == 2, after lepton" << std::endl;
                             // auto lepton = g3.filter_3values("pdg", 11, 13, 15);
                             // auto antilepton = g4.filter_3values("pdg", -11, -13, -15); 
                             auto antilepton = merge(g4.filter_equal("pdg", -11), g4.filter_equal("pdg", -13), g4.filter_equal("pdg", -15));
//                             std::cout << "lhe_event indexer ngluon == 2, after antilepton" << std::endl;
                             //auto antilepton = g4.filter_equal("pdg", -11, g4.filter_equal("pdg", -13, g4.filter_equal("pdg", -15)));
                             //int empty = 0;
//                             std::cout << "top size:" << top.size() << std::endl;
                             if (top.size() != 1 or antitop.size() != 1 or lepton.size() != 1 or antilepton.size() != 1) {
                               return{};}
                             else {std::cout << "es geht!!" << std::endl;}
                             return {{top[0], antitop[0], lepton[0], antilepton[0]}};
                          }
//                          std::cout << "lhe_event indexer end" << std::endl;

                          return {};
                       });
//  std::cout << "after lhe_indexer" << std::endl;
  


  // having specified all the branches we are interested in, we associate the collections with the dataset
  // this is done by the call below, where the arguments are simply all the collections we are considering
  // this call is equivalent to SetBranchAddress(...) etc steps in a more traditional flat tree analyses
  // be sure to include all the collections in the call, as step-wise association is currently not supported
  dat.associate(metadata, gen_particle, lhe_particle);
  
  //calc_weight_version calc_weight_version = calc_weight_version::juan_paper;
  //higgs_type_t higgs_type = higgs_type_t::pseudo_scalar;
  //res_int_t res_int = res_int_t::both;
  //float mass = 400;
  //float width = 20;

  lhe_event.add_attribute("weight_float", calculate_weight<float,float>(calc_weight_variant, higgs_type, mass, width, res_int), 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass",
                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass");
  
  lhe_event.add_attribute("weight_double", calculate_weight<float, double>(calc_weight_variant, higgs_type, mass, width, res_int), 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass",
                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass");
  
//  lhe_event.add_attribute("weight_juan_code_float", calculate_weight<float,float>(calc_weight_version::juan_code, higgs_type), 
//                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
//                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
//                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass",
//                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass");
//  
//  lhe_event.add_attribute("weight_juan_code_double", calculate_weight<float, double>(calc_weight_version::juan_code, higgs_type), 
//                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
//                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
//                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass",
//                             "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass");
//
  // calculate how precise the weights are by comparing double and float precision
  lhe_event.transform_attribute("weight_precision", 
                                   [] (float weight_float, double weight_double) -> double {
                                     double precision = (weight_double - weight_float) / weight_double;
                                     return precision;
                                   }, "weight_float", "weight_double");


  lhe_event.add_attribute("cpTTT_LHE", spin_correlation<>("cpTTT"), "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                                    "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                                    "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass",
                                                                    "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass");
  
  lhe_event.add_attribute("cpTTT_LHE_self_calc", cpTTT, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass");

  lhe_event.add_attribute("t_CosTheta_LHE", cos_theta_old, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                           "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass");

  lhe_event.add_attribute("cpTP_LHE", spin_correlation<>("cpTP"), "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                                  "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                                  "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass",
                                                                  "lhe_particle::pt", "lhe_particle::eta", "lhe_particle::phi", "lhe_particle::default_mass");
  
  lhe_event.add_attribute("diff_cos_theta_cpTP_LHE", diff_cos_theta_cpTP, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass",
                                                                          "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::default_mass");

  // now we move to the case of attributes that are well-defined only for some selection of elements from the collections
  // for example, the invariant mass of the system of final top quark pair is relevant only for gen_particle with attribute dileptonic_ttbar == 1 or 6
  // for this we have aggregate, which is used to combine multiple not necessarily distinct groups according to some arbitrary indexing rule
  // the arguments for the constructor are:
  // 1- the name of the aggregate
  // 2- the number of attributes it is expected to contain
  // 3- the number of elements each array is expected to contain
  // 4- the collections that contribute an index to the aggregate (specified once per index)
  // here we want the 2-particle system made of the finaection{What I do} top quark pair, so we provide the gen_particle twice, once for top and once for antitop
  // it is worth noting that currently an aggregate can only be made from groups of the same type
  Aggregate gen_ttbar("gen_ttbar", 7, 1, gen_particle, gen_particle);

  std::cout << "before gen_ttbar_indexer" << std::endl;

  // when using aggregates, one must specify the indexing rule i.e. how the elements from the underlying groups are to be combined
  // this is done by providing a function, whose arguments are references to the groups
  // the return type of the function is a vector of array of indices; the array size corresponds to the number of underlying index
  // the first argument is identified with the first group given to the aggregate constructor and so on
  gen_ttbar.set_indexer([&g = gen_particle] (const Group<boolean, int, float, double> &g1, const Group<boolean, int, float, double> &g2)
                       -> std::vector<std::array<int, 2>> {
                          // we use the tags defined above to find the top and antitop
                          // filter_XXX returns a vector of indices of elements fullfilling the criteria
                          // list of currently supported filter operations are in the group header file
                          auto top = g1.filter_equal("dileptonic_ttbar", 1);
                          auto antitop = g2.filter_equal("dileptonic_ttbar", 6);

//                          auto tops = g.sort_ascending("pdg", merge(g.filter_equal("pdg", 6), g.filter_equal("pdg", -6)));
//                          auto tops_22 = g.sort_ascending("pdg", merge(g.filter_equal("pdg", 6).filter_equal("status", 22), g.filter_equal("pdg", -6).filter_equal("status", 22)));
//                          std::cout << "tops.size" << tops.size() << std::endl;
//                          std::cout << "tops_22.size" << tops_22.size() << std::endl;
//                          float eta_top = g1.template get<float>("eta")[tops_22[0]]; 
//                          float eta_antitop = g2.template get<float>("eta")[tops_22[1]];
//                          float phi_top = g1.template get<float>("phi")[tops_22[0]];
//                          float phi_antitop = g2.template get<float>("phi")[tops_22[1]];
//                          float pt_top = g1.template get<float>("pt")[tops_22[0]];
//                          float pt_antitop = g2.template get<float>("pt")[tops_22[1]];
//                          float m_top = g1.template get<float>("default_mass")[tops_22[0]];
//                          float m_antitop = g2.template get<float>("default_mass")[tops_22[1]];
//                            
//                          TLorentzVector top_vec, antitop_vec;
//                          top_vec.SetPtEtaPhiM(pt_top, eta_top, phi_top, m_top);
//                          antitop_vec.SetPtEtaPhiM(pt_antitop, eta_antitop, phi_antitop, m_antitop);
//                          double_t top_X = top_vec.X();
//                          double_t top_Y = top_vec.Y();
//                          double_t top_Z = top_vec.Z();
//                          double_t top_T = top_vec.T();
//                          double_t antitop_X = antitop_vec.X();
//                          double_t antitop_Y = antitop_vec.Y();
//                          double_t antitop_Z = antitop_vec.Z();
//                          double_t antitop_T = antitop_vec.T();
                          
 //                         std::cout << "top_p4 = ( " << top_X << ", " << top_Y << ", " << top_Z << ", " << top_T << " )" << std::endl;  
 //                         std::cout << "antitop_p4 = ( " << antitop_X << ", " << antitop_Y << ", " << antitop_Z << ", " << antitop_T << " )"<< std::endl;  

                          // check that the collection contains exactly one top and one antitop
                          // if not the case, return an empty index list
                          if (top.size() != 1 or antitop.size() != 1){
//                            std::cout << "tops are empty in gen_ttbar" << std::endl;
                            return {};
                          }
//                          else { 
//                            int top_status = g1.template get<int>("status")[tops[0]];
//                            int antitop_status = g2.template get<int>("status")[tops[1]];
//                            std::cout << "tops are there, status top: " << top_status << ", antitop status: " << antitop_status << std::endl;
//                            int top_22 = g.filter_equal("pdg", 6).filter_equal("status", 22).size();
//                            int top_44 = g.filter_equal("pdg", 6).filter_equal("status", 44).size();
//                            int top_62 = g.filter_equal("pdg", 6).filter_equal("status", 62).size();
//                            int antitop_22 = g.filter_equal("pdg", -6).filter_equal("status", 22).size();
//                            int antitop_44 = g.filter_equal("pdg", -6).filter_equal("status", 44).size();
//                            int antitop_62 = g.filter_equal("pdg", -6).filter_equal("status", 62).size();
//                            std::cout << "top_22: " << top_22 << ", top_44: " << top_44 << ", top_62: " << top_62 << std::endl; 
//                            std::cout << "antitop_22: " << antitop_22 << ", antitop_44: " << antitop_44 << ", antitop_62: " << antitop_62 << std::endl; 
//                          }
                          return {{top[0], antitop[0]}};
                       });
  std::cout << "after gen_ttbar_indexer";

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
  
  // to make adding attributes shorter
  const auto add_attribute = [&](const std::string &funcname, attr_func_type_ttbar func) {
    gen_ttbar.add_attribute(funcname, func, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
   					    "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");  
  };


  gen_ttbar.add_attribute("ttbar_mass", system_invariant_mass, 
                          "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                          "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_ttbar.add_attribute("ttbar_pt", pt_vector_sum, 
                          "gen_particle::pt", "gen_particle::phi", "gen_particle::pt", "gen_particle::phi");

  gen_ttbar.add_attribute("ttbar_angle", my_angle, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
						   "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_ttbar.add_attribute("t_CosTheta", cos_theta_old, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
  						        "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_ttbar.add_attribute("cpTTT_gen_self_calc", cpTTT, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
					      "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  add_attribute("ttbar_rapidity", get_attr_function([](TLorentzVector sum){return sum.Rapidity();}));
  add_attribute("ttbar_beta", get_attr_function([](TLorentzVector sum){return sum.Beta();}));
  add_attribute("ttbar_pseudo_rapidity", get_attr_function([](TLorentzVector sum){return sum.PseudoRapidity();}));
  add_attribute("ttbar_gamma", get_attr_function([](TLorentzVector sum){return sum.Gamma();}));
  add_attribute("ttbar_e", get_attr_function([](TLorentzVector sum){return sum.E();}));
  add_attribute("ttbar_energy", get_attr_function([](TLorentzVector sum){return sum.Energy();}));
  add_attribute("ttbar_et", get_attr_function([](TLorentzVector sum){return sum.Et();}));
  add_attribute("ttbar_et2", get_attr_function([](TLorentzVector sum){return sum.Et2();}));
  add_attribute("ttbar_m2", get_attr_function([](TLorentzVector sum){return sum.M2();}));
  add_attribute("ttbar_mag", get_attr_function([](TLorentzVector sum){return sum.Mag();}));
  add_attribute("ttbar_mag2", get_attr_function([](TLorentzVector sum){return sum.Mag2();}));
  add_attribute("ttbar_minus", get_attr_function([](TLorentzVector sum){return sum.Minus();}));
  add_attribute("ttbar_mt", get_attr_function([](TLorentzVector sum){return sum.Mt();}));
  add_attribute("ttbar_mt2", get_attr_function([](TLorentzVector sum){return sum.Mt2();}));
  add_attribute("ttbar_p", get_attr_function([](TLorentzVector sum){return sum.P();}));
  add_attribute("ttbar_perp", get_attr_function([](TLorentzVector sum){return sum.Perp();}));
  add_attribute("ttbar_perp2", get_attr_function([](TLorentzVector sum){return sum.Perp2();}));
  add_attribute("ttbar_phi", get_attr_function([](TLorentzVector sum){return sum.Phi();}));
  add_attribute("ttbar_plus", get_attr_function([](TLorentzVector sum){return sum.Plus();}));
  add_attribute("ttbar_px", get_attr_function([](TLorentzVector sum){return sum.Px();}));
  add_attribute("ttbar_py", get_attr_function([](TLorentzVector sum){return sum.Py();}));
  add_attribute("ttbar_pz", get_attr_function([](TLorentzVector sum){return sum.Pz();}));
  add_attribute("ttbar_rho", get_attr_function([](TLorentzVector sum){return sum.Rho();}));
  add_attribute("ttbar_t", get_attr_function([](TLorentzVector sum){return sum.T();}));
  add_attribute("ttbar_theta", get_attr_function([](TLorentzVector sum){return sum.Theta();}));
  //add_attribute("ttbar_x", get_attr_function([](TLorentzVector sum){return sum.X();}));
  //add_attribute("ttbar_y", get_attr_function([](TLorentzVector sum){return sum.Y();}));
  //add_attribute("ttbar_z", get_attr_function([](TLorentzVector sum){return sum.Z();}));



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
  
  // t and tbar mass
  gen_ttbar.add_attribute("top_mass", return_first, "gen_particle::mass", "gen_particle::mass");
  gen_ttbar.add_attribute("antitop_mass", return_second, "gen_particle::mass", "gen_particle::mass");

  // t and tbar eta
  gen_ttbar.add_attribute("top_eta", return_first, "gen_particle::eta", "gen_particle::eta");
  gen_ttbar.add_attribute("antitop_eta", return_second, "gen_particle::eta", "gen_particle::eta");

  // and redundantly obtain the ttbar pt again
  gen_ttbar.transform_attribute("ttbar_pt_transform", pt_vector_sum, 
                                "top_pt", "top_phi", "antitop_pt", "antitop_phi");

  // now clearly the index masking examples above are rather artificial, serving only to highlight features
  // rather than something one might actually be interested in doing in an actual analysis
  // so now let us consider a more typical example of a dileptonic ttbar system
  // where we are interested in six particles: tops, charged leptons and bottoms
  Aggregate gen_tt_ll_bb("gen_tt_ll_bb", 22, 1, gen_particle, gen_particle, gen_particle, gen_particle, gen_particle, gen_particle);

  std::cout << "before gen_tt_ll_bb_indexer" << std::endl;

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

  std::cout << "after gen_tt_ll_bb_indexer" << std::endl;

  // having to define a masking for each case can be cumbersome
  // so we will be using the function_util plugin to simplify the work somewhat
  // there a function identity<type> is defined that returns a copy of its only argument
  // i.e. identity<float> do the same thing as return_first above, except it takes only one argument
  // furthermore, a function apply_to<N> is provided, that takes one argument that is a function
  // let ff be a function taking F arguments
  // then a call of apply_to<N>(ff) returns a function that takes (N + 1) * F arguments
  // whose result is the same as calling ff on the last F arguments, ignoring the first NF arguments
  // e.g. apply_to<1>(identity<float>) returns a function that is identical to return_second above
  gen_tt_ll_bb.add_attribute("lepton_pt", apply_to<2>(identity<>), 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("lepton_eta", apply_to<2>(identity<>), 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  gen_tt_ll_bb.add_attribute("antilepton_pt", apply_to<3>(identity<>), 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("antilepton_eta", apply_to<3>(identity<>), 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  gen_tt_ll_bb.add_attribute("bottom_pt", apply_to<4>(identity<>), 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("bottom_eta", apply_to<4>(identity<>), 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  gen_tt_ll_bb.add_attribute("antibottom_pt", apply_to<5>(identity<>), 
                             "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt", "gen_particle::pt");
  gen_tt_ll_bb.add_attribute("antibottom_eta", apply_to<5>(identity<>), 
                             "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta", "gen_particle::eta");

  // let's also consider some invariant masses
  gen_tt_ll_bb.add_attribute("ttbar_mass", system_invariant_mass, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("llbar_mass", apply_to<1>(system_invariant_mass), 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("bbbar_mass", apply_to<2>(system_invariant_mass), 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  // in the cases below where we can't divide the arguments into first NF and last F arguments
  // we need to use regular masking for the moment
  gen_tt_ll_bb.add_attribute("lbbar_mass", [] (float , float , float , float , 
                                               float , float , float , float ,
                                               float pt1, float eta1, float phi1, float m1,
                                               float , float , float , float , 
                                               float , float , float , float ,
                                               float pt2, float eta2, float phi2, float m2)
                             { return system_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
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
                             { return system_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
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
                             { return system_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
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
                             { return system_invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  
  gen_tt_ll_bb.add_attribute("llbarbbbar_mass", [] (float , float , float , float , 
                                                  float , float , float , float ,
                                                  float pt1, float eta1, float phi1, float m1,
                                                  float pt2, float eta2, float phi2, float m2,
                                                  float pt3, float eta3, float phi3, float m3,
                                                  float pt4, float eta4, float phi4, float m4)
                             { return system_invariant_mass_four_particles(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2, pt3, eta3, phi3, m3, pt4, eta4, phi4, m4); }, 
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");

  gen_tt_ll_bb.add_attribute("llbar_phi", llbar_phi,
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");
  
  gen_tt_ll_bb.add_attribute("llbar_phi_diff", llbar_dphi,
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
                             "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");


  // adding spin correlation 
  const auto spin_add_attribute = [&](const std::string &funcname, attr_func_type_spin func) {
    gen_tt_ll_bb.add_attribute(funcname, func, "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",
					    "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",  
   					    "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass",  
   					    "gen_particle::pt", "gen_particle::eta", "gen_particle::phi", "gen_particle::mass");  
  };

  spin_add_attribute("cHel", spin_correlation<>("cHel"));
  spin_add_attribute("cHel", spin_correlation<>("cHel"));
  spin_add_attribute("cLab", spin_correlation<>("cLab"));
  spin_add_attribute("ckk", spin_correlation<>("ckk"));
  spin_add_attribute("crr", spin_correlation<>("crr"));
  spin_add_attribute("cnn", spin_correlation<>("cnn"));
  spin_add_attribute("cHan", spin_correlation<>("cHan"));
  spin_add_attribute("cSca", spin_correlation<>("cSca"));
  spin_add_attribute("cTra", spin_correlation<>("cTra"));
  spin_add_attribute("phi0", spin_correlation<>("phi0"));
  spin_add_attribute("phi1", spin_correlation<>("phi1"));
  spin_add_attribute("cpTTT", spin_correlation<>("cpTTT"));
  spin_add_attribute("cpTP", spin_correlation<>("cpTP"));
  
  spin_add_attribute("dPhi", spin_correlation<>("dPhi"));
  spin_add_attribute("dEta", spin_correlation<>("dEta"));

  spin_add_attribute("kdx", spin_correlation<>("kdx"));
  spin_add_attribute("kdy", spin_correlation<>("kdy"));
  spin_add_attribute("kdz", spin_correlation<>("kdz"));

  spin_add_attribute("rdx", spin_correlation<>("rdx"));
  spin_add_attribute("rdy", spin_correlation<>("rdy"));
  spin_add_attribute("rdz", spin_correlation<>("rdz"));

  spin_add_attribute("ndx", spin_correlation<>("ndx"));
  spin_add_attribute("ndy", spin_correlation<>("ndy"));
  spin_add_attribute("ndz", spin_correlation<>("ndz"));

  spin_add_attribute("b1k", spin_correlation<>("b1k"));
  spin_add_attribute("b2k", spin_correlation<>("b2k"));

  spin_add_attribute("b1j", spin_correlation<>("b1j"));
  spin_add_attribute("b2j", spin_correlation<>("b2j"));

  spin_add_attribute("b1r", spin_correlation<>("b1r"));
  spin_add_attribute("b2r", spin_correlation<>("b2r"));

  spin_add_attribute("b1q", spin_correlation<>("b1q"));
  spin_add_attribute("b2q", spin_correlation<>("b2q"));

  spin_add_attribute("b1n", spin_correlation<>("b1n"));
  spin_add_attribute("b2n", spin_correlation<>("b2n"));

  spin_add_attribute("b1x", spin_correlation<>("b1x"));
  spin_add_attribute("b2x", spin_correlation<>("b2x"));

  spin_add_attribute("b1y", spin_correlation<>("b1y"));
  spin_add_attribute("b2y", spin_correlation<>("b2y"));

  spin_add_attribute("b1z", spin_correlation<>("b1z"));
  spin_add_attribute("b2z", spin_correlation<>("b2z"));

  spin_add_attribute("bPkk", spin_correlation<>("bPkk"));
  spin_add_attribute("bMkk", spin_correlation<>("bMkk"));

  spin_add_attribute("bPjj", spin_correlation<>("bPjj"));
  spin_add_attribute("bMjj", spin_correlation<>("bMjj"));

  spin_add_attribute("bPrr", spin_correlation<>("bPrr"));
  spin_add_attribute("bMrr", spin_correlation<>("bMrr"));

  spin_add_attribute("bPqq", spin_correlation<>("bPqq"));
  spin_add_attribute("bMqq", spin_correlation<>("bMqq"));

  spin_add_attribute("bPnn", spin_correlation<>("bPnn"));
  spin_add_attribute("bMnn", spin_correlation<>("bMnn"));

  spin_add_attribute("bPxx", spin_correlation<>("bPxx"));
  spin_add_attribute("bMxx", spin_correlation<>("bMxx"));

  spin_add_attribute("bPyy", spin_correlation<>("bPyy"));
  spin_add_attribute("bMyy", spin_correlation<>("bMyy"));

  spin_add_attribute("bPzz", spin_correlation<>("bPzz"));
  spin_add_attribute("bMzz", spin_correlation<>("bMzz"));

  spin_add_attribute("crk", spin_correlation<>("crk"));
  spin_add_attribute("ckr", spin_correlation<>("ckr"));

  spin_add_attribute("cnr", spin_correlation<>("cnr"));
  spin_add_attribute("crn", spin_correlation<>("crn"));

  spin_add_attribute("cnk", spin_correlation<>("cnk"));
  spin_add_attribute("ckn", spin_correlation<>("ckn"));

  spin_add_attribute("cPrk", spin_correlation<>("cPrk"));
  spin_add_attribute("cMrk", spin_correlation<>("cMrk"));

  spin_add_attribute("cPnr", spin_correlation<>("cPnr"));
  spin_add_attribute("cMnr", spin_correlation<>("cMnr"));

  spin_add_attribute("cPnk", spin_correlation<>("cPnk"));
  spin_add_attribute("cMnk", spin_correlation<>("cMnk"));

  spin_add_attribute("cxx", spin_correlation<>("cxx"));
  spin_add_attribute("cyy", spin_correlation<>("cyy"));
  spin_add_attribute("czz", spin_correlation<>("czz"));

  spin_add_attribute("cyx", spin_correlation<>("cyx"));
  spin_add_attribute("cxy", spin_correlation<>("cxy"));

  spin_add_attribute("czy", spin_correlation<>("czy"));
  spin_add_attribute("cyz", spin_correlation<>("cyz"));

  spin_add_attribute("czx", spin_correlation<>("czx"));
  spin_add_attribute("cxz", spin_correlation<>("cxz"));

  spin_add_attribute("cPyx", spin_correlation<>("cPyx"));
  spin_add_attribute("cMyx", spin_correlation<>("cMyx"));

  spin_add_attribute("cPzy", spin_correlation<>("cPzy"));
  spin_add_attribute("cMzy", spin_correlation<>("cMzy"));

  spin_add_attribute("cPzx", spin_correlation<>("cPzx"));
  spin_add_attribute("cMzx", spin_correlation<>("cMzx"));

  spin_add_attribute("crkP", spin_correlation<>("crkP"));
  spin_add_attribute("cnrP", spin_correlation<>("cnrP"));
  spin_add_attribute("cnkP", spin_correlation<>("cnkP"));

  spin_add_attribute("crkM", spin_correlation<>("crkM"));
  spin_add_attribute("cnrM", spin_correlation<>("cnrM"));
  spin_add_attribute("cnkM", spin_correlation<>("cnkM"));

  spin_add_attribute("cXxx", spin_correlation<>("cXxx"));
  spin_add_attribute("cYyy", spin_correlation<>("cYyy"));
  spin_add_attribute("cZzz", spin_correlation<>("cZzz"));

  spin_add_attribute("cyxP", spin_correlation<>("cyxP"));
  spin_add_attribute("czyP", spin_correlation<>("czyP"));
  spin_add_attribute("czxP", spin_correlation<>("czxP"));

  spin_add_attribute("cyxM", spin_correlation<>("cyxM"));
  spin_add_attribute("czyM", spin_correlation<>("czyM"));
  spin_add_attribute("czxM", spin_correlation<>("czxM"));

  spin_add_attribute("kNorm", spin_correlation<>("kNorm"));
  spin_add_attribute("rNorm", spin_correlation<>("rNorm"));
  spin_add_attribute("nNorm", spin_correlation<>("nNorm"));

  std::cout << "after all add attributes" << std::endl;


  // let's histogram the attributes we defined above
  // this is done through the histogram class, which handles a group of histograms sharing the same weights and to be filled at the same time
  // in this example we will have two instances, before and after some acceptance cuts
  Histogram hist_no_cut;

  // we define an argument-less impure function that computes the weight for each histogram entry
  // here we only take the per-event weight from the metadata collection
  // we can see here that internally a non-array collection is in fact an array collection of size 1
  // if no weighter is defined, histograms are filled with weight 1
  //hist_no_cut.set_weighter([&weight = metadata.get<float>("weight")] () { return weight[0]; });

  hist_no_cut.set_weighter([&weight = metadata.get<float>("weight"), &reweighting_weight = lhe_event.get<float>("weight_float")] () { return (weight[0] * reweighting_weight[0]); });
  //hist_no_cut.set_weighter([&reweighting_weight = lhe_event.get<float>("weight_float")] () { return (reweighting_weight[0]); });
  //hist_no_cut.set_weighter([&weight = lhe_event.get<float>("weight_float")] () { return (weight[0]); });

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
  
  // z aus Juans Paper
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "t_CosTheta"), "t_CosTheta_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "cpTTT_gen_self_calc"), "cpTTT_gen_self_calc_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE"), "cpTTT_LHE_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE_self_calc"), "cpTTT_LHE_self_calc_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTP_LHE"), "cpTP_LHE_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "t_CosTheta_LHE"), "t_CosTheta_LHE_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "diff_cos_theta_cpTP_LHE"), "diff_cos_theta_LHE_no_cut", "", 100, -1.f, 1.f);
  
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
  // war vorher bis 1500
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ttbar_mass"), "ttbar_mass_2_no_cut", "", 120, 300.f, 600.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_mass"), "llbar_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bbbar_mass"), "bbbar_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lb_mass"), "lb_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarbbar_mass"), "lbarbbar_mass_no_cut", "", 120, 0.f, 1200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbbar_mass"), "lbbar_mass_no_cut", "", 100, 0.f, 200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarb_mass"), "lbarb_mass_no_cut", "", 100, 0.f, 200.f);


  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi"), "llbar_phi_no_cut", "", 100, 0.f, 5.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi_diff"), "llbar_phi_diff_no_cut", "", 100, 0.f, 5.f);
  
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_mass"), "top_mass_no_cut", "", 100, 150.f, 200.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_mass"), "antitop_mass_no_cut", "", 100, 150.f, 200.f);

  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_phi"), "top_phi_no_cut", "", 160, 0.f, 5.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_phi"), "antitop_phi_no_cut", "", 160, 0.f, 5.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_eta"), "top_eta_no_cut", "", 100, -5.f, 5.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_eta"), "antitop_eta_no_cut", "", 100, -5.f, 5.f);

  // ttbar things
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rapidity"), "ttbar_rapidity_no_cut", "", 100, -4.f, 4.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pseudo_rapidity"), "ttbar_pseudo_rapidity_no_cut", "", 100, -10.f, 10.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_angle"), "ttbar_angle_no_cut", "", 100, -1.f, 4.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_beta"), "ttbar_beta_no_cut", "", 100, 0.f, 1.f);
  // looks like this is pseudo rapidity as well
  //hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_eta"), "ttbar_eta_no_cut", "", 100, -10.f, 10.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_gamma"), "ttbar_gamma_no_cut", "", 100, 0.f, 10.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_e"), "ttbar_gamma_no_cut", "", 100, 0.f, 10.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_energy"), "ttbar_gamma_no_cut", "", 100, 0.f, 10.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et"), "ttbar_et_no_cut", "", 100, -1.f, 100.f);
  // looks strange
  //hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et2"), "ttbar_et2_no_cut", "", 100, 0.f, 100.f);
  // m2, mag, mag2, mt, mt2 and t are always zero
  // mag(2) and m(2) are the same
  // I dont know what that is
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_p"), "ttbar_p_no_cut", "", 100, 0.f, 500.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp"), "ttbar_perp_no_cut", "", 100, 0.f, 100.f);
  //hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp2"), "ttbar_perp2_no_cut", "", 100, 0.f, 100.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_phi"), "ttbar_phi_no_cut", "", 100, -4.f, 4.f);
  // plus and minus are T +/- Z and has sth to do with light cones in relativity
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_plus"), "ttbar_plus_no_cut", "", 100, 0.f, 100.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_minus"), "ttbar_minus_no_cut", "", 100, 0.f, 500.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_px"), "ttbar_px_no_cut", "", 100, -100.f, 100.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_py"), "ttbar_py_no_cut", "", 100, -100.f, 100.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pz"), "ttbar_pz_no_cut", "", 100, -100.f, 100.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rho"), "ttbar_rho_no_cut", "", 100, 0.f, 100.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_theta"), "ttbar_theta_no_cut", "", 100, -1.f, 4.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_x"), "ttbar_x_no_cut", "", 100, -100.f, 100.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_y"), "ttbar_y_no_cut", "", 100, -100.f, 100.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_z"), "ttbar_z_no_cut", "", 100, -100.f, 100.f);


  //spin correlation
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHel"), "spin_cHel_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cLab"), "spin_cLab_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crr"), "spin_crr_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnn"), "spin_cnn_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckk"), "spin_ckk_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cSca"), "spin_cSca_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cTra"), "spin_cTra_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi0"), "spin_phi0_no_cut", "", 100, 0.f, 4.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi1"), "spin_phi1_no_cut", "", 100, 0.f, 7.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTP"), "spin_cpTP_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTTT"), "spin_cpTTT_no_cut", "", 100, -1.f, 1.f);
  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHan"), "spin_cHan_no_cut", "", 100, -1.f, 1.f);
  
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "dPhi"), "spin_dPhi", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "dEta"), "spin_dEta", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kdx"), "spin_kdx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kdy"), "spin_kdy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kdz"), "spin_kdz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rdx"), "spin_rdx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rdy"), "spin_rdy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rdz"), "spin_rdz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ndx"), "spin_ndx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ndy"), "spin_ndy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ndz"), "spin_ndz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1k"), "spin_b1k", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2k"), "spin_b2k", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1j"), "spin_b1j", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2j"), "spin_b2j", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1r"), "spin_b1r", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2r"), "spin_b2r", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1q"), "spin_b1q", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2q"), "spin_b2q", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1n"), "spin_b1n", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2n"), "spin_b2n", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1x"), "spin_b1x", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2x"), "spin_b2x", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1y"), "spin_b1y", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2y"), "spin_b2y", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1z"), "spin_b1z", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2z"), "spin_b2z", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPkk"), "spin_bPkk", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMkk"), "spin_bMkk", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPjj"), "spin_bPjj", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMjj"), "spin_bMjj", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPrr"), "spin_bPrr", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMrr"), "spin_bMrr", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPqq"), "spin_bPqq", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMqq"), "spin_bMqq", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPnn"), "spin_bPnn", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMnn"), "spin_bMnn", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPxx"), "spin_bPxx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMxx"), "spin_bMxx", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPyy"), "spin_bPyy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMyy"), "spin_bMyy", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPzz"), "spin_bPzz", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMzz"), "spin_bMzz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crk"), "spin_crk", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckr"), "spin_ckr", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnr"), "spin_cnr", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crn"), "spin_crn", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnk"), "spin_cnk", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckn"), "spin_ckn", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPrk"), "spin_cPrk", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMrk"), "spin_cMrk", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPnr"), "spin_cPnr", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMnr"), "spin_cMnr", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPnk"), "spin_cPnk", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMnk"), "spin_cMnk", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cxx"), "spin_cxx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyy"), "spin_cyy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czz"), "spin_czz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyx"), "spin_cyx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cxy"), "spin_cxy", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czy"), "spin_czy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyz"), "spin_cyz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czx"), "spin_czx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cxz"), "spin_cxz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPyx"), "spin_cPyx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMyx"), "spin_cMyx", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPzy"), "spin_cPzy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMzy"), "spin_cMzy", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPzx"), "spin_cPzx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMzx"), "spin_cMzx", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crkP"), "spin_crkP", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnrP"), "spin_cnrP", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnkP"), "spin_cnkP", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crkM"), "spin_crkM", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnrM"), "spin_cnrM", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnkM"), "spin_cnkM", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cXxx"), "spin_cXxx", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cYyy"), "spin_cYyy", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cZzz"), "spin_cZzz", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyxP"), "spin_cyxP", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czyP"), "spin_czyP", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czxP"), "spin_czxP", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyxM"), "spin_cyxM", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czyM"), "spin_czyM", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czxM"), "spin_czxM", "", 100, -1.f, 1.f);
//
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kNorm"), "spin_kNorm", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rNorm"), "spin_rNorm", "", 100, -1.f, 1.f);
//  hist_no_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "nNorm"), "spin_nNorm", "", 100, -1.f, 1.f);
  
  // let's define another histogram instance but now with acceptance cuts
  // we can, but don't need to, define the cuts in the filling function themselves
  // so the histogram instance is defined identically as above except the histogram names
  Histogram hist_cut;
  //hist_cut.set_weighter([&weight = metadata.get<float>("weight")] () { return weight[0]; });
  hist_cut.set_weighter([&weight = metadata.get<float>("weight"), &reweighting_weight = lhe_event.get<float>("weight_float")] () { return (weight[0] * reweighting_weight[0]); });
  //hist_cut.set_weighter([&reweighting_weight = lhe_event.get<float>("weight_float")] () { return (reweighting_weight[0]); });
  //hist_cut.set_weighter([&weight = lhe_event.get<float>("weight_float")] () { return (weight[0]); });
  
  // z aus Juans Paper
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "t_CosTheta"), "t_CosTheta_no_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "cpTTT_gen_self_calc"), "cpTTT_gen_self_calc_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE"), "cpTTT_lhe_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE_self_calc"), "cpTTT_lhe_cut_self_calc", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTP_LHE"), "cpTP_lhe_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(lhe_event, "t_CosTheta_LHE"), "t_CosTheta_lhe_no_cut", "", 100, -1.f, 1.f);

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

  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi"), "llbar_phi_cut", "", 100, 0.f, 5.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi_diff"), "llbar_phi_diff_cut", "", 100, 0.f, 5.f);
  
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_mass"), "top_mass_cut", "", 100, 150.f, 200.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_mass"), "antitop_mass_cut", "", 100, 150.f, 200.f);

  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_phi"), "top_phi_cut", "", 160, 0.f, 5.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_phi"), "antitop_phi_cut", "", 160, 0.f, 5.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_eta"), "top_eta_cut", "", 100, -5.f, 5.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_eta"), "antitop_eta_cut", "", 100, -5.f, 5.f);
 
  // ttbar things
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rapidity"), "ttbar_rapidity_cut", "", 100, -4.f, 4.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_angle"), "ttbar_angle_cut", "", 100, -1.f, 4.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_beta"), "ttbar_beta_cut", "", 100, 0.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pseudo_rapidity"), "ttbar_pseudo_rapidity_cut", "", 100, -10.f, 10.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_gamma"), "ttbar_gamma_cut", "", 100, 0.f, 10.f);
  
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_e"), "ttbar_gamma_cut", "", 100, 0.f, 10.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_energy"), "ttbar_gamma_cut", "", 100, 0.f, 10.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et"), "ttbar_et_cut", "", 100, -1.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et2"), "ttbar_et2_cut", "", 100, 0.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_minus"), "ttbar_minus_cut", "", 100, 0.f, 500.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_p"), "ttbar_p_cut", "", 100, 0.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp"), "ttbar_perp_cut", "", 100, 0.f, 100.f);
  //hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp2"), "ttbar_perp2_cut", "", 100, 0.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_phi"), "ttbar_phi_cut", "", 100, -4.f, 4.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_plus"), "ttbar_plus_cut", "", 100, 0.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_px"), "ttbar_px_cut", "", 100, -100.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_py"), "ttbar_py_cut", "", 100, -100.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pz"), "ttbar_pz_cut", "", 100, -100.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rho"), "ttbar_rho_cut", "", 100, 0.f, 100.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_theta"), "ttbar_theta_cut", "", 100, -1.f, 4.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_x"), "ttbar_x_cut", "", 100, -100.f, 100.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_y"), "ttbar_y_cut", "", 100, -100.f, 100.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_z"), "ttbar_z_cut", "", 100, -100.f, 100.f);
 

  //spin correlation
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHel"), "spin_cHel_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cLab"), "spin_cLab_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crr"), "spin_crr_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnn"), "spin_cnn_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckk"), "spin_ckk_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cSca"), "spin_cSca_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cTra"), "spin_cTra_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi0"), "spin_phi0_cut", "", 100, 0.f, 4.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi1"), "spin_phi1_cut", "", 100, 0.f, 7.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTP"), "spin_cpTp_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTTT"), "spin_cpTTT_cut", "", 100, -1.f, 1.f);
  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHan"), "spin_cHan_cut", "", 100, -1.f, 1.f);


//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "dPhi"), "spin_dPhi", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "dEta"), "spin_dEta", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kdx"), "spin_kdx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kdy"), "spin_kdy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kdz"), "spin_kdz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rdx"), "spin_rdx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rdy"), "spin_rdy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rdz"), "spin_rdz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ndx"), "spin_ndx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ndy"), "spin_ndy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ndz"), "spin_ndz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1k"), "spin_b1k", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2k"), "spin_b2k", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1j"), "spin_b1j", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2j"), "spin_b2j", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1r"), "spin_b1r", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2r"), "spin_b2r", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1q"), "spin_b1q", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2q"), "spin_b2q", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1n"), "spin_b1n", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2n"), "spin_b2n", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1x"), "spin_b1x", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2x"), "spin_b2x", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1y"), "spin_b1y", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2y"), "spin_b2y", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b1z"), "spin_b1z", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "b2z"), "spin_b2z", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPkk"), "spin_bPkk", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMkk"), "spin_bMkk", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPjj"), "spin_bPjj", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMjj"), "spin_bMjj", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPrr"), "spin_bPrr", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMrr"), "spin_bMrr", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPqq"), "spin_bPqq", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMqq"), "spin_bMqq", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPnn"), "spin_bPnn", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMnn"), "spin_bMnn", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPxx"), "spin_bPxx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMxx"), "spin_bMxx", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPyy"), "spin_bPyy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMyy"), "spin_bMyy", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bPzz"), "spin_bPzz", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bMzz"), "spin_bMzz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crk"), "spin_crk", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckr"), "spin_ckr", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnr"), "spin_cnr", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crn"), "spin_crn", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnk"), "spin_cnk", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckn"), "spin_ckn", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPrk"), "spin_cPrk", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMrk"), "spin_cMrk", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPnr"), "spin_cPnr", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMnr"), "spin_cMnr", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPnk"), "spin_cPnk", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMnk"), "spin_cMnk", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cxx"), "spin_cxx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyy"), "spin_cyy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czz"), "spin_czz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyx"), "spin_cyx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cxy"), "spin_cxy", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czy"), "spin_czy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyz"), "spin_cyz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czx"), "spin_czx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cxz"), "spin_cxz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPyx"), "spin_cPyx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMyx"), "spin_cMyx", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPzy"), "spin_cPzy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMzy"), "spin_cMzy", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cPzx"), "spin_cPzx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cMzx"), "spin_cMzx", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crkP"), "spin_crkP", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnrP"), "spin_cnrP", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnkP"), "spin_cnkP", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crkM"), "spin_crkM", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnrM"), "spin_cnrM", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnkM"), "spin_cnkM", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cXxx"), "spin_cXxx", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cYyy"), "spin_cYyy", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cZzz"), "spin_cZzz", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyxP"), "spin_cyxP", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czyP"), "spin_czyP", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czxP"), "spin_czxP", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cyxM"), "spin_cyxM", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czyM"), "spin_czyM", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "czxM"), "spin_czxM", "", 100, -1.f, 1.f);
//
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "kNorm"), "spin_kNorm", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "rNorm"), "spin_rNorm", "", 100, -1.f, 1.f);
//  hist_cut.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "nNorm"), "spin_nNorm", "", 100, -1.f, 1.f);

  Histogram hist_positive;

  hist_positive.set_weighter([&weight = metadata.get<float>("weight"), &reweighting_weight = lhe_event.get<float>("weight_float")] () { return (weight[0] * reweighting_weight[0]); });
  //hist_positive.set_weighter([&reweighting_weight = lhe_event.get<float>("weight_float")] () { return (reweighting_weight[0]); });
  //hist_positive.set_weighter([&weight = lhe_event.get<float>("weight_float")] () { return (weight[0]); });

  // next we define the histograms, where the histogram type are given inside the <> bracket
  // all histogram types supported by ROOT are supported
  // the first argument is an impure function instructing how the histograms should be filled
  // the function takes two arguments, a histogram pointer and a weight
  // the remaining arguments are those expected by ROOT histogram constructor of the type being used
  hist_positive.make_histogram<TH1F>([&gen_ttbar] (TH1F *hist, double weight) {
      // do not fill if the aggregate has no elements
      if (gen_ttbar.n_elements() != 1)
        return;

      auto &mass = gen_ttbar.get<float>("ttbar_mass");
      hist->Fill(mass[0], weight);
    }, "ttbar_mass_positive", "", 120, 300.f, 1500.f);

  // in many cases we will be filling the histograms in similar ways
  // e.g. check for presence, and if yes, fill the first/all elements
  // and having to write out the filling function every time can be cumbersome
  // so in the plugins some utility functions are provided for these commonly used functions
  
  // z aus Juans Paper
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "t_CosTheta"), "t_CosTheta_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "cpTTT_gen_self_calc"), "cpTTT_gen_self_calc_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE"), "cpTTT_LHE_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE_self_calc"), "cpTTT_LHE_self_calc_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTP_LHE"), "cpTP_LHE_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(lhe_event, "t_CosTheta_LHE"), "t_CosTheta_LHE_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(lhe_event, "diff_cos_theta_cpTP_LHE"), "diff_cos_theta_LHE_positive", "", 100, -1.f, 1.f);
  
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt"), "ttbar_pt_positive", "", 120, 0.f, 1200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt_transform"), "ttbar_pt_transform_positive", "", 120, 0.f, 1200.f);

  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_pt"), "lepton_pt_positive", "", 100, 0.f, 400.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_eta"), "lepton_eta_positive", "", 100, -5.f, 5.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_pt"), "antilepton_pt_positive", "", 100, 0.f, 400.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_eta"), "antilepton_eta_positive", "", 100, -5.f, 5.f);

  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_pt"), "bottom_pt_positive", "", 100, 0.f, 400.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_eta"), "bottom_eta_positive", "", 100, -5.f, 5.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_pt"), "antibottom_pt_positive", "", 100, 0.f, 400.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_eta"), "antibottom_eta_positive", "", 100, -5.f, 5.f);
  // war vorher bis 1500
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ttbar_mass"), "ttbar_mass_2_positive", "", 120, 300.f, 600.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_mass"), "llbar_mass_positive", "", 120, 0.f, 1200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bbbar_mass"), "bbbar_mass_positive", "", 120, 0.f, 1200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lb_mass"), "lb_mass_positive", "", 120, 0.f, 1200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarbbar_mass"), "lbarbbar_mass_positive", "", 120, 0.f, 1200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbbar_mass"), "lbbar_mass_positive", "", 100, 0.f, 200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarb_mass"), "lbarb_mass_positive", "", 100, 0.f, 200.f);


  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi"), "llbar_phi_positive", "", 100, 0.f, 5.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi_diff"), "llbar_phi_diff_positive", "", 100, 0.f, 5.f);
  
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_mass"), "top_mass_positive", "", 100, 150.f, 200.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_mass"), "antitop_mass_positive", "", 100, 150.f, 200.f);

  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_phi"), "top_phi_positive", "", 160, 0.f, 5.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_phi"), "antitop_phi_positive", "", 160, 0.f, 5.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_eta"), "top_eta_positive", "", 100, -5.f, 5.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_eta"), "antitop_eta_positive", "", 100, -5.f, 5.f);

  // ttbar things
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rapidity"), "ttbar_rapidity_positive", "", 100, -4.f, 4.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pseudo_rapidity"), "ttbar_pseudo_rapidity_positive", "", 100, -10.f, 10.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_angle"), "ttbar_angle_positive", "", 100, -1.f, 4.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_beta"), "ttbar_beta_positive", "", 100, 0.f, 1.f);
  // looks like this is pseudo rapidity as well
  //hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_eta"), "ttbar_eta_positive", "", 100, -10.f, 10.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_gamma"), "ttbar_gamma_positive", "", 100, 0.f, 10.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_e"), "ttbar_gamma_positive", "", 100, 0.f, 10.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_energy"), "ttbar_gamma_positive", "", 100, 0.f, 10.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et"), "ttbar_et_positive", "", 100, -1.f, 100.f);
  // looks strange
  //hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et2"), "ttbar_et2_positive", "", 100, 0.f, 100.f);
  // m2, mag, mag2, mt, mt2 and t are always zero
  // mag(2) and m(2) are the same
  // I dont know what that is
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_p"), "ttbar_p_positive", "", 100, 0.f, 500.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp"), "ttbar_perp_positive", "", 100, 0.f, 100.f);
  //hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp2"), "ttbar_perp2_positive", "", 100, 0.f, 100.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_phi"), "ttbar_phi_positive", "", 100, -4.f, 4.f);
  // plus and minus are T +/- Z and has sth to do with light cones in relativity
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_plus"), "ttbar_plus_positive", "", 100, 0.f, 100.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_minus"), "ttbar_minus_positive", "", 100, 0.f, 500.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_px"), "ttbar_px_positive", "", 100, -100.f, 100.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_py"), "ttbar_py_positive", "", 100, -100.f, 100.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pz"), "ttbar_pz_positive", "", 100, -100.f, 100.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rho"), "ttbar_rho_positive", "", 100, 0.f, 100.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_theta"), "ttbar_theta_positive", "", 100, -1.f, 4.f);
//  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_x"), "ttbar_x_positive", "", 100, -100.f, 100.f);
//  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_y"), "ttbar_y_positive", "", 100, -100.f, 100.f);
//  hist_positive.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_z"), "ttbar_z_positive", "", 100, -100.f, 100.f);


  //spin correlation
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHel"), "spin_cHel_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cLab"), "spin_cLab_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crr"), "spin_crr_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnn"), "spin_cnn_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckk"), "spin_ckk_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cSca"), "spin_cSca_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cTra"), "spin_cTra_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi0"), "spin_phi0_positive", "", 100, 0.f, 4.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi1"), "spin_phi1_positive", "", 100, 0.f, 7.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTP"), "spin_cpTP_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTTT"), "spin_cpTTT_positive", "", 100, -1.f, 1.f);
  hist_positive.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHan"), "spin_cHan_positive", "", 100, -1.f, 1.f);



  Histogram hist_negative;

  hist_negative.set_weighter([&weight = metadata.get<float>("weight"), &reweighting_weight = lhe_event.get<float>("weight_float")] () { return ( - weight[0] * reweighting_weight[0]); });
  //hist_negative.set_weighter([&reweighting_weight = lhe_event.get<float>("weight_float")] () { return (reweighting_weight[0]); });
  //hist_negative.set_weighter([&weight = lhe_event.get<float>("weight_float")] () { return (weight[0]); });

  // next we define the histograms, where the histogram type are given inside the <> bracket
  // all histogram types supported by ROOT are supported
  // the first argument is an impure function instructing how the histograms should be filled
  // the function takes two arguments, a histogram pointer and a weight
  // the remaining arguments are those expected by ROOT histogram constructor of the type being used
  hist_negative.make_histogram<TH1F>([&gen_ttbar] (TH1F *hist, double weight) {
      // do not fill if the aggregate has no elements
      if (gen_ttbar.n_elements() != 1)
        return;

      auto &mass = gen_ttbar.get<float>("ttbar_mass");
      hist->Fill(mass[0], weight);
    }, "ttbar_mass_negative", "", 120, 300.f, 1500.f);

  // in many cases we will be filling the histograms in similar ways
  // e.g. check for presence, and if yes, fill the first/all elements
  // and having to write out the filling function every time can be cumbersome
  // so in the plugins some utility functions are provided for these commonly used functions
  
  // z aus Juans Paper
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "t_CosTheta"), "t_CosTheta_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "cpTTT_gen_self_calc"), "cpTTT_gen_self_calc_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE"), "cpTTT_LHE_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTTT_LHE_self_calc"), "cpTTT_LHE_self_calc_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(lhe_event, "cpTP_LHE"), "cpTP_LHE_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(lhe_event, "t_CosTheta_LHE"), "t_CosTheta_LHE_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(lhe_event, "diff_cos_theta_cpTP_LHE"), "diff_cos_theta_LHE_negative", "", 100, -1.f, 1.f);
  
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt"), "ttbar_pt_negative", "", 120, 0.f, 1200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pt_transform"), "ttbar_pt_transform_negative", "", 120, 0.f, 1200.f);

  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_pt"), "lepton_pt_negative", "", 100, 0.f, 400.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lepton_eta"), "lepton_eta_negative", "", 100, -5.f, 5.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_pt"), "antilepton_pt_negative", "", 100, 0.f, 400.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antilepton_eta"), "antilepton_eta_negative", "", 100, -5.f, 5.f);

  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_pt"), "bottom_pt_negative", "", 100, 0.f, 400.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bottom_eta"), "bottom_eta_negative", "", 100, -5.f, 5.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_pt"), "antibottom_pt_negative", "", 100, 0.f, 400.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "antibottom_eta"), "antibottom_eta_negative", "", 100, -5.f, 5.f);
  // war vorher bis 1500
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ttbar_mass"), "ttbar_mass_2_negative", "", 120, 300.f, 600.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_mass"), "llbar_mass_negative", "", 120, 0.f, 1200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "bbbar_mass"), "bbbar_mass_negative", "", 120, 0.f, 1200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lb_mass"), "lb_mass_negative", "", 120, 0.f, 1200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarbbar_mass"), "lbarbbar_mass_negative", "", 120, 0.f, 1200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbbar_mass"), "lbbar_mass_negative", "", 100, 0.f, 200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "lbarb_mass"), "lbarb_mass_negative", "", 100, 0.f, 200.f);


  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi"), "llbar_phi_negative", "", 100, 0.f, 5.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "llbar_phi_diff"), "llbar_phi_diff_negative", "", 100, 0.f, 5.f);
  
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_mass"), "top_mass_negative", "", 100, 150.f, 200.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_mass"), "antitop_mass_negative", "", 100, 150.f, 200.f);

  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_phi"), "top_phi_negative", "", 160, 0.f, 5.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_phi"), "antitop_phi_negative", "", 160, 0.f, 5.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "top_eta"), "top_eta_negative", "", 100, -5.f, 5.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "antitop_eta"), "antitop_eta_negative", "", 100, -5.f, 5.f);

  // ttbar things
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rapidity"), "ttbar_rapidity_negative", "", 100, -4.f, 4.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pseudo_rapidity"), "ttbar_pseudo_rapidity_negative", "", 100, -10.f, 10.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_angle"), "ttbar_angle_negative", "", 100, -1.f, 4.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_beta"), "ttbar_beta_negative", "", 100, 0.f, 1.f);
  // looks like this is pseudo rapidity as well
  //hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_eta"), "ttbar_eta_negative", "", 100, -10.f, 10.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_gamma"), "ttbar_gamma_negative", "", 100, 0.f, 10.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_e"), "ttbar_gamma_negative", "", 100, 0.f, 10.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_energy"), "ttbar_gamma_negative", "", 100, 0.f, 10.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et"), "ttbar_et_negative", "", 100, -1.f, 100.f);
  // looks strange
  //hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_et2"), "ttbar_et2_negative", "", 100, 0.f, 100.f);
  // m2, mag, mag2, mt, mt2 and t are always zero
  // mag(2) and m(2) are the same
  // I dont know what that is
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_p"), "ttbar_p_negative", "", 100, 0.f, 500.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp"), "ttbar_perp_negative", "", 100, 0.f, 100.f);
  //hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_perp2"), "ttbar_perp2_negative", "", 100, 0.f, 100.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_phi"), "ttbar_phi_negative", "", 100, -4.f, 4.f);
  // plus and minus are T +/- Z and has sth to do with light cones in relativity
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_plus"), "ttbar_plus_negative", "", 100, 0.f, 100.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_minus"), "ttbar_minus_negative", "", 100, 0.f, 500.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_px"), "ttbar_px_negative", "", 100, -100.f, 100.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_py"), "ttbar_py_negative", "", 100, -100.f, 100.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_pz"), "ttbar_pz_negative", "", 100, -100.f, 100.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_rho"), "ttbar_rho_negative", "", 100, 0.f, 100.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_theta"), "ttbar_theta_negative", "", 100, -1.f, 4.f);
//  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_x"), "ttbar_x_negative", "", 100, -100.f, 100.f);
//  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_y"), "ttbar_y_negative", "", 100, -100.f, 100.f);
//  hist_negative.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_z"), "ttbar_z_negative", "", 100, -100.f, 100.f);


  //spin correlation
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHel"), "spin_cHel_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cLab"), "spin_cLab_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "crr"), "spin_crr_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cnn"), "spin_cnn_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "ckk"), "spin_ckk_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cSca"), "spin_cSca_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cTra"), "spin_cTra_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi0"), "spin_phi0_negative", "", 100, 0.f, 4.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "phi1"), "spin_phi1_negative", "", 100, 0.f, 7.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTP"), "spin_cpTP_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cpTTT"), "spin_cpTTT_negative", "", 100, -1.f, 1.f);
  hist_negative.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHan"), "spin_cHan_negative", "", 100, -1.f, 1.f);

//  Histogram hist_cut_juan_paper;
//  //hist_cut_juan_paper.set_weighter([&weight = metadata.get<float>("weight")] () { return weight[0]; });
//  hist_cut_juan_paper.set_weighter([&weight = metadata.get<float>("weight"), &weight_new =lhe_event.get<float>("weight_juan_paper_float")] () { return (weight[0] * weight_new[0]); });
//  hist_cut_juan_paper.make_histogram<TH1F>(filler_first_of(gen_ttbar, "ttbar_mass"), "ttbar_mass_cut", "", 120, 300.f, 1500.f);
//  hist_cut_juan_paper.make_histogram<TH1F>(filler_first_of(gen_tt_ll_bb, "cHel"), "spin_cHel", "", 100, -1.f, 1.f);
//

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
  auto f_analyze = [&metadata, &lhe_particle, &lhe_event, &gen_particle, &gen_ttbar, &gen_tt_ll_bb, &hist_no_cut, &hist_cut, &hist_positive, &hist_negative, &tree_gen] (long long entry) {
    // first we start by populating the collections
    // this is essentially equivalent of the tree->GetEntry(entry)
    // with the (compulsory) freedom of timing the call separately for each group
//    std::cout << "populate metadata" << std::endl;
    metadata.populate(entry);
//    std::cout << "populate gen_particle" << std::endl;
    gen_particle.populate(entry);
//    std::cout << "populate lhe_particle" << std::endl;
    lhe_particle.populate(entry);

    // since the collections serve as input to the aggregates, they need to be populated first
//    std::cout << "populate gen_ttbar" << std::endl;
    gen_ttbar.populate(entry);
//    std::cout << "populate gen_tt_ll_bb" << std::endl;
    gen_tt_ll_bb.populate(entry);
//    std::cout << "populate lhe_event" << std::endl;
    lhe_event.populate(entry);

//    std::cout << "populates finished" << std::endl;

    //printing
    //std::cout << "normal weight: ";
   // metadata.iterate(printer_normal, metadata.ref_to_indices(), "weight");
   // std::cout << "float weight: ";
   // lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "weight_float");
   // std::cout << "double weight: ";
   // lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "weight_double");
    std::cout << "cpTP LHE: ";
    lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "cpTP_LHE");
    std::cout << "t cos theta LHE: ";
    lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "t_CosTheta_LHE");
    std::cout << "cos theta - cpTP LHE: ";
    lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "diff_cos_theta_cpTP_LHE");
    std::cout << "cpTP gen: ";
    gen_tt_ll_bb.iterate(printer_normal, gen_tt_ll_bb.ref_to_indices(), "cpTP");
    std::cout << "t cos theta gen: ";
    gen_ttbar.iterate(printer_normal, gen_ttbar.ref_to_indices(), "t_CosTheta");

    //std::cout << "normal * double: " << metadata.get<float>("weight") * lhe_event.get<float>("weight_float") << std::endl;

//    std::cout << "ttbar mass: ";
//    lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "ttbar_mass");
    //std::cout << "antitop mass: ";
    //lhe_event.iterate(printer_normal, lhe_event.ref_to_indices(), "tbar_mass");
    //lhe_event.iterate(printer_reweighted, lhe_event.ref_to_indices(), "weight_precision");
    
    // we make an oversimplification here, considering only the events where gen_tt_ll_bb contain an element
    // this is because in the above, we have grouped the gen_ttbar and gen_tt_ll_bb histograms together
    // despite the fact that the requirements of gen_tt_ll_bb is strictly tighter than gen_ttbar
    // without this restriction, when filling the hist_no_cut histograms we get spurious entries in the gen_tt_ll_bb histograms
    // when this aggregate is empty e.g. when we have taus in the event
    if (!gen_tt_ll_bb.n_elements())
      return;

    if (!lhe_event.n_elements())
      {
        std::cout << "No elements in lhe_event" << std::endl;
      	return;
      }

    // printing stuff
    // metadata.iterate(logger(std::cout), -1, -1, "lumi", "event");
    // gen_tt_ll_bb.iterate(printer, -1, -1, "lbbar_mass");
    //gen_particle.iterate(printer, -1, -1, "pt");
    //gen_particle.iterate(printer, -1, -1, "pt");

    // fill the no (acceptance) cut histograms
    hist_no_cut.fill();

    // as promised above we would like some acceptance cuts 
    // which we impose on the charged leptons and bottom quarks
    // we could have used the daughters that we transferred to the gen_tt_ll_bb aggregate
    // but to highlight some additional features we will instead do it through the gen_particle collections instead
    // begin by selecting the daughters among all the gen particles using a generic filter method
    // which needs a function that evaluates to true or false based on a list of attributes
//    auto lepton_bottom_passing_pt_eta_cut = gen_particle.filter([] (int tag, float pt, float eta) {
//        // not lepton or bottom, reject
//        if (tag != 9 and tag != 4 and tag != 3 and tag != 8)
//          return false;
//
//        if (pt > 20.f and std::abs(eta) < 2.4f)
//          return true;
//        else
//          return false;
//      }, "dileptonic_ttbar", "pt", "eta");
//
//    // recall that filter methods return a list of indices
//    // to overwrite the indices list of the group, we use the update_indices method
//    gen_particle.update_indices(lepton_bottom_passing_pt_eta_cut);
//
//    // if all four objects pass the cut, then gen_particle will have 4 elements left
//    // fill also our tree at this point
//    if (gen_particle.n_elements() == 4) {
//      hist_cut.fill();
//      tree_gen.fill();
//    }

    auto passll = gen_tt_ll_bb.filter_greater("lepton_pt", 20.f, 
                                              gen_tt_ll_bb.filter_in("lepton_eta", -2.4f, 2.4f,
                                                                     gen_tt_ll_bb.filter_greater("antilepton_pt", 20.f,
                                                                                                 gen_tt_ll_bb.filter_in("antilepton_eta", -2.4f, 2.4f))));
    auto passllbb = gen_tt_ll_bb.filter_greater("bottom_pt", 20.f, 
                                               gen_tt_ll_bb.filter_in("bottom_eta", -2.4f, 2.4f,
                                                                      gen_tt_ll_bb.filter_greater("antibottom_pt", 20.f,
                                                                                                  gen_tt_ll_bb.filter_in("antibottom_eta", -2.4f, 2.4f, 
                                                                                                                         passll)))).size();


    if (passllbb) {
      hist_cut.fill();
//      hist_cut_juan_paper.fill();
      tree_gen.fill();
    }
    
    auto positive_event = lhe_event.filter_greater("weight_float", 0.f);
    if (positive_event){
       hist_positive.fill();
    }
    else{
       hist_negative.fill(); 
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
 
  std::string suffix = "N_and_P_real";
 
  std::string filename_cut = create_filename("hist_ttbarlo_reweighting", higgs_type, mass, width, calc_weight_variant, res_int, (std::string) "cut_after_reordering_" + suffix);
  std::string filename_nocut = create_filename("hist_ttbarlo_reweighting", higgs_type, mass, width, calc_weight_variant, res_int, (std::string) "no_cut_after_reordering_" + suffix);
  std::string filename_positive = create_filename("hist_ttbarlo_reweighting", higgs_type, mass, width, calc_weight_variant, res_int, (std::string) "no_cut_positive_after_reordering_" + suffix);
  std::string filename_negative = create_filename("hist_ttbarlo_reweighting", higgs_type, mass, width, calc_weight_variant, res_int, "no_cut_negative_after_reordering_" + suffix);

  std::cout << "Saving as: " << filename_cut << std::endl;

  hist_no_cut.save_as(filename_nocut);
  hist_cut.save_as(filename_cut);
  hist_positive.save_as(filename_positive);
  hist_negative.save_as(filename_negative);
  tree_gen.save();

  return 0;
}

