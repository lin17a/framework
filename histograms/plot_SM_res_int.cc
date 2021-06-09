#include<algorithm>
#include<vector>
#include<sys/stat.h>
#include<limits>

{
  struct histogram {
    TH1F *hist;
    std::string name;
    Color_t line_color;
  };
  
  struct histogram_raw {
    TH1F *hist;
    bool status;
  };
  
  auto SetUpHistograms = [](std::vector<histogram *> *hists, TLegend *legend, int line_width = 3, int rebin = 2) {
    for (histogram *hist : *hists) {
      hist->hist->SetLineColor(hist->line_color);
      hist->hist->SetLineWidth(line_width);
      hist->hist->Rebin(rebin);
      int nbins = hist->hist->GetXaxis()->GetNbins();
      hist->hist->Scale(1. / hist->hist->Integral());
      for (int i=0; i<nbins ;i++){
        Float_t v = hist->hist->GetBinContent(i);
        //if (abs(v) < 10){
        //  v = 0;
        //}
        hist->hist->SetBinContent(i, abs(v));
      }
    } 

    std::vector<double_t> maxs;
    std::transform(hists->begin(), hists->end(), std::back_inserter(maxs),
                   [](auto h) -> double_t { std::cout <<  std::to_string(h->hist->GetMaximum()) << "\n"; return h->hist->GetMaximum(); });
    double_t max_value = *(std::max_element(maxs.begin(), maxs.end()));
    std::vector<double_t> mins;
    std::transform(hists->begin(), hists->end(), std::back_inserter(mins),
                   [](auto h) -> double_t { std::cout <<  std::to_string(h->hist->GetMinimum()) << "\n"; return h->hist->GetMinimum(); });
    double_t min_value = *(std::min_element(mins.begin(), mins.end()));


   
    if (min_value <= 0) {
       std::cout << "min too small" << std::endl;
       min_value = 0.0000001;
    }
 
    for (histogram *hist : *hists) { 
      hist->hist->SetMaximum(max_value + max_value );
      hist->hist->SetMinimum(min_value);
      //hist->GetYaxis()->SetRangeUser(min_value, max_value + max_value / 8);
      legend->AddEntry(hist->hist, hist->name.c_str(), "l");
      std::cout << std::to_string(hist->hist->GetMaximum()) << std::endl;
    }
  };

  auto GetHistogram = [](TFile *file, std::string hist_name){
    bool status = true;
    TH1F *h1;
    std::set<std::string> set_hist_names;
    TIter iter(file->GetListOfKeys());
    TKey *key;
    //std::cout << "Keys in h3:" << std::endl;
    while ((key=(TKey*)iter.Next()) ) {  
       set_hist_names.insert(key->GetName());
       //std::cout << key->GetName() << std::endl;
    }
    if (set_hist_names.find(hist_name) != set_hist_names.end()){
      h1 = (TH1F *) file->Get(hist_name.c_str());
      std::cout << "Habe histogram " << hist_name << " auch in h3 gefunden!" << std::endl;
    }
    else { 
      size_t pos_1 = hist_name.find("_no_cut");
      std::string attr = hist_name.substr(0, pos_1);
      std::string hist_name_positive = attr + "_positive";
      std::string hist_name_negative = attr + "_negative";
      if (set_hist_names.find(attr) != set_hist_names.end()){
        h1 = (TH1F *) file->Get(attr.c_str());
        std::cout << "attr: " << attr << std::endl; 
      }
      else if (set_hist_names.find(hist_name_positive) != set_hist_names.end()){
        h1 = (TH1F *) file->Get(hist_name_positive.c_str());
      }
      else if (set_hist_names.find(hist_name_negative) != set_hist_names.end()){
        h1 = (TH1F *) file->Get(hist_name_negative.c_str());
      }
      else{
        std::cout << "Habe histogram " << hist_name << " nicht in h3 gefunden!" << std::endl;
        status = false;
      }
    }
    histogram_raw hist_struct = {h1, status};
    return hist_struct; 
  };

  auto MakeTitle = []( auto variable, std::string res_int, std::string higgs_type, std::string mass, std::string width){
          std::string higgs_name;
          std::string subscript; 
          if (higgs_type == "pseudo_scalar"){
            higgs_name = "pseudo scalar";
            subscript = "A";
          }
          else {
            higgs_name = "scalar";
            subscript = "H";
          }

          //std::string variable;
          //if (hist_name == "ttbar_mass_no_cut"){
          //  variable = "m_{t#bar{t}}";
          //}
          //std::string string_variable = (std::string) variable;

          std::string hist_title = variable + "\\text{ for\\ " + higgs_name + "\\ heavyhiggs,\\ " + res_int + ", } m_{" + subscript + "}=\\SI{" + mass + "}{\\giga\\electronvolt}, " + "\\Gamma_{" + subscript + "}=\\SI{" + width + "}{\\giga\\electronvolt}";
          return hist_title;
   };


  // map to get the correct unit for every attribute
  std::unordered_map<std::string, std::string> units = {
    {"pt","GeV"},
    {"mass","GeV"},
    {"et","GeV"},
    {"et2","GeV"},
    {"mass_2","GeV"},
    {"m2","GeV"},
    {"mt","GeV"},
    {"mt2","GeV"},
    {"mag","GeV"},
    {"mag2","GeV"},
    {"mt","GeV"},
    {"mt2","GeV"},
    {"px","GeV"},
    {"py","GeV"},
    {"pz","GeV"}
  };

  // map to get the correct identifier attribute
  std::unordered_map<std::string, std::string> variables = {
    {"ttbar_mass_no_cut", "m_{t\\bar{t}}"},
    {"antibottom_eta_no_cut", "\\eta_{\\bar{b}}"},
    {"antilepton_pt_no_cut", "p^{\\perp}_{\\bar{l}}"},
    {"lbarb_mass_no_cut", "m_{\\bar{l}b}"},
    {"llbar_mass_no_cut", "m_{l\\bar{l}}"},
    {"spin_cTra_no_cut", "spin_cTra_no_cut"},
    {"spin_crr_no_cut", "\\, c_{rr}"},
    {"bbar_mass_no_cut", "m_{b\\bar{b}}"},
    {"lbarbbar_mass_no_cut", "m_{\\bar{l}\\bar{b}}"},
    {"spin_cHan_no_cut", "c_{\\text{han}}"},
    {"spin_ckk_no_cut", "\\, c_{kk}"},
    {"spin_phi0_no_cut", "\\varphi_0"},
    {"ttbar_pt_no_cut", "p^{\\perp}_{t\\bar{t}}"},
    {"antibottom_eta_no_cut", "\\eta_{\\bar{b}}"},
    {"bottom_eta_no_cut", "\\eta_{b}"},
    {"lbbar_mass_no_cut", "m_{l\\bar{b}}"},
    {"spin_cHel_no_cut", "c_{\\text{hel}}"},
    {"spin_cnn_no_cut", "\\, c_{nn}"},
    {"spin_phi1_no_cut", "\\varphi_1"},
    {"ttbar_pt_transform_no_cut", "p^{\\perp}_{t\\bar{t}}"},
    {"antibottom_pt_no_cut", "p^{\\perp}_{\\bar{b}}"},
    {"bottom_pt_no_cut", "p^{\\perp}_{b}"},
    {"lepton_eta_no_cut", "\\eta_{l}"},
    {"spin_cLab_no_cut", "c_{\\text{lab}}"},
    {"spin_cpTP_no_cut", "\\, c_{t, p}"},
    {"t_CosTheta_no_cut", "cos\\Theta_{t}"},
    {"antilepton_eta_no_cut", "\\eta_{\\bar{l}}"},
    {"lb_mass_no_cut", "\\, m_{lb}"},
    {"lepton_pt_no_cut", "p^{\\perp}_{l}"},
    {"spin_cSca_no_cut", "spin_cSca_no_cut"},
    {"spin_cpTTT_no_cut", "c_{t, t\\overline{t}}"},
    {"ttbar_mass_2_no_cut", "m_{t\\bar{t}}"},
    {"spin_b1k_no_cut", "\\, b^1_k"},
    {"spin_b2k_no_cut", "\\, b^2_k"},
    {"spin_b1r_no_cut", "\\, b^1_r"},
    {"spin_b2r_no_cut", "\\, b^2_r"},
    {"spin_b1n_no_cut", "\\, b^1_n"},
    {"spin_b1n_no_cut", "\\, b^2_n"}
  };

  auto DrawHistograms = [ SetUpHistograms, GetHistogram, MakeTitle, units, variables ](std::string filename_generated_res, std::string filename_generated_int, std::string filename_generated_SM, std::string filename_reweighted, std::string folder, std::string filename_histogram, std::string res_int, std::string higgs_type, std::string mass, std::string width){
   
    TFile *f_res = TFile::Open( filename_generated_res.c_str() );
    TFile *f_int = TFile::Open( filename_generated_int.c_str() );
    TFile *f_SM = TFile::Open( filename_generated_SM.c_str() );
    TFile *f_rew = TFile::Open( filename_reweighted.c_str() );

    if (!f_res->IsOpen()  or !f_int->IsOpen() or !f_SM->IsOpen() or !f_rew->IsOpen()){
      std::cout << "File does not exist." << std::endl;
      return 0;
    }


    int n_histograms = f_res->GetListOfKeys()->GetSize();
    std::cout << "number of histograms: " << n_histograms << std::endl;
    TIter iter1(f_res->GetListOfKeys());
    TKey *key1;

    int i;

    // iterate over all histograms in the first file
    while ((key1=(TKey*)iter1.Next()) ) {
      
       std::string hist_name = key1->GetName();
       TH1F *h_res = (TH1F *) f_res->Get( hist_name.c_str());   
       std::string title = hist_name;     

      // try to find histograms with the same name 
      try {

          std::cout << hist_name << std::endl;
          
          // get histograms with the same attribute from the other files 
          histogram_raw hist_struct_int = GetHistogram(f_int, hist_name);   
          histogram_raw hist_struct_SM = GetHistogram(f_SM, hist_name);   
          histogram_raw hist_struct_rew = GetHistogram(f_rew, hist_name);   
          TH1F *h_int;
          TH1F *h_SM;
          TH1F *h_rew;
          if (hist_struct_int.status and hist_struct_SM.status and hist_struct_rew.status ) {
            h_int = hist_struct_int.hist;
            h_SM = hist_struct_SM.hist;
            h_rew = hist_struct_rew.hist;
            std::cout << "Habe alle anderen hists gesetzt" << std::endl;
          }
          else {
            continue;
          }

          // defining histogram structs with a histogram, color and name to hand it all together to the SetUp function
         // histogram hist1;
         // histogram hist3;
         // if (higgs_type == "pseudo_scalar"){
         //   hist1 = {h1, "\\text{generated } A", kRed + 1};
         //   hist3 = {h3, "\\text{reweighted } A", kBlue+2};
         // }
         // else {
         //   hist1 = {h1, "\\text{generated } H", kGreen + 2};
         //   hist3 = {h3, "\\text{reweighted } H", kBlue+2};
         // }

          h_res->Scale(1.389 / h_res->Integral());
          h_int->Scale(-0.147 / h_int->Integral());
          h_SM->Scale(420.3 / h_SM->Integral());

          TH1F *h_gen = (TH1F*) h_res->Clone();
          h_gen->Add(h_int);
          //h_gen->Add(h_SM);
         
          histogram hist_gen = {h_gen, "generated", kRed};
          histogram hist_rew = {h_rew, "reweighted", kBlue};

          //histogram hist1 = {h1, legend_1, kGreen-3};
          //histogram hist2 = {h2, "int pseudo", kGreen+3};
          
          //double_t chi2 = h1->TH1::Chi2Test(h3, "WW");
          //std::cout << hist_name << ", chi^2 = " << chi2 << std::endl;

         // add a new histogram containing the sum of ttbar and the pseudo scalar higgs histograms
          //TH1F *h_sum_res_int_pseudo = (TH1F*) h1->Clone();
        //  h_sum_res_int_pseudo->Add(h2);
          //h_sum_ttbarlo_res_int_pseudo->Add(h4);
          //h_sum_res_int_pseudo->SetLineColor(kAzure+3);

         // histogram hist_sum_pseudo = {h_sum_res_int_pseudo, "sum pseudo", kRed+2};

          // initialize legend, numbers are the coordinates of the bottom left corner und upper right corner
          // fits exactly under the info panel
          //auto legend = new TLegend(0.78,0.69,0.98,0.77);
          //auto legend = new TLegend(0.78,0.79,0.98,0.87);
          auto legend = new TLegend(0.78,0.79,1.05,0.91);
          legend->SetTextSize(0.02);        
  
          std::vector<histogram *> hists{&hist_gen, &hist_rew};

          // adjust colors, scales, ... 
          SetUpHistograms(&hists, legend, 2, 1);
          std::cout << "after setup histograms" << std::endl;       
 
          auto variable2 = variables.find(hist_name);

          std::pair<std::string,std::string> variable;
          if (variable2 == variables.end()){
            variable = std::make_pair(hist_name, hist_name);
          } else {
            variable = *variable2;
          }

          // find and set the right title for the x axis
          size_t pos_1 = hist_name.find('_');
          size_t pos_2 = hist_name.find("_no_cut", pos_1 + 1);
          std::string attr = hist_name.substr(pos_1 + 1, pos_2 - pos_1 - 1);
          std::string x_title; 
          auto unit = units.find(attr);
   
          //if (hist_name == "ttbar_mass_no_cut"){
          //  attr = "m_{t#bar{t}}";
          //}
          if (unit != units.end()) {
            x_title = variable.second + "\\ [\\text{" + unit->second + "}]"; 
            }
          else {
            x_title = variable.second;
          }
 
          //h1->GetXaxis()->SetTitle(l);
          h_gen->GetXaxis()->SetTitle(x_title.c_str());
          h_gen->GetYaxis()->SetTitle("normalized number of events");
          h_gen->GetYaxis()->SetTitleOffset(4);
          //h1->GetYaxis()->SetNdivisions();
          //h1->GetYaxis()->CenterTitle(true);
          h_gen->SetTitleSize(.05, "XY");
          h_gen->SetTitleSize(.08, "t");

          std::string hist_title = MakeTitle(variable.second, res_int, higgs_type, mass, width);

          h_gen->SetNameTitle(hist_name.c_str(), hist_title.c_str());
          h_gen->SetLabelSize(.05, "XY");
          //gStyle->SetTitleFontSize(20);
          
          TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
          //TCanvas can("canvas",
          if (attr == "mass") 
             can->SetLogy();
          if (hist_name == "ttbar_mass_no_cut") 
             can->SetLogy();
          can->cd();
             
          
          // draw all in the same histogram
          h_gen->Draw("pe");
//          h2->Draw("pesame");
          h_rew->Draw("pesame");
          //gStyle->SetTitleFontSize(0.7);

          can->Update();
          can->Clear();
         
          gStyle->SetLabelSize(.04, "XY");
          //gStyle->SetTitleFontSize(.11);
          gStyle->SetTitleSize(0.05, "t");
          gStyle->SetTitleSize(.04, "XY");
          gStyle->SetOptStat(0);
 
          // add a sub pad with the ratio 
          auto rp1 = new TRatioPlot(h_gen, h_rew);
          rp1->Draw();
          rp1->GetLowerRefYaxis()->SetTitle("\\frac{\\text{generated}}{\\text{reweighted}}");
          rp1->GetLowerRefYaxis()->SetRangeUser(0,2);
          rp1->GetLowerRefYaxis()->SetTitleOffset(2.5);
          rp1->GetLowerRefXaxis()->SetTitleOffset(2.5);
          rp1->GetLowerRefYaxis()->CenterTitle(true);
          //rp1->Draw();
          //rp1->GetYaxis->SetNdivisions(503);
          //rp1->GetLowerRefXaxis()->SetTitleSize(20);
          rp1->GetXaxis()->SetLabelSize(0.023);
          //rp1->SetMinimum(-10);
          //rp1->SetMaximum(10); 
          can->Update();


          legend->Draw();
 
          // save with correct name
          std::string filename_single_hist = folder + "/" + std::string(hist_name) + ".pdf";
          if (i == 0){
            filename_histogram = filename_histogram + "[";          
          }
          else if (i == n_histograms){
            filename_histogram = filename_histogram + "]";
          }
          std::string title = std::string(hist_name);
          can->Print(filename_histogram.c_str(), title.c_str());
          can->SaveAs( filename_single_hist.c_str() );
          delete can;
          delete rp1; 
       }
       catch (const std::exception&) { 
           std:cout << "exception gecatcht\n"; 
       }
       i++;
    }
    return 0;
};





  //std::vector<std::string> mass_width_combis = { "m400_w20", "m1000_w25", "m600_w30", "m800_w20" };
  //std::vector<std::string> masses = { "400",  "1000", "600", "800" };
  //std::vector<std::string> widths = { "20", "25", "30", "20" };
  std::vector<std::string> masses = { "600" };
  std::vector<std::string> widths = { "30" };
  std::vector<std::string> zeros = {"00", "000", "00", "000"};
  std::vector<std::string> res_int =  {"res"}; //, "int"}; //
  std::vector<std::string> RES_INT = {"RES"}; //, "INT"}; //
  std::vector<std::string> r_i = {"r"}; //, "i"}; //
  std::vector<std::string> resonance_interference = {"resonance"}; //, "interference"}; //
  std::vector<std::string> PSEUDO_SCALAR = { "PSEUDO"}; //,  {"SCALAR",
  std::vector<std::string> pseudo_scalar = { "pseudo_scalar"}; // "scalar",
  std::vector<std::string> p_s = {"p"}; // "s", 

  std::vector<std::string> positive_negative; 


  for (int i = 0; i < masses.size(); i++){
    for (int j = 0; j < pseudo_scalar.size(); j++){
      for (int l = 0; l < res_int.size(); l++){

        if (res_int[l] == "int"){
          positive_negative = {"_positive", "_negative"}; 
        }
        else {
          positive_negative = {""};
        }

        for (int k = 0; k < positive_negative.size(); k++){
 

          std::string filename_generated_res = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_madspin_m" + masses[i] + "_w" + widths[i] + "_" + zeros[i]  + "_" + "RES" +"_" + PSEUDO_SCALAR[j] + "_generated/hist_madspin_m" + masses[i] + "_w" + widths[i] + "_" + zeros[i] +  "_" + "RES" + "_" + PSEUDO_SCALAR[j] + "_no_cut" + positive_negative[k] + "_madspin.root";
          std::string filename_generated_int = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_madspin_m" + masses[i] + "_w" + widths[i] + "_" + zeros[i]  + "_" + "INT" +"_" + PSEUDO_SCALAR[j] + "_generated/hist_madspin_m" + masses[i] + "_w" + widths[i] + "_" + zeros[i] +  "_" + "INT" + "_" + PSEUDO_SCALAR[j] + "_no_cut" + positive_negative[k] + "_madspin.root";
          std::string filename_generated_SM = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[l] + "_m" + masses[i] + "_w" + widths[i] + "_reweighted_juan_paper_fixed_width/hist_ttbarlo_no_reweighting_" + pseudo_scalar[j] + "_m" + masses[i] + "_w" + widths[i] + "_juan_paper_" + resonance_interference[l] + "_no_cut" + positive_negative[k] + "_after_reordering_with_madspin_more_spin.root";

          //std::cout << "filename_generated: " << filename_generated << std::endl;

          std::string filename_reweighted = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + "ri" + "_m" + masses[i] + "_w" + widths[i] + "_reweighted_juan_code_fixed_width/hist_ttbarlo_reweighting_" + pseudo_scalar[j] + "_m" + masses[i] + "_w" + widths[i] + "_juan_code_" + "res_int" + "_no_cut" + positive_negative[k] + "_after_reordering_with_madspin_juan_code_test.root";
          //std::string filename_reweighted = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[l] + "_m" + masses[i] + "_w" + widths[i] + "_reweighted_juan_paper_fixed_width/hist_ttbarlo_no_reweighting_" + pseudo_scalar[j] + "_m" + masses[i] + "_w" + widths[i] + "_juan_paper_" + resonance_interference[l] + "_no_cut" + positive_negative[k] + "_after_reordering_with_madspin_more_spin.root";
          //std::string filename_reweighted = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[l] + "_m" + masses[i] + "_w" + widths[i] + "_reweighted_juan_paper_fixed_width/hist_ttbarlo_reweighting_" + pseudo_scalar[j] + "_m" + masses[i] + "_w" + widths[i] + "_juan_paper_" + resonance_interference[l] + "_no_cut" + positive_negative[k] + "_after_reordering_with_madspin_juan_code_test.root";

          std::cout << "filename_reweighted: " << filename_reweighted << std::endl;
 
          std::string folder =  "m" + masses[i] + "_w" + widths[i] + "_" + pseudo_scalar[j] + "_" + res_int[l] + positive_negative[k] + "_generated_vs_reweighted_res_int_SM";
          //std::string folder =  "m" + masses[i] + "_w" + widths[i] + "_" + pseudo_scalar[j] + "_" + res_int[l] + positive_negative[k] + "_madspin_juan_paper_generated_vs_reweighted";
          //std::string folder = "test"; 
          if (mkdir(folder.c_str(), 0777) != 0){
            std::cout << "Couldn't create directory." << std::endl;
          }
          
          std::string filename_histogram = folder + "/" +  std::string("all") + std::string(".pdf");
          std::cout << "filename_histogram: " << filename_histogram << std::endl;
          
          DrawHistograms(filename_generated_res, filename_generated_int, filename_generated_SM, filename_reweighted, folder, filename_histogram, resonance_interference[l], pseudo_scalar[j], masses[i], widths[i]);
        }
      }
    }
  }

  //std::string suffix("_m400_w20_pseudo_generated_res_int_vs_reweighting_after_reordering_juan_code.pdf");
  // std::string folder("m_400_w20_pseudo_scalar_res_generated_vs_reweighted");
  // std::string filename = folder + "/" +  std::string("all") + std::string(".pdf");
  // if (mkdir(folder.c_str(), 0777) != 0){
  //   std::cout << "Couldn't create directory." << std::endl;
  // }
  // add all root files
  
 // .!ls;
  //.!mkdir suffix
  //.!mv *suffix suffix
  //.!cd suffix
  //.!pdfunite *suffix
  //.!mkdir suffix
  //.!mv *suffix suffix
  //.!cd suffix
  //.!pdfunite *suffix '"all_" + suffix' 
}
//.!ls

