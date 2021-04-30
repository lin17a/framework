#include<algorithm>
#include<vector>
#include <sys/stat.h>

{
  struct histogram {
    TH1F *hist;
    std::string name;
    Color_t line_color;
  };
  
  auto SetUpHistograms = [](std::vector<histogram *> *hists, TLegend *legend, int line_width = 3, int rebin = 4) {
    for (histogram *hist : *hists) {
      hist->hist->SetLineColor(hist->line_color);
      hist->hist->SetLineWidth(line_width);
      hist->hist->Rebin(rebin);
      // set all bins under 400 GeV to zero
      //for (int bin = 0; bin < 10; bin++){
      //  hist->hist->SetBinContent(bin, 0);
      //}
      hist->hist->Scale(1. / hist->hist->Integral());
    } 

    std::vector<double_t> maxs;
    std::transform(hists->begin(), hists->end(), std::back_inserter(maxs),
                   [](auto h) -> double_t { std::cout <<  std::to_string(h->hist->GetMaximum()) << "\n"; return h->hist->GetMaximum(); });
    double_t max_value = *(std::max_element(maxs.begin(), maxs.end()));
    std::vector<double_t> mins;
    std::transform(hists->begin(), hists->end(), std::back_inserter(mins),
                   [](auto h) -> double_t { std::cout <<  std::to_string(h->hist->GetMinimum()) << "\n"; return h->hist->GetMinimum(); });
    double_t min_value = *(std::min_element(mins.begin(), mins.end()));

    std::cout << "max: " << max_value << std::endl;
   
    if (min_value <= 0) {
       std::cout << "min too small" << std::endl;
       min_value = 0.0000001;
    }
 
    for (histogram *hist : *hists) { 
      hist->hist->SetMaximum(max_value + max_value/ 8);
      hist->hist->SetMinimum(min_value);
      //hist->GetYaxis()->SetRangeUser(min_value, max_value + max_value / 8);
      legend->AddEntry(hist->hist, hist->name.c_str(), "l");
      std::cout << std::to_string(hist->hist->GetMaximum()) << std::endl;
    }
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


  std::vector<std::string> log_scale_hists = {"ttbar_mass_cut"};

  auto DrawHistograms = [ SetUpHistograms, units, log_scale_hists ](std::vector<std::string> filenames_higgs_reweighted, std::string filename_reweighted_ttbarlo, std::string folder, std::string filename_histogram, std::string res_int, std::string higgs_type, std::string mass_width){
   
    std::vector<TFile *> higgs_files;
    std::vector<TH1F *> higgs_hs;


    higgs_files.clear();
    for (int i = 0; i < filenames_higgs_reweighted.size(); i++){
        higgs_files.push_back( TFile::Open(filenames_higgs_reweighted[i].c_str()) );
    }

    TFile *f1 = TFile::Open( filename_reweighted_ttbarlo.c_str() );
    //TFile *f2 = TFile::Open( "hist_ttbarlo_reweighting_pseudo_scalar_m1000.000000_w25.000000_juan_paper_resonance_cut_after_reordering.root" );
    //TFile *f2 = TFile::Open( filename_reweighted_higgs_1.c_str() );
    //TFile *f3 = TFile::Open( filename_reweighted_higgs_2.c_str() );
    //TFile *f4 = TFile::Open( filename_reweighted_higgs_3.c_str() );
    //TFile *f5 = TFile::Open( filename_reweighted_higgs_4.c_str() );

    if (!f1->IsOpen()){ // or !f3->IsOpen()){
      std::cout << "File does not exist." << std::endl;
      return 0;
    }

    int n_histograms = f1->GetListOfKeys()->GetSize();
    std::cout << "number of histograms: " << n_histograms << std::endl;
    TIter iter1(f1->GetListOfKeys());
    TKey *key1;

    int i;

    // iterate over all histograms in the first file
    while ((key1=(TKey*)iter1.Next()) ) {
      
       std::string hist_name = key1->GetName();
       // get the first histogram from the first file
       TH1F *h1 = (TH1F *) f1->Get(hist_name.c_str());
       std::string title = hist_name;
       std::cout << "hist name ttbar: " << hist_name << std::endl;     

      // try to find histograms with the same name 
      try {
          
          // get histograms with the same attribute from the other files 
        higgs_hs.clear();
        for (int i = 0; i < higgs_files.size(); i++){
          higgs_hs.push_back( (TH1F *) higgs_files[i]->Get(hist_name.c_str()) );
          std::cout << "hist name heavyhiggs: " << higgs_files[i]->Get(hist_name.c_str())->GetName() << std::endl;
        }
        std::cout << "hier sollten jetzt " << higgs_files.size() << " viele Titel stehen!" << std::endl;
      //    TH1F *h2 = (TH1F *) f2->Get( hist_name.c_str());   
          //TH1F *h2 = (TH1F *) f2->Get( hist_name.c_str());   
          //TH1F *h3 = (TH1F *) f3->Get( hist_name.c_str());   
          //TH1F *h4 = (TH1F *) f4->Get( hist_name.c_str());   
          //TH1F *h5 = (TH1F *) f5->Get( hist_name.c_str());   

          // defining histogram structs with a histogram, color and name to hand it all together to the SetUp function
          histogram hist1 = {h1, "reweighted ttbarlo", kGreen-3};
        //  histogram hist2 = {h2, "int pseudo", kGreen+3};
        //  histogram hist3 = {h3, "reweighted heavyhiggs", kRed}; //kBlue+2
          
          // add a new histogram containing the sum of ttbar and the pseudo scalar higgs histograms
          TH1F *h_sum_heavyhiggs = (TH1F*) higgs_hs[0]->Clone();
       
          for (int i = 1; i < higgs_hs.size(); i++){
            h_sum_heavyhiggs->Add(higgs_hs[i]);
          }
          //h_sum_heavyhiggs->Add(h3);
          //h_sum_heavyhiggs->Add(h4);
          //h_sum_heavyhiggs->Add(h5);
          //h_sum_res_int_pseudo->SetLineColor(kAzure+3);

          histogram hist_sum_heavyhiggs = {h_sum_heavyhiggs, "sum heavyhiggs", kRed+2};

          // initialize legend, numbers are the coordinates of the bottom left corner und upper right corner
          // fits exactly under the info panel
          auto legend = new TLegend(0.78,0.69,0.98,0.77);
          
          std::vector<histogram *> hists{&hist1, &hist_sum_heavyhiggs};

          // adjust colors, scales, ... 
          SetUpHistograms(&hists, legend, 2, 1);
          
          // find and set the right title for the x axis
          size_t pos_1 = hist_name.find('_');
          size_t pos_2 = hist_name.find("_no_cut", pos_1 + 1);
          std::string attr = hist_name.substr(pos_1 + 1, pos_2 - pos_1 - 1);
          std::string x_title; 
          auto unit = units.find(attr);
          if (unit != units.end()) {
            x_title = attr + " in " + unit->second; 
            }
          else {
            x_title = attr;
          }
          
          h1->GetXaxis()->SetTitle(x_title.c_str());
          h1->GetYaxis()->SetTitle("number of events");
 
          std::string hist_title = higgs_type + " " + res_int + " " + mass_width;
          h1->SetNameTitle(hist_name.c_str(), hist_title.c_str());
          
          TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
          //if the hist name is in the log scale list, make a log scale
          if (find(log_scale_hists.begin(), log_scale_hists.end(), hist_name) != log_scale_hists.end()){ 
            std::cout << "Making log scale for " << hist_name << std::endl;
            can->SetLogy();
            can->cd();
          }

          // draw all in the same histogram
          h1->Draw("pe");
//          h2->Draw("pesame");
          h_sum_heavyhiggs->Draw("pesame");

          can->Clear();
          
          // add a sub pad with the ratio 
          auto rp1 = new TRatioPlot(h1, h_sum_heavyhiggs);
          rp1->Draw();
          rp1->GetLowerRefYaxis()->SetLabelSize(0.02);
          rp1->GetLowerRefYaxis()->SetTitleSize(0.02);
          rp1->GetLowerRefYaxis()->SetTitle("reweighted ttbar / sum heavyhiggs");
          rp1->GetLowerRefYaxis()->SetRangeUser(0,2);
//          rp1->SetMinimum(-10);
//          rp1->SetMaximum(10); 
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



// leaving out m600_w90 and m1000_w150

  std::vector<std::string> mass_width_combis = { "m400_w20", "m1000_w25", "m600_w30", "m800_w20", "m1000_w100", "m1000_w50", "m400_w10", "m400_w40", "m400_w60", "m500_w25", "m500_w50", "m500_w75", "m600_w15", "m600_w60", "m800_w120", "m800_w40", "m800_w80"};
  std::vector<std::string> zeros = { "00", "000", "00", "000"};
  std::vector<std::string> res_int = {"res"}; //, "int"};
  std::vector<std::string> RES_INT = {"RES"}; //, "INT"};
  std::vector<std::string> r_i = {"r"}; //, "i"};
  std::vector<std::string> resonance_interference = {"resonance"}; //, "interference"};
  std::vector<std::string> PSEUDO_SCALAR = {"PSEUDO", "SCALAR"};
  std::vector<std::string> pseudo_scalar = {"pseudo_scalar", "scalar"};
  std::vector<std::string> p_s = {"p", "s"};
  std::vector<std::string> cut = {"cut", "no_cut"};

  std::vector<std::string> filenames_reweighted_heavyhiggs;

  //for (int i = 0; i < mass_width_combis.size(); i++){
    //for (int j = 0; j < 1; j++){
  for (int j = 0; j < pseudo_scalar.size(); j++){
    for (int k = 0; k < cut.size(); k++){
      
      filenames_reweighted_heavyhiggs.clear();
      
      for (int i = 0; i < mass_width_combis.size(); i++){
         filenames_reweighted_heavyhiggs.push_back("/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_reweighting_" + p_s[j] + "_" + r_i[0] + "_" + mass_width_combis[i] + "_reweighted_juan_paper_fixed_width/" + "hist_heavyhiggs_reweighting_" + pseudo_scalar[j] + "_" + mass_width_combis[i] + "_juan_paper_" + resonance_interference[0]+ "_" + cut[k] + "_after_reordering_" + res_int[0] + "_reweighting_only.root");  // _different_beta_in_QCD_and_BSM.root");
         std::cout << "filiename_reweighted_heavyhiggs: " << filenames_reweighted_heavyhiggs[i] << std::endl;
        }
 
      std::string filename_reweighted_ttbarlo = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[0] + "_m400_w20_reweighted_juan_paper_fixed_width/" + "hist_ttbarlo_reweighting_" + pseudo_scalar[j] + "_m400_w20_juan_paper_" + resonance_interference[0] +  "_" + cut[k] + "_after_reordering_std_reweighting_only_2.root"; //_different_beta_in_QCD_and_BSM.root";
      
      std::cout << "filename_reweighted_ttbarlo: " << filename_reweighted_ttbarlo << std::endl;
 
      std::string folder = pseudo_scalar[j] + "_" + res_int[0] + "_" + cut[k] +"_more_higgs_summed_vs_ttbarlo"; //_different_beta_in_QCD_and_BSM";
      if (mkdir(folder.c_str(), 0777) != 0){
        std::cout << "Couldn't create directory." << std::endl;
      }
      
      std::string filename_histogram = folder + "/" +  std::string("all") + std::string(".pdf");
      std::cout << "filename_histogram: " << filename_histogram << std::endl;
      
      DrawHistograms(filenames_reweighted_heavyhiggs, filename_reweighted_ttbarlo, folder, filename_histogram, resonance_interference[0], pseudo_scalar[j], mass_width_combis[0]);
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

