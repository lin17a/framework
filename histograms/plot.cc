#include<algorithm>
#include<vector>
#include <sys/stat.h>

{
  struct histogram {
    TH1F *hist;
    std::string name;
    Color_t line_color;
  };
  
  auto SetUpHistograms = [](std::vector<histogram *> *hists, TLegend *legend, int line_width = 3, int rebin = 2) {
    for (histogram *hist : *hists) {
      hist->hist->SetLineColor(hist->line_color);
      hist->hist->SetLineWidth(line_width);
      hist->hist->Rebin(rebin);
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

  auto DrawHistograms = [ SetUpHistograms, units ](std::string filename_generated, std::string filename_reweighted, std::string folder, std::string filename_histogram, std::string res_int, std::string higgs_type, std::string mass_width){
   
    TFile *f1 = TFile::Open( filename_generated.c_str() );
    //TFile *f2 = TFile::Open( "hist_ttbarlo_reweighting_pseudo_scalar_m1000.000000_w25.000000_juan_paper_resonance_cut_after_reordering.root" );
    TFile *f3 = TFile::Open( filename_reweighted.c_str() );

    if (!f1->IsOpen() or !f3->IsOpen()){
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

      // try to find histograms with the same name 
      try {
          
          // get histograms with the same attribute from the other files 
      //    TH1F *h2 = (TH1F *) f2->Get( hist_name.c_str());   
          TH1F *h3 = (TH1F *) f3->Get( hist_name.c_str());   

          // defining histogram structs with a histogram, color and name to hand it all together to the SetUp function
          histogram hist1 = {h1, "generated", kGreen-3};
        //  histogram hist2 = {h2, "int pseudo", kGreen+3};
          histogram hist3 = {h3, "reweighted", kBlue+2};
          
          // add a new histogram containing the sum of ttbar and the pseudo scalar higgs histograms
          TH1F *h_sum_res_int_pseudo = (TH1F*) h1->Clone();
        //  h_sum_res_int_pseudo->Add(h2);
          //h_sum_ttbarlo_res_int_pseudo->Add(h4);
          h_sum_res_int_pseudo->SetLineColor(kAzure+3);

         // histogram hist_sum_pseudo = {h_sum_res_int_pseudo, "sum pseudo", kRed+2};

          // initialize legend, numbers are the coordinates of the bottom left corner und upper right corner
          // fits exactly under the info panel
          auto legend = new TLegend(0.78,0.69,0.98,0.77);
          
          std::vector<histogram *> hists{&hist1, &hist3};

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
          //TCanvas can("canvas", 
          can->SetLogy();
          can->cd();
             
          // draw all in the same histogram
          h1->Draw("pe");
//          h2->Draw("pesame");
          h3->Draw("pesame");

          can->Clear();
          
          // add a sub pad with the ratio 
          auto rp1 = new TRatioPlot(h1, h3);
          rp1->Draw();
          rp1->GetLowerRefYaxis()->SetTitle("ratio");
          rp1->GetLowerRefYaxis()->SetRangeUser(0,2);
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





  std::vector<std::string> mass_width_combis = { "m400_w20", "m1000_w25", "m600_w30", "m800_w20" };
  std::vector<std::string> zeros = { "00", "000", "00", "000"};
  std::vector<std::string> res_int = {"res"}; //, "int"};
  std::vector<std::string> RES_INT = {"RES"}; //, "INT"};
  std::vector<std::string> r_i = {"r"}; //, "i"};
  std::vector<std::string> resonance_interference = {"resonance"}; //, "interference"};
  std::vector<std::string> PSEUDO_SCALAR = {"PSEUDO", "SCALAR"};
  std::vector<std::string> pseudo_scalar = {"pseudo_scalar", "scalar"};
  std::vector<std::string> p_s = {"p", "s"};


  for (int i = 0; i < mass_width_combis.size(); i++){
    for (int j = 0; j < pseudo_scalar.size(); j++){
  
      std::string filename_generated = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_" + mass_width_combis[i] + "_" + zeros[i]  + "_"+ RES_INT[0] +"_" + PSEUDO_SCALAR[j] + "_generated/hist_" + mass_width_combis[i] + "_" + zeros[i] +  "_" + RES_INT[0] + "_" + PSEUDO_SCALAR[j] + "_cut.root";
      std::cout << "filname_generated: " << filename_generated << std::endl;
      std::string filename_reweighted = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[0] + "_" + mass_width_combis[i] + "_reweighted/hist_ttbarlo_reweighting_" + pseudo_scalar[j] + "_" + mass_width_combis[i]  + "_juan_paper_" + resonance_interference[0] + "_cut_after_reordering.root";
      std::cout << "filname_reweighted: " << filename_reweighted << std::endl;
 
      std::string folder = mass_width_combis[i] + "_" + pseudo_scalar[j] + "_" + res_int[0] + "_generated_vs_reweighted";
      if (mkdir(folder.c_str(), 0777) != 0){
        std::cout << "Couldn't create directory." << std::endl;
      }
      
      std::string filename_histogram = folder + "/" +  std::string("all") + std::string(".pdf");
      std::cout << "filename_histogram: " << filename_histogram << std::endl;
      
      DrawHistograms(filename_generated, filename_reweighted, folder, filename_histogram, resonance_interference[0], pseudo_scalar[j], mass_width_combis[i]);
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

