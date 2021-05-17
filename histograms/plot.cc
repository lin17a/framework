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
      hist->hist->Scale(1. / hist->hist->Integral());
      int nbins = hist->hist->GetXaxis()->GetNbins();
      for (int i=0; i<nbins ;i++){
        Float_t v = hist->hist->GetBinContent(i);
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
    while ((key=(TKey*)iter.Next()) ) {  
       set_hist_names.insert(key->GetName());
    }
    if (set_hist_names.find(hist_name) != set_hist_names.end()){
      h1 = (TH1F *) file->Get(hist_name.c_str());
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
        status = false;
      }
    }
    histogram_raw hist_struct = {h1, status};
    return hist_struct; 
  };

  auto MakeTitle = []( std::string hist_name, std::string res_int, std::string higgs_type, std::string mass, std::string width){
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

          std::string variable;
          if (hist_name == "ttbar_mass_no_cut"){
            variable = "m_{t#bar{t}}";
          }

          std::string hist_title = variable + " for " + higgs_name + " heavyhiggs, " + res_int + ", m_{" + subscript + "}=" + mass + " GeV, " + "#Gamma_{" + subscript + "}=" + width + " GeV";
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
    {"ttbar_mass_no_cut", "m_{t#bar{t}}"}
  };

  auto DrawHistograms = [ SetUpHistograms, GetHistogram, MakeTitle, units ](std::string filename_generated, std::string filename_reweighted, std::string folder, std::string filename_histogram, std::string res_int, std::string higgs_type, std::string mass, std::string width){
   
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
       TH1F *h1 = (TH1F *) f1->Get( hist_name.c_str());   
       std::string title = hist_name;     
       //size_t pos_1 = hist_name.find("_no_cut");
       //size_t pos_2 = hist_name.find("_no_cut", pos_1 + 1);
       //std::string attr = hist_name.substr(0, pos_1);
       //std::string hist_name_positive = attr + "_positive";
       //std::cout << hist_name_positive << std::endl;
      
      // try to find histograms with the same name 
      try {

          std::cout << hist_name << std::endl;
          
          // get histograms with the same attribute from the other files 
          histogram_raw hist_struct = GetHistogram(f3, hist_name);   
          TH1F *h3;
          if (hist_struct.status) {
            h3 = hist_struct.hist;
          }
          else {
            continue;
          }
          // defining histogram structs with a histogram, color and name to hand it all together to the SetUp function
          histogram hist1 = {h1, "generated", kGreen-3};
        //  histogram hist2 = {h2, "int pseudo", kGreen+3};
          histogram hist3 = {h3, "reweighted", kBlue+2};
          
          // add a new histogram containing the sum of ttbar and the pseudo scalar higgs histograms
          //TH1F *h_sum_res_int_pseudo = (TH1F*) h1->Clone();
        //  h_sum_res_int_pseudo->Add(h2);
          //h_sum_ttbarlo_res_int_pseudo->Add(h4);
          //h_sum_res_int_pseudo->SetLineColor(kAzure+3);

         // histogram hist_sum_pseudo = {h_sum_res_int_pseudo, "sum pseudo", kRed+2};

          // initialize legend, numbers are the coordinates of the bottom left corner und upper right corner
          // fits exactly under the info panel
          //auto legend = new TLegend(0.78,0.69,0.98,0.77);
          auto legend = new TLegend(0.78,0.79,0.98,0.87);
          legend->SetTextSize(0.02);        
  
          std::vector<histogram *> hists{&hist1, &hist3};

          // adjust colors, scales, ... 
          SetUpHistograms(&hists, legend, 2, 1);
          
          // find and set the right title for the x axis
          size_t pos_1 = hist_name.find('_');
          size_t pos_2 = hist_name.find("_no_cut", pos_1 + 1);
          std::string attr = hist_name.substr(pos_1 + 1, pos_2 - pos_1 - 1);
          std::string x_title; 
          auto unit = units.find(attr);
          if (hist_name == "ttbar_mass_no_cut"){
            attr = "m_{t#bar{t}}";
          }
          if (unit != units.end()) {
            x_title = attr + " [" + unit->second + "]"; 
            }
          else {
            x_title = attr;
          }
 
          
          h1->GetXaxis()->SetTitle(x_title.c_str());
          h1->GetYaxis()->SetTitle("normalized number of events");
          //h1->GetYaxis()->CenterTitle(true);
          h1->SetTitleSize(.02, "XY");
          h1->SetTitleSize(.02, "t");

          std::string hist_title = MakeTitle(hist_name, res_int, higgs_type, mass, width);

          h1->SetNameTitle(hist_name.c_str(), hist_title.c_str());
          h1->SetLabelSize(.02, "XY");
          //gStyle->SetTitleFontSize(20);
          
          TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
          //TCanvas can("canvas",
          if (attr == "mass") 
             can->SetLogy();
          if (hist_name == "ttbar_mass_no_cut") 
             can->SetLogy();
          can->cd();
             
          
          // draw all in the same histogram
          h1->Draw("pe");
//          h2->Draw("pesame");
          h3->Draw("pesame");
          //gStyle->SetTitleFontSize(0.7);

          can->Update();
          can->Clear();
         
          gStyle->SetLabelSize(.02, "XY");
          //gStyle->SetTitleFontSize(.11);
          gStyle->SetTitleSize(0.1, "t");
          gStyle->SetTitleSize(.02, "XY");
          gStyle->SetOptStat(0);
 
          // add a sub pad with the ratio 
          auto rp1 = new TRatioPlot(h1, h3);
          rp1->Draw();
          rp1->GetLowerRefYaxis()->SetTitle("#frac{generated}{reweighted}");
          rp1->GetLowerRefYaxis()->SetRangeUser(0,2);
          rp1->GetLowerRefYaxis()->SetTitleOffset(2);
          rp1->GetLowerRefYaxis()->CenterTitle(true);
          //rp1->GetLowerRefXaxis()->SetTitleSize(20);
          rp1->GetXaxis()->SetLabelSize(0.015);
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
  std::vector<std::string> masses = { "400", "1000", "600", "800" };
  std::vector<std::string> widths = { "20", "25", "30", "20" };
  std::vector<std::string> zeros = { "00", "000", "00", "000"};
  std::vector<std::string> res_int =  {"res"}; //, "int"}; //
  std::vector<std::string> RES_INT = {"RES"}; //, "INT"}; //
  std::vector<std::string> r_i = {"r"}; //, "i"}; //
  std::vector<std::string> resonance_interference = {"resonance"}; //, "interference"}; //
  std::vector<std::string> PSEUDO_SCALAR = {"PSEUDO", "SCALAR"}; 
  std::vector<std::string> pseudo_scalar = {"pseudo_scalar", "scalar"};
  std::vector<std::string> p_s = {"p", "s"};

  std::vector<std::string> positive_negative; 


  for (int i = 0; i < mass_width_combis.size(); i++){
    for (int j = 0; j < pseudo_scalar.size(); j++){
      for (int l = 0; l < res_int.size(); l++){

        if (res_int[l] == "int"){
          positive_negative = {"_positive", "_negative"}; 
        }
        else {
          positive_negative = {""};
        }

        for (int k = 0; k < positive_negative.size(); k++){
 

          std::string filename_generated = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_m" + masses[i] + "_w" + widths[i] + "_" + zeros[i]  + "_"+ RES_INT[l] +"_" + PSEUDO_SCALAR[j] + "_generated/hist_m" + masses[i] + "_w" + widths[i] + "_" + zeros[i] +  "_" + RES_INT[l] + "_" + PSEUDO_SCALAR[j] + "_no_cut" + positive_negative[k] + "_with_z.root";

          std::cout << "filename_generated: " << filename_generated << std::endl;

          std::string filename_reweighted = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[l] + "_m" + masses[i] + "_w" + widths[i] + "_reweighted_juan_paper_fixed_width/hist_ttbarlo_reweighting_" + pseudo_scalar[j] + "_m" + masses[i] + "_w" + widths[i] + "_juan_paper_" + resonance_interference[l] + "_no_cut" + positive_negative[k] + "_after_reordering_with_z.root";
          //std::string filename_reweighted = "/nfs/dust/cms/user/meyerlin/ba/framework/runs/" + p_s[j] + "_" + r_i[l] + "_" + mass_width_combis[i] + "_reweighted_juan_paper_fixed_width/hist_ttbarlo_reweighting_" + pseudo_scalar[j] + "_" + mass_width_combis[i]  + "_juan_paper_" + resonance_interference[l] + "_no_cut_positive_after_reordering_running_width_with_s.root";

          std::cout << "filename_reweighted: " << filename_reweighted << std::endl;
 
          std::string folder =  "m" + masses[i] + "_w" + widths[i] + "_" + pseudo_scalar[j] + "_" + res_int[l] + positive_negative[k] + "_original_generated_vs_reweighted";
          //std::string folder = "test"; 
          if (mkdir(folder.c_str(), 0777) != 0){
            std::cout << "Couldn't create directory." << std::endl;
          }
          
          std::string filename_histogram = folder + "/" +  std::string("all") + std::string(".pdf");
          std::cout << "filename_histogram: " << filename_histogram << std::endl;
          
          DrawHistograms(filename_generated, filename_reweighted, folder, filename_histogram, resonance_interference[l], pseudo_scalar[j], masses[i], widths[i]);
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

