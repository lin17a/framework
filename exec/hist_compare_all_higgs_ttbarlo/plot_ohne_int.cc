#include<algorithm>
#include<vector>

{
  auto SetUpHistogram = [](TH1F *hist, TLegend *legend, std::string name, Color_t line_color, double_t max_value, double_t min_value, int line_width = 3, int rebin = 2) {
    hist->SetLineColor(line_color);
    hist->SetLineWidth(line_width);
    hist->Rebin(rebin);
    hist->Scale(1. / hist->Integral());
    hist->SetMaximum(max_value + max_value/ 8);
    std::cout << std::to_string(min_value) << std::lineend;    
//hist->SetMinimum(min_value);
    legend->AddEntry(hist, name.c_str(), "l");
  };


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

  // add all root files
  TFile *f1 = TFile::Open( "./hist_ttbarlo_spin_no_cut.root" );
  TFile *f2 = TFile::Open( "./hist_RES_PSEUDO_spin_no_cut.root" );
  TFile *f3 = TFile::Open( "./hist_RES_SCALAR_spin_no_cut.root" );
  //TFile *f4 = TFile::Open( "./hist_INT_PSEUDO_spin_no_cut.root" );
  //TFile *f5 = TFile::Open( "./hist_INT_SCALAR_spin_no_cut.root" );

  TIter iter1(f1->GetListOfKeys());
  TKey *key1;

  // iterate over all histograms in the first file
  while ((key1=(TKey*)iter1.Next()) ) {
    
     std::string hist_name = key1->GetName();
     // get the first histogram from the first file
     TH1F *h1 = (TH1F *) f1->Get(hist_name.c_str());
     std::string title = hist_name;      

    // try to find histograms with the same name 
    try {
	
        // get histograms with the same attribute from the other files 
        TH1F *h2 = (TH1F *) f2->Get( hist_name.c_str());   
        TH1F *h3 = (TH1F *) f3->Get( hist_name.c_str());   
        //TH1F *h4 = (TH1F *) f4->Get( hist_name.c_str());   
        //TH1F *h5 = (TH1F *) f5->Get( hist_name.c_str());   
        
        // initialize legend, numbers are the coordinates of the bottom left corner und upper right corner
        // fits exactly under the info panel
        auto legend = new TLegend(0.78,0.69,0.98,0.77);
        // find the maximum of the maxima to adjust the y axis 
        std::vector<double_t> maxs{h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum(), h4->GetMaximum(), h5->GetMaximum()};
        double_t max_value = *( std::max_element(maxs.begin(), maxs.end()) ); 
        // find the minimum of the minima to adjust the y axis 
        std::vector<double_t> mins{h1->GetMinimum(), h2->GetMinimum(), h3->GetMinimum(), h4->GetMinimum(), h5->GetMinimum()};
        double_t min_value = *( std::min_element(mins.begin(), mins.end()) );
        // adjust colors, scales, ... 
        SetUpHistogram(h1, legend, "ttbarlo", kRed+2, max_value, min_value);
        SetUpHistogram(h2, legend, "res pseudo", kGreen+3, max_value, min_value);
        SetUpHistogram(h3, legend, "res scalar", kGreen-3, max_value, min_value);
        //SetUpHistogram(h4, legend, "int pseudo", kBlue+2, max_value, min_value);
        //SetUpHistogram(h5, legend, "int scalar", kAzure+1, max_value, min_value);
        
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
 
        h1->SetNameTitle(hist_name.c_str(), "");
        
        TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
        can->cd();
           
        // draw all in the same histogram
        h1->Draw("pe");
        h2->Draw("pesame");
        h3->Draw("pesame");
        //h4->Draw("pesame");
        //h5->Draw("pesame");
    
        legend->Draw();
 
        // save with correct name
        std::string filename = std::string(hist_name) + std::string("_ttbarlo_all_higgs.pdf");
        can->SaveAs( filename.c_str() ); 
     }
        
     catch (const std::exception&) { 
	 std:cout << "exception gecatched\n"; 
     }
  }
}

