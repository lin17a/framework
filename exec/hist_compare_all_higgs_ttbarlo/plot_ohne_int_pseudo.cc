#include<algorithm>
#include<vector>

{
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
//  TFile *f4 = TFile::Open( "./hist_INT_PSEUDO_spin_no_cut.root" );
  TFile *f5 = TFile::Open( "./hist_INT_SCALAR_spin_no_cut.root" );

  TIter iter1(f1->GetListOfKeys());
  TKey *key1;

  while ((key1=(TKey*)iter1.Next()) ) {
    
     TH1F *h1 = (TH1F *) f1->Get( key1->GetName());
     std::string hist_name1 = key1->GetName();
     //std::string no_cut = "_no_cut";
     //size_t pos = hist_name1.find(no_cut);
//      if (pos != std::string::npos)
//      {
//         // If found then erase it from string
//         hist_name1.erase(pos, no_cut.length());
//      }     
//     std::string hist_name2 = hist_name1 + std::string("_cut");
     std::string title = hist_name1; // + std::string("_cut_vs_no_cut");
     //std::cout << hist_name1 << "   "<< hist_name2 << "\n";
     
     try {
	
        TH1F *h2 = (TH1F *) f2->Get( hist_name1.c_str());   
        TH1F *h3 = (TH1F *) f3->Get( hist_name1.c_str());   
//        TH1F *h4 = (TH1F *) f4->Get( hist_name1.c_str());   
        TH1F *h5 = (TH1F *) f5->Get( hist_name1.c_str());   
        h1->SetLineColor(kRed + 2);
        h2->SetLineColor(kGreen + 3); 
        h3->SetLineColor(kGreen - 3); 
//        h4->SetLineColor(kBlue + 2); 
        h5->SetLineColor(kAzure + 1); 
        
        size_t pos_1 = hist_name1.find('_');
	size_t pos_2 = hist_name1.find("_no_cut", pos_1 + 1);
        std::string attr = hist_name1.substr(pos_1 + 1, pos_2 - pos_1 - 1);
	std::cout << attr << " \n";
       
        std::string x_title; 
        auto unit = units.find(attr);
	if (unit != units.end()) {
          x_title = attr + " in " + unit->second; 
          //std::cout << unit->second << "\n"; 
          }
        else {
          x_title = attr;
        }
        
        std::cout << x_title << "\n";
        
        h1->GetXaxis()->SetTitle(x_title.c_str());
        h1->GetYaxis()->SetTitle("number of events");
 
        h1->SetNameTitle(hist_name1.c_str(), "");
        h1->Scale(1. / h1->Integral());
        h2->Scale(1. / h2->Integral());
        h3->Scale(1. / h3->Integral());
//        h4->Scale(1. / h4->Integral());
        h5->Scale(1. / h5->Integral());

        std::vector<double_t> maxs{h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum(), h5->GetMaximum()}; // h4->GetMaximum()};
	std::vector<double_t>::iterator max_p;
        max_p = std::max_element(maxs.begin(), maxs.end());
        double_t max_value = *max_p;

	h1->SetMaximum(max_value + max_value/ 8);
	h2->SetMaximum(max_value + max_value/ 8);
	h3->SetMaximum(max_value + max_value/ 8);
//	h4->SetMaximum(max_value + max_value/ 8);
	h5->SetMaximum(max_value + max_value/ 8);

	TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
        can->cd();
           
        h1->Draw("pe");
        h2->Draw("pesame");
        h3->Draw("pesame");
//        h4->Draw("pesame");
        h5->Draw("pesame");
    
        auto legend = new TLegend(0.78,0.60,0.98,0.77);
        //auto legend = new TLegend(0.78,0.69,0.98,0.77);
        //legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(h1,"ttbarlo","l");
        legend->AddEntry(h2,"res pseudo","l");
        legend->AddEntry(h3,"res scalar","l");
//        legend->AddEntry(h4,"int pseudo","l");
        legend->AddEntry(h5,"int scalar","l");
        legend->Draw();
 
        std::string filename = std::string(hist_name1) + std::string("_all_wo_int_pseudo.pdf");
        can->SaveAs( filename.c_str() ); 
        }
     
     catch (const std::exception&) { 
	 std:cout << "exception gecatched\n"; }
    }
}
