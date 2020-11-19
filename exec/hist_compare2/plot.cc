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


  // no_cut histogram
  TFile *f1 = TFile::Open( "./hist_RES_PSEUDO_no_cut.root" );
  // cut histogram
  TFile *f2 = TFile::Open( "./hist_RES_PSEUDO_cut.root" );

  TIter iter1(f1->GetListOfKeys());
  TKey *key1;

  while ((key1=(TKey*)iter1.Next()) ) {
    
     TH1F *h1 = (TH1F *) f1->Get( key1->GetName());
     std::string hist_name1 = key1->GetName();
     std::string no_cut = "_no_cut";
     size_t pos = hist_name1.find(no_cut);
     if (pos != std::string::npos)
     {
        // If found then erase it from string
        hist_name1.erase(pos, no_cut.length());
     }     
     std::string hist_name2 = hist_name1 + std::string("_cut");
     std::string title = hist_name1 + std::string("_cut_vs_no_cut");
     //std::cout << hist_name1 << "   "<< hist_name2 << "\n";
     
     try {
	
        TH1F *h2 = (TH1F *) f2->Get( hist_name2.c_str());   
        h1->SetLineColor(kPink - 1);
        h2->SetLineColor(kAzure - 1); 
        
        size_t pos_1 = hist_name1.find('_');
	size_t pos_2 = hist_name1.find("_RES_PSEUDO", pos_1 + 1);
        std::string attr = hist_name1.substr(pos_1 + 1, pos_2 - pos_1);
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

        double_t max_h1 = h1->GetMaximum();
        double_t max_h2 = h2->GetMaximum();
	if (max_h1 > max_h2 ) 
		h2->SetMaximum(max_h1 + max_h1 / 8);
        else
                h1->SetMaximum(max_h2 + max_h2 / 8);

	TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
        can->cd();
           
        h1->Draw("pe");
        h2->Draw("pesame");
    
        auto legend = new TLegend(0.78,0.69,0.98,0.77);
        //legend->SetHeader("Legend","C"); // option "C" allows to center the header
        legend->AddEntry(h1,"no cut","l");
        legend->AddEntry(h2,"cut","l");
        legend->Draw();
 
        std::string filename = std::string(hist_name1) + std::string("_RES_PSEUDO_cut_vs_no_cut.pdf");
        can->SaveAs( filename.c_str() ); 
        }
     
     catch (const std::exception&) { 
	 std:cout << "exception gecatched\n"; }
    }
}
