{
  TFile *f1 = TFile::Open( "hist_no_cut.root" );
  TFile *f2 = TFile::Open( "hist_cut.root" );

   TIter iter1(f1->GetListOfKeys());
   TKey *key1;
   TIter iter2(f2->GetListOfKeys());
   TKey *key2;

  //while ((key1=(TKey*)iter1.Next()) && (key2=(TKey*)iter2.Next())) {
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
     
     std::cout << hist_name1 << "   "<< hist_name2 << "\n";
     
     try {
        TH1F *h2 = (TH1F *) f2->Get( hist_name2.c_str());   
        h1->SetLineColor(kPink - 1);
        h2->SetLineColor(kAzure - 1); 
    
        TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
        can->cd();
           
        h1->Draw("pe");
        h2->Draw("pesame");
     
        std::string filename = std::string(hist_name1) + std::string("_cut_vs_no_cut.pdf");
        can->SaveAs( filename.c_str() ); }
     
     catch (const std::exception&) { 
	 std:cout << "exception gecatched"; }
    }
}
