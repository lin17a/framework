{
 TFile *f1 = TFile::Open( "hist_no_cut.root" );
 TFile *f2 = TFile::Open( "hist_cut.root" );

 TH1F *h1 = (TH1F *) f1->Get( "angle_no_cut" );
 TH1F *h2 = (TH1F *) f2->Get( "angle_cut" );

 h1->SetLineColor(kPink - 1);
 h2->SetLineColor(kAzure - 1);

 TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
 can->cd();
 h1->Draw("pe");
 h2->Draw("pesame");

 can->SaveAs("ttbar_angle_cut_vs_no_cut.pdf");
}
