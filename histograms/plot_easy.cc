
#include<algorithm>
#include<vector>
//#include <sys/stat.h>


{
 std::vector<std::string> mass_width_combis = { "m400_w20_00", "m1000_w25_000", "m600_w30_00", "m800_w20_000" };
 //for (int i=0; i < mass_width_combis.size(); i++){
   //TFile *f1 = TFile::Open( "/nfs/dust/cms/user/meyerlin/ba/framework/runs/p_r_m400_w20_reweighted_juan_paper_fixed_width/hist_ttbarlo_no_reweighting_pseudo_scalar_m400_w20_juan_paper_resonance_no_cut_after_reordering_with_different_spins.root" );
   TFile *f1 = TFile::Open( "/nfs/dust/cms/user/meyerlin/ba/framework/runs/s_r_m800_w20_reweighted_juan_paper_fixed_width/hist_ttbarlo_no_reweighting_scalar_m800_w20_juan_paper_resonance_no_cut_after_reordering_with_different_spins.root" );
   //std::cout << mass_width_combis[0] << std::endl; 
  
   //std::string filename =  "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_" + mass_width_combis[i] + "_RES_PSEUDO_generated/hist_" + mass_width_combis[i] + "_RES_PSEUDO_no_cut_with_z.root";
   //std::cout << filename << std::endl;
   //TFile *f1 = TFile::Open( "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_" + mass_width_combis[i] + "_RES_PSEUDO_generated/hist_" + mass_width_combis[i] + "_RES_PSEUDO_no_cut_with_z.root" );
   //TFile *f1 = TFile::Open( "/nfs/dust/cms/user/meyerlin/ba/framework/runs/heavyhiggs_m800_w20_000_RES_SCALAR_generated/hist_m800_w20_000_RES_SCALAR_no_cut_with_z.root" );
   //TFile *f2 = TFile::Open( "hist_cut.root" );
  
   TH1F *h1 = (TH1F *) f1->Get( "t_CosTheta_LHE_no_cut" );
   TH1F *h2 = (TH1F *) f1->Get( "cpTP_LHE_no_cut" );
   TH1F *h_diff = (TH1F*) h1->Clone();
   h_diff->Add(h2, -1.);
  
   h1->SetLineColor(kGreen-3);
   h2->SetLineColor(kAzure - 1);
   h_diff->SetLineColor(kRed);
  
   h1->Scale(1. / h1->Integral());
   h2->Scale(1. / h2->Integral());
   //h_diff->Scale(1. / h_diff->Integral());
  
   h1->SetLineWidth(2);
   h2->SetLineWidth(2);
   h_diff->SetLineWidth(2);
         
   h1->GetXaxis()->SetTitle("");
   h1->GetYaxis()->SetTitle("number of events");
  
   auto legend = new TLegend(0.78,0.69,0.98,0.77);
   //auto legend = new TLegend(0.78,0.76,0.98,0.83);
   //legend->AddEntry(h1, "t_CosTheta LHE", "l");
   //legend->AddEntry(h2, "cpTP LHE", "l");
   legend->AddEntry(h_diff, "cpTP - t_cos_theta, LHE", "l");
  
   //h1->SetNameTitle("t_cos_theta vs cpTP LHE", "SM scalar t_CosTheta vs cpTP, LHE,  m=800, w=20");
   h_diff->SetNameTitle("t_cos_theta - cpTP ", "SM scalar t_CosTheta minus cpTP, LHE,  m=800, w=20");
  
   TCanvas *can = new TCanvas("canvas", "canvas", 200, 10, 1000, 1000);
   can->cd();
   //h1->Draw("pe");
   //h2->Draw("pesame");
   h_diff->Draw("pe");
  
   //can->Clear();
  
   gStyle->SetTitleFontSize(.05);
   gStyle->SetLabelSize(.025, "XY");
   gStyle->SetTitleSize(.03, "XY");
  
  
 //  auto rp1 = new TRatioPlot(h1, h2);
 //  rp1->Draw();
 //  rp1->GetLowerRefYaxis()->SetTitle("t_cos_theta / cpTP");
 // // rp1->GetLowerRefYaxis()->SetRangeUser(0,2); 
   can->Update();
  
   legend->Draw();
  
   can->SaveAs("t_CosTheta_minus_cpTP_LHE_s_r_SM_m800_w20_no_cut.pdf");
   delete can;
   //delete rp1;

}
