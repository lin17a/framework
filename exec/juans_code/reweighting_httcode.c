template <typename double>
double weight_ggHtt(const HTT_Input<double>& httInput) {
/*
struct HTT_Input {
       unsigned int HIGGS_OPTION; // 0: SCALAR only; 1: PSEUDOSCALAR only; 2: BOTH
       double MH; // mass of scalar Higgs
       double MA; // mass of pseudo-scalar Higgs
       double GH; // width of scalar Higgs (only used if WIDTH_OPTION<2)
       double GA; // width of pseudoscalar Higgs (only used if WIDTH_OPTION<2)
       double YTOP; // times the SM Httbar coupling (1/tan(beta) in 2HDM aligned scenarios); also used to determine widths for WIDTH_OPTION>=2; note that when WIDTH_OPTION<2 this parameter is interpreted as a scale factor w.r.t the effective ggH (ggA) SM Yukawa coupling, and we assume that only the H->ttbar (A->ttbar) decay channel has sizable contributions to the total H (A) width
       unsigned int WIDTH_OPTION; // 0: fixed widths taken from GH and GA; 1: sqrt(s) running widths taken from GH and GA; 2: determine FIXED widths from YTOP; 3: RUNNING widths determined from YTOP
       TLorentzVector P4gen_t[2]; // index 0: top; index 1: antitop
       TLorentzVector P4gen_d[2]; // index 0: antidown-type fermion from W+ decay; index 1: down-type fermion from W- decay
};
*/
      // Check that the chosen options are valid
      if (httInput.HIGGS_OPTION>2 || httInput.WIDTH_OPTION>3) {
            printf("What? HIGGS_OPTION=%d, WIDTH_OPTION=%d; returning weight=1 !!\n", httInput.HIGGS_OPTION, httInput.WIDTH_OPTION);
            return 1.;
      }

      double MH2_REF = pow(httInput.MH,2);
      double MA2_REF = pow(httInput.MA,2);
      double ytop2 = pow(httInput.YTOP,2);

      TLorentzVector P4gen_t1 = httInput.P4gen_t[0]; 
      TLorentzVector P4gen_t2 = httInput.P4gen_t[1]; 
      TLorentzVector P4gen_d1 = httInput.P4gen_d[0]; 
      TLorentzVector P4gen_d2 = httInput.P4gen_d[1]; 

      const double MT_REF = 172.5;
      const double MT2_REF = MT_REF*MT_REF;
      const double GF = 1.16637876e-5;
      const double pi = TMath::Pi();
      const double pi2 = pi*pi;
      const double sqrt2 = sqrt(2.);

      TLorentzVector P4gen_tt = P4gen_t1+P4gen_t2;
      double sqrts = P4gen_tt.M();
      double sqrts2 = sqrts*sqrts;
      double mt1 = httInput.P4gen_t[0].M();
      double mt2 = httInput.P4gen_t[1].M();
      double beta2 = (1-(mt1-mt2)*(mt1-mt2)/sqrts2)*(1-(mt1+mt2)*(mt1+mt2)/sqrts2);
      //double beta = sqrt(beta2);
      //double beta3 = beta*beta2;
      double mtrefOverE = 2*MT_REF/sqrts;
      double beta2Ref = 1-mtrefOverE*mtrefOverE;

      TVector3 zBoost(0.,0.,-P4gen_tt.Pz()/P4gen_tt.E()); P4gen_tt.Boost(zBoost);
      TVector3 transverseBoost(-P4gen_tt.Px()/P4gen_tt.E(),-P4gen_tt.Py()/P4gen_tt.E(),0.); P4gen_tt.Boost(transverseBoost);
      P4gen_t1.Boost(zBoost); P4gen_t1.Boost(transverseBoost);
      P4gen_d1.Boost(zBoost); P4gen_d1.Boost(transverseBoost);
      P4gen_t2.Boost(zBoost); P4gen_t2.Boost(transverseBoost);
      P4gen_d2.Boost(zBoost); P4gen_d2.Boost(transverseBoost);
      double z1 = P4gen_t1.CosTheta();
      double z2 = z1*z1;

      double minusPhi1Plane = -P4gen_t1.Phi();
      TVector3 zVect(0.,0.,1.);
      P4gen_t1.Rotate(minusPhi1Plane,zVect);
      P4gen_d1.Rotate(minusPhi1Plane,zVect);
      P4gen_t2.Rotate(minusPhi1Plane,zVect);
      P4gen_d2.Rotate(minusPhi1Plane,zVect);
      double t1Yangle = -atan2(P4gen_t1.Px(),P4gen_t1.Pz());
      TVector3 yVect(0.,1.,0.);
      P4gen_t1.Rotate(t1Yangle,yVect);
      P4gen_d1.Rotate(t1Yangle,yVect);
      P4gen_t2.Rotate(t1Yangle,yVect);
      P4gen_d2.Rotate(t1Yangle,yVect);

      TVector3 t1Boost(-P4gen_t1.Px()/P4gen_t1.E(),-P4gen_t1.Py()/P4gen_t1.E(),-P4gen_t1.Pz()/P4gen_t1.E()); 
      P4gen_d1.Boost(t1Boost);
      TLorentzVector P4gen_d1_trest = P4gen_d1;
      TVector3 t2Boost(-P4gen_t2.Px()/P4gen_t2.E(),-P4gen_t2.Py()/P4gen_t2.E(),-P4gen_t2.Pz()/P4gen_t2.E()); 
      P4gen_d2.Boost(t2Boost);
      TLorentzVector P4gen_d2_trest = P4gen_d2;

      double common_factor_for_M2_gg_QCD = pi2/12*(7+9*beta2*z2)/pow(1-beta2*z2,2);
      double common_factor_for_M2_gg_ss_nointerf = pi2/12*(5+9*beta2*z2)/pow(1-beta2*z2,2);
      double common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta2*z2);

      using namespace std::literals::complex_literals;
      using namespace std::literals;
      const std::complex<double> i_cmplx = std::complex<double>(0., 1.L);

      std::complex<double> common_factor_for_M_H = 0.;
      double mh_gh_partial = 0.;
      double mh_gh = 0.;
      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
            if (httInput.MH>2*MT_REF) {
                  double beta3_tdecay = pow(1-4*MT2_REF/MH2_REF,1.5);
                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
                        mh_gh_partial += 3*GF*MT2_REF*sqrts2/(4*pi*sqrt2)*beta3_tdecay*ytop2;
                  } else {
                        mh_gh_partial += 3*GF*MT2_REF*MH2_REF/(4*pi*sqrt2)*beta3_tdecay*ytop2;
                  }
            }

            if (httInput.WIDTH_OPTION==1) {
                  mh_gh = httInput.GH*sqrts2/httInput.MH;
            } else if (httInput.WIDTH_OPTION==0) {
                  mh_gh = httInput.GH*httInput.MH;
            } else {
                  mh_gh = mh_gh_partial;
            }

            double NB_real, NB_imag;
            double auxh = 1.5*(1-beta2Ref)*ytop2;
            if (beta2Ref>0) {
                  double betaRef = sqrt(beta2Ref);
                  NB_real = auxh*(1-beta2Ref/4*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2));
                  NB_imag = auxh*(beta2Ref/2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  NB_real = auxh*(1+beta2Ref*pow(asin(sqrts/2/MT_REF),2));
                  NB_imag = 0.;
            }
            std::complex<double> NB = NB_real + NB_imag*i_cmplx;
            std::complex<double> denomH = sqrts2 - MH2_REF + mh_gh*i_cmplx;

            common_factor_for_M_H = NB / denomH;
            common_factor_for_M_H *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
      }

      std::complex<double> common_factor_for_M_A = 0.;
      double ma_ga_partial = 0.;
      double ma_ga = 0.;
      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
            if (httInput.MA>2*MT_REF) {
                  double beta_tdecay = pow(1-4*MT2_REF/MA2_REF,0.5);
                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
                        ma_ga_partial += 3*GF*MT2_REF*sqrts2/(4*pi*sqrt2)*beta_tdecay*ytop2;
                  } else {
                        ma_ga_partial += 3*GF*MT2_REF*MA2_REF/(4*pi*sqrt2)*beta_tdecay*ytop2;
                  }
            }
            if (httInput.WIDTH_OPTION==1) {
                        ma_ga = httInput.GA*sqrts2/httInput.MA;
            } else if (httInput.WIDTH_OPTION==0) {
                        ma_ga = httInput.GA*httInput.MA;
            } else {
                        ma_ga = ma_ga_partial;
            }

            double PB_real, PB_imag;
            double auxa = -1.5*(1-beta2Ref)/4.*ytop2;
            if (beta2Ref>0) {
                  double betaRef = sqrt(beta2Ref);
                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2);
                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  PB_real = -4*auxa*pow(asin(sqrts/2/MT_REF),2);
                  PB_imag = 0.;
            }
            std::complex<double> PB = PB_real + PB_imag*i_cmplx;
            std::complex<double> denomA = sqrts2 - MA2_REF + ma_ga*i_cmplx;

            common_factor_for_M_A = PB / denomA;
            common_factor_for_M_A *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
      }

      // Printing global input information (once)
      static bool printOnce = true;
      if (printOnce) {
            printOnce = false;
            printf("\nHTT_INPUT >>>> (printed only once)\n");
            printf("\tYTOP=%.3f times the SM Higgs-ttbar coupling\n", httInput.YTOP);
            if (httInput.HIGGS_OPTION==0) {
                  printf("\tOnly SCALAR CASE, MH=%.3f GeV\n", httInput.MH);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed width: GH=%.3f GeV\n", httInput.GH);
                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning width: GH(MH)=%.3f GeV\n", httInput.GH);
                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/sqrts2*httInput.MH);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed width: GH(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", mh_gh/sqrts2*httInput.MH);
                  }
            } else if (httInput.HIGGS_OPTION==1) {
                  printf("\tOnly PSEUDO-SCALAR CASE, MA=%.3f GeV\n", httInput.MA);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed width: GA=%.3f GeV\n", httInput.GA);
                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning width: GA(MA)=%.3f GeV\n", httInput.GA);
                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/sqrts2*httInput.MA);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed width: GA(top+bottom)=%.3f GeV\n", ma_ga/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", ma_ga/sqrts2*httInput.MA);
                  }
            } else {
                  printf("\tSCALAR + PSEUDO-SCALAR CASE, MH=%.3f GeV, MA=%.3f GeV\n", httInput.MH, httInput.MA);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed widths: GH=%.3f GeV, GA=%.3f GeV\n", httInput.GH, httInput.GA);
                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH, ma_ga_partial/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning widths: GH(MH)=%.3f GeV, GA(MA)=%.3f GeV\n", httInput.GH, httInput.GA);
                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/sqrts2*httInput.MH, ma_ga_partial/sqrts2*httInput.MA);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH, ma_ga/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/sqrts2*httInput.MH, ma_ga/sqrts2*httInput.MA);
                  }
            }
            printf("\n");
      }

      double c1 = P4gen_d1_trest.CosTheta();
      double c2 = P4gen_d2_trest.CosTheta();
      double s1 = sqrt(1-c1*c1);
      double s2 = sqrt(1-c2*c2);
      double phi1 = P4gen_d1_trest.Phi();
      double phi2 = P4gen_d2_trest.Phi();
      double os_factor = (1+z2+(1-z2)*(1-beta2Ref))*(1.-c1*c2) + 2*(1-beta2Ref)*(1-z2)*c1*c2
            - beta2Ref*(1-z2)*s1*s2*cos(phi1+phi2)
            -2*(1-beta2Ref)*(1-z2)*s1*s2*cos(phi1)*cos(phi2)
            +2*sqrt(1-beta2Ref)*z1*sqrt(1-z2)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
      double ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
      double ss_factor_beta = beta2 * (1+c1*c2 - s1*s2*cos(phi1-phi2));

      double M2_QCD = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
                     + common_factor_for_M2_gg_QCD*(1-beta2Ref)*(ss_factor_1+ss_factor_beta);

      double M2 = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
                + common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
                + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
                + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;
      //double M2 = common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
      //             + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
      //             + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;

      // Final weight
      double weight = M2 / M2_QCD;

      std::cout << "printf(\"M2_QCD = " << M2_QCD << "\n\");" << std::endl;
      std::cout << "printf(\"M2_BSM = " << M2 << "\n\");" << std::endl;
      std::cout << "printf(\"weight = " << weight << "\n\");" << std::endl;

      if (isnan(weight)) {
            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
            return 0.;
      }

      return weight;
}
