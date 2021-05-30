#include <complex>
#include "TLorentzVector.h"
template <typename Number2>
struct HTT_Input {
       unsigned int HIGGS_OPTION; // 0: SCALAR only; 1: PSEUDOSCALAR only; 2: BOTH
       Number2 MH; // mass of scalar Higgs
       Number2 MA; // mass of pseudo-scalar Higgs
       Number2 GH; // width of scalar Higgs (only used if WIDTH_OPTION<2)
       Number2 GA; // width of pseudoscalar Higgs (only used if WIDTH_OPTION<2)
       Number2 YTOP; // times the SM Httbar coupling (1/tan(beta) in 2HDM aligned scenarios); also used to determine widths for WIDTH_OPTION>=2; note that when WIDTH_OPTION<2 this parameter is interpreted as a scale factor w.r.t the effective ggH (ggA) SM Yukawa coupling, and we assume that only the H->ttbar (A->ttbar) decay channel has sizable contributions to the total H (A) width
       unsigned int WIDTH_OPTION; // 0: fixed widths taken from GH and GA; 1: sqrt(s) running widths taken from GH and GA; 2: determine FIXED widths from YTOP; 3: RUNNING widths determined from YTOP
       TLorentzVector P4gen_t[2]; // index 0: top; index 1: antitop
       TLorentzVector P4gen_d[2]; // index 0: antidown-type fermion from W+ decay; index 1: down-type fermion from W- decay
};

// Reweighting gg-> H/A -> ttbar events
// Input is an HTT_Input structure; output is the event weight
template <typename Number2>
Number2 weight_ggHtt(const HTT_Input<Number2>& httInput);

// Simple code to decide whether the initial state should be gg or not
// pdgId1 is the the PDG id of first parton, along the +Z direction
// pdgId2 is the the PDG id of other parton, along the -Z direction
// pz_hard_radiation is the pz of the system made of the hard radiated partons recoiling 
template <typename Number2>
bool is_gg_initial_state(const int& pdgId1, const int& pdgId2, const Number2& pz_hard_radiation);

template <typename Number2>
Number2 weight_ggHtt(const HTT_Input<Number2>& httInput) {
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

      // define all important variables
      Number2 m_H_ref_sq = pow(httInput.MH,2);
      Number2 m_A_ref_sq = pow(httInput.MA,2);
      // coupling constant
      Number2 ytop_sq = pow(httInput.YTOP,2);

      TLorentzVector vec_4_t = httInput.P4gen_t[0]; 
      TLorentzVector vec_4_tbar = httInput.P4gen_t[1]; 
      TLorentzVector vec_4_lbar = httInput.P4gen_d[0]; 
      TLorentzVector vec_4_l = httInput.P4gen_d[1]; 

      const Number2 m_t_ref = 172.5;
      const Number2 m_t_ref_sq = m_t_ref*m_t_ref;
      //Fermi-Konstante
      const Number2 G_F = 1.16637876e-5;
      const Number2 alpha_s = 0.12;
      const Number2 pi = TMath::Pi();
      const Number2 pi_sq = pi*pi;
      const Number2 sqrt_2 = sqrt(2.);

      TLorentzVector vec_4_ttbar = vec_4_t+vec_4_tbar;
      Number2 sqrt_s = vec_4_ttbar.M();
      Number2 s = sqrt_s*sqrt_s;
      Number2 m_t = httInput.P4gen_t[0].M();
      Number2 m_tbar = httInput.P4gen_t[1].M();
      Number2 beta_sq = (1-(m_t-m_tbar)*(m_t-m_tbar)/s)*(1-(m_t+m_tbar)*(m_t+m_tbar)/s);
      Number2 beta = sqrt(beta_sq);
      //Number2 beta3 = beta*beta_sq;
      std::cout << "top mass: " << m_t << std::endl;
      Number2 mtrefOverE = 2 * m_t_ref / sqrt_s;
      std::cout << "square root of s: " << sqrt_s << std::endl;
      std::cout << "mt over E: " << mtrefOverE << std::endl; 
      Number2 beta_sq_ref = 1 - mtrefOverE * mtrefOverE;
      std::cout << "beta squared ref: " << beta_sq_ref << std::endl;
      Number2 beta_ref = 0;
      if (beta_sq_ref > 0) {
          beta_ref = sqrt(beta_sq_ref);
      }

      TVector3 zBoost(0.,0.,-vec_4_ttbar.Pz()/vec_4_ttbar.E()); vec_4_ttbar.Boost(zBoost);
      TVector3 transverseBoost(-vec_4_ttbar.Px()/vec_4_ttbar.E(),-vec_4_ttbar.Py()/vec_4_ttbar.E(),0.); vec_4_ttbar.Boost(transverseBoost);
      vec_4_t.Boost(zBoost); vec_4_t.Boost(transverseBoost);
      vec_4_lbar.Boost(zBoost); vec_4_lbar.Boost(transverseBoost);
      vec_4_tbar.Boost(zBoost); vec_4_tbar.Boost(transverseBoost);
      vec_4_l.Boost(zBoost); vec_4_l.Boost(transverseBoost);
      
      Number2 z = vec_4_t.CosTheta();
      Number2 z_sq = z*z;

      Number2 minusPhi1Plane = -vec_4_t.Phi();
      TVector3 zVect(0.,0.,1.);
      vec_4_t.Rotate(minusPhi1Plane,zVect);
      vec_4_lbar.Rotate(minusPhi1Plane,zVect);
      vec_4_tbar.Rotate(minusPhi1Plane,zVect);
      vec_4_l.Rotate(minusPhi1Plane,zVect);
      Number2 t1Yangle = -atan2(vec_4_t.Px(),vec_4_t.Pz());
      TVector3 yVect(0.,1.,0.);
      vec_4_t.Rotate(t1Yangle,yVect);
      vec_4_lbar.Rotate(t1Yangle,yVect);
      vec_4_tbar.Rotate(t1Yangle,yVect);
      vec_4_l.Rotate(t1Yangle,yVect);

      TVector3 t1Boost(-vec_4_t.Px()/vec_4_t.E(),-vec_4_t.Py()/vec_4_t.E(),-vec_4_t.Pz()/vec_4_t.E()); 
      vec_4_lbar.Boost(t1Boost);
      TLorentzVector P4gen_d1_trest = vec_4_lbar;
      TVector3 t2Boost(-vec_4_tbar.Px()/vec_4_tbar.E(),-vec_4_tbar.Py()/vec_4_tbar.E(),-vec_4_tbar.Pz()/vec_4_tbar.E()); 
      vec_4_l.Boost(t2Boost);
      TLorentzVector vec_4_l_trest = vec_4_l;


      Number2 common_factor_for_M2_gg_QCD = pi_sq/12*(7+9*beta_sq*z_sq)/pow(1-beta_sq*z_sq,2);
      Number2 common_factor_for_M2_gg_ss_nointerf = pi_sq/12*(5+9*beta_sq*z_sq)/pow(1-beta_sq*z_sq,2);
      Number2 common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta_sq*z_sq);


      using namespace std::complex_literals;
      //const std::complex<Number2> i_cmplx = 1i;
      //const std::complex<Number2>{0.0, static_cast<double>(d)};
      const std::complex i_cmplx = std::complex<Number2>(0., 1.L);
      //const std::complex<Number2> i_cmplx = static_cast<Number2>(1);


      std::complex<Number2> NB;
      std::complex<Number2> denomH;

      std::complex<Number2> common_factor_for_M_H = 0.;
      Number2 mh_gh_partial = 0.;
      Number2 mh_gh = 0.;
      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
            if (httInput.MH>2*m_t_ref_sq) {
                  Number2 beta3_tdecay = pow(1-4*m_t_ref_sq/m_H_ref_sq,1.5);
                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
                        mh_gh_partial += 3*G_F*m_t_ref_sq*s/(4*pi*sqrt_2)*beta3_tdecay*ytop_sq;
                  } else {
                        mh_gh_partial += 3*G_F*m_t_ref_sq*m_H_ref_sq/(4*pi*sqrt_2)*beta3_tdecay*ytop_sq;
                  }
            }

            if (httInput.WIDTH_OPTION==1) {
                  mh_gh = httInput.GH*s/httInput.MH;
            } else if (httInput.WIDTH_OPTION==0) {
                  mh_gh = httInput.GH*httInput.MH;
            } else {
                  mh_gh = mh_gh_partial;
            }

            Number2 NB_real, NB_imag;
            Number2 auxh = 1.5 * (1-beta_sq_ref) * ytop_sq;
            if (beta_sq_ref>0) {
                  Number2 betaRef = sqrt(beta_sq_ref);
                  NB_real = auxh * (1 - beta_sq_ref / 4 * (pow( log((1+betaRef) / (1-betaRef)), 2) - pi_sq ));
                  NB_imag = auxh * (beta_sq_ref / 2 * pi * (log((1+betaRef)/(1-betaRef))));
            } else {
                  NB_real = auxh * (1 + beta_sq_ref * pow(asin(sqrt_s/2/m_t_ref_sq), 2));
                  NB_imag = 0.;
            }
            NB = NB_real + NB_imag * i_cmplx;
            denomH = s - m_H_ref_sq + mh_gh * i_cmplx;

            common_factor_for_M_H = NB / denomH;
            common_factor_for_M_H *= pow(s,1.5)*m_t_ref_sq*G_F/6.*sqrt(3.)/4./pi;
      }

      std::complex<Number2> PB;
      std::complex<Number2> denomA;

      std::complex<Number2> common_factor_for_M_A = 0.;
      Number2 ma_ga_partial = 0.;
      Number2 ma_ga = 0.;
      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
            if (httInput.MA>2*m_t_ref_sq) {
                  Number2 beta_tdecay = pow(1-4*m_t_ref_sq/m_A_ref_sq,0.5);
                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
                        ma_ga_partial += 3*G_F*m_t_ref_sq*s/(4*pi*sqrt_2)*beta_tdecay*ytop_sq;
                  } else {
                        ma_ga_partial += 3*G_F*m_t_ref_sq*m_A_ref_sq/(4*pi*sqrt_2)*beta_tdecay*ytop_sq;
                  }
            }
            if (httInput.WIDTH_OPTION==1) {
                        ma_ga = httInput.GA*s/httInput.MA;
            } else if (httInput.WIDTH_OPTION==0) {
                        ma_ga = httInput.GA*httInput.MA;
            } else {
                        ma_ga = ma_ga_partial;
            }

            Number2 PB_real, PB_imag;
            Number2 auxa = -1.5*(1-beta_sq_ref)/4.*ytop_sq;
            if (beta_sq_ref>0) {
                  Number2 betaRef = sqrt(beta_sq_ref);
                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi_sq);
                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  PB_real = -4*auxa*pow(asin(sqrt_s/2/m_t_ref_sq),2);
                  PB_imag = 0.;
            }
            PB = PB_real + PB_imag * i_cmplx;
            denomA = s - m_A_ref_sq + ma_ga * i_cmplx;

            common_factor_for_M_A = PB / denomA;
            common_factor_for_M_A *= pow(s,1.5)*m_t_ref_sq*G_F/6.*sqrt(3.)/4./pi;
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
                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/s*httInput.MH);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed width: GH(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", mh_gh/s*httInput.MH);
                  }
            } else if (httInput.HIGGS_OPTION==1) {
                  printf("\tOnly PSEUDO-SCALAR CASE, MA=%.3f GeV\n", httInput.MA);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed width: GA=%.3f GeV\n", httInput.GA);
                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning width: GA(MA)=%.3f GeV\n", httInput.GA);
                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/s*httInput.MA);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed width: GA(top+bottom)=%.3f GeV\n", ma_ga/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", ma_ga/s*httInput.MA);
                  }
            } else {
                  printf("\tSCALAR + PSEUDO-SCALAR CASE, MH=%.3f GeV, MA=%.3f GeV\n", httInput.MH, httInput.MA);
                  if (httInput.WIDTH_OPTION==0) {
                        printf("\tNon-running fixed widths: GH=%.3f GeV, GA=%.3f GeV\n", httInput.GH, httInput.GA);
                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH, ma_ga_partial/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==1) {
                        printf("\tRunning widths: GH(MH)=%.3f GeV, GA(MA)=%.3f GeV\n", httInput.GH, httInput.GA);
                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/s*httInput.MH, ma_ga_partial/s*httInput.MA);
                  } else if (httInput.WIDTH_OPTION==2) {
                        printf("\tNon-running fixed widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH, ma_ga/httInput.MA);
                  } else if (httInput.WIDTH_OPTION==3) {
                        printf("\tRunning widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/s*httInput.MH, ma_ga/s*httInput.MA);
                  }
            }
            printf("\n");
      }

      Number2 c1 = P4gen_d1_trest.CosTheta();
      //std::cout << "c1: " << c1 << std::endl;
      Number2 c2 = vec_4_l_trest.CosTheta();
      Number2 s1 = sqrt(1-c1*c1);
      Number2 s2 = sqrt(1-c2*c2);
      Number2 phi1 = P4gen_d1_trest.Phi();
      Number2 phi2 = vec_4_l_trest.Phi();
     
      // os factor deleted one sqrt in last line 
      Number2 os_factor = (1 + z_sq + (1 - z_sq) * (1 - beta_sq_ref)) * (1. - c1*c2)
            + 2 * (1 - beta_sq_ref) * (1 - z_sq) * c1 * c2
            - beta_sq_ref * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
            - 2 * (1 - beta_sq_ref) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
            + 2 * sqrt(1 - beta_sq_ref) * z * (1-z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));


      // old os factor
      //Number2 os_factor = (1+z_sq+(1-z_sq)*(1-beta_sq_ref))*(1.-c1*c2) + 2*(1-beta_sq_ref)*(1-z_sq)*c1*c2
      //      - beta_sq_ref*(1-z_sq)*s1*s2*cos(phi1+phi2)
      //      -2*(1-beta_sq_ref)*(1-z_sq)*s1*s2*cos(phi1)*cos(phi2)
      //      +2*sqrt(1-beta_sq_ref)*z*sqrt(1-z_sq)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
      //std::cout << "os_factor: " << os_factor << std::endl;
      Number2 ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
      //std::cout << "ss_factor_1: " << ss_factor_1 << std::endl;
      Number2 ss_factor_beta = beta_sq * (1+c1*c2 - s1*s2*cos(phi1-phi2));
      
      Number2 additional_factor_for_nointerf = 8 * pow(m_t, 2) * pi * alpha_s * beta / pow(s, 2);
      
      //std::cout << "ss_factor_beta" << ss_factor_beta << std::endl;
      //Number2 M2_QCD = common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor
      //              + common_factor_for_M2_gg_QCD*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta);
      
      // just to try it out, new M for QCD, looking at the other changes I made 
      // first part like the qcd part in BSM term and second part as it is, but (1-beta^2) becomes the same factor as in the no interference part
      Number2 M2_QCD = common_factor_for_M2_gg_QCD * pow(alpha_s, 2) * pow(beta, 3) / (8 * pi * s) * (1 - z_sq) * os_factor
                    + common_factor_for_M2_gg_QCD * additional_factor_for_no_interf * (ss_factor_1 + ss_factor_beta);

      Number2 beta_z_factor = 1 - beta_sq_ref * z_sq;
      Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
      Number2 factor_interference = 4 / 3 * pow(m_t, 2) * pow(pi, 3) * pow(alpha_s, 2) / (pow(s, 2) * beta_z_factor_sq);
      Number2 factor_in_norm_interference = sqrt_2 * pow(s, 2) * G_F * beta_z_factor / (16 * pi_sq);
     
      Number2 factor_interference_pseudo = 0;
      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
          factor_interference_pseudo = pow(std::norm(std::complex<Number2>(1., 0.) - ( factor_in_norm_interference * PB / denomA )), 2);
      }
      std::cout << "N = " << NB << std::endl;
      std::cout << "denom H: " << denomH << std::endl;
      Number2 factor_interference_scalar = 0; 
      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
          factor_interference_scalar = pow(std::norm(std::complex<Number2>(1., 0.) - ( factor_in_norm_interference * NB / denomH )), 2);
      }
      // terms to add for the BSM cross section
      Number2 qcd_term = pow(alpha_s, 2) * pow(beta, 3) / (8 * s * pi) * common_factor_for_M2_gg_QCD * (1 - z_sq) * os_factor;
      Number2 no_interference_term = common_factor_for_M2_gg_ss_nointerf * additional_factor_for_nointerf * (ss_factor_1 + ss_factor_beta);

      // if beta_sq_ref < 0, calculate the old way, because there's no decription of what to do in that case in the paper
      Number2 interference_term_pseudo;
      if (beta_sq_ref>0) {
            interference_term_pseudo = factor_interference * beta_ref * factor_interference_pseudo * ss_factor_1;
      }
      else {
            std::cout << "oh no! beta_sq_ref < 0, I'll use the old definition." << std::endl;
            interference_term_pseudo = + std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_A) * ss_factor_1;
      }

      Number2 interference_term_scalar;
      if (beta_sq_ref>0) {
            interference_term_scalar = factor_interference * pow(beta_ref, 3) * factor_interference_scalar * ss_factor_beta;
      }
      else {
            std::cout << "oh no! beta_sq_ref < 0, I'll use the old definition." << std::endl;
            interference_term_scalar = std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_H) * ss_factor_beta;
      }

      std::cout << "qcd term: " << qcd_term << std::endl;
      std::cout << "no interference term: " << no_interference_term << std::endl;
      std::cout << "interference term pseudo scalar: " << interference_term_pseudo << std::endl;
      std::cout << "interference term scalar: " << interference_term_scalar << std::endl;  

      // BSM cross section
      Number2 M2 = qcd_term + no_interference_term + interference_term_pseudo + interference_term_scalar;
      
      std::cout << "M2_gg_QCD factor: " << common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor << std::endl;
      std::cout << "M_gg_ss_nointerf factor: " << common_factor_for_M2_gg_ss_nointerf*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta) << std::endl;
      std::cout << "M_gg_ss_interf M_A factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1 << std::endl;
      std::cout << "M_gg_ss_interf M_H factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta << std::endl;

      // Final weight
      Number2 weight = M2 / M2_QCD;
      if (isnan(weight)) {
            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
            return 0.;
      }

     // return weight;
     return weight - 1;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool is_gg_initial_state(const int& pdgId1, const int& pdgId2, const double& pz_hard_radiation) {
/*
 * pdgId1 is the the PDG id of first parton, along the +Z direction
 * pdgId2 is the the PDG id of other parton, along the -Z direction
 * pz_hard_radiation is the pz of the system made of the hard radiated partons recoiling 
 *    against the ttbar system, measured in the LAB system
*/
      
      // Ultra-simple assignment ito gg or qqbar for qg initial states
      if (pdgId1==21 && pdgId2==21) {
            // gg
            return true;
      } else if (abs(pdgId1)<6 && pdgId2==-pdgId1) {
            // qqbar
            return false;
      } else if (pdgId1==21 || pdgId2==21) {
            // qg
            if ((pdgId1!=21 && pz_hard_radiation>0) || (pdgId2!=21 && pz_hard_radiation<0)) {
                  // radiated quark more collinear with initial quark: assume it is gg
                  return true;
            } else {
                  // radiated quark more collinear with initial gluon: assume it is qqbar
                  return false;
            }
      } 

      // If we got here then it is a rare initial state with multiple hard radiated partons: assume gg
      return true;
}

template <typename Number1, typename Number2>
Number2 compute_weight(Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
                       Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
                       Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
                       Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m)
                       //Number ini1_pt, Number ini1_eta, Number ini1_phi, Number ini1_m,
                       //Number ini2_pt, Number ini2_eta, Number ini2_phi, Number ini2_m,
                       //int pdg1, int pdg2)
{
  // what's that?
  //double pz_hard_radiation = 0.;

 // bool is_gg = is_gg_initial_state(pdg1, pdg2, pz_hard_radiation); 
	
  HTT_Input<Number2> httInput;
  httInput.HIGGS_OPTION = 1;
  httInput.WIDTH_OPTION = 0;
  httInput.YTOP = 1.0;
  httInput.MH = 1100.;
  //httInput.MH = 2100.;
  httInput.GH = 55.176;
  //httInput.GH = 100;
  httInput.MA = 400.;
  //httInput.MA = 1300.;
  httInput.GA = 20;
  //httInput.GA = 100;

  // std::cout << "pTop_pt: " << pTop_pt << ", pTop_eta: " << pTop_eta << ", pTop_phi: " << pTop_phi << ", pTop_m: " << pTop_m << std::endl;
  // std::cout << "aTop_pt: " << aTop_pt << ", aTop_eta: " << aTop_eta << ", aTop_phi: " << aTop_phi << ", aTop_m: " << aTop_m << std::endl;
  // std::cout << "pLep_pt: " << pLep_pt << ", pLep_eta: " << pLep_eta << ", pLep_phi: " << pLep_phi << ", pLep_m: " << pLep_m << std::endl;
  // std::cout << "aLep_pt: " << aLep_pt << ", aLep_eta: " << aLep_eta << ", aLep_phi: " << aLep_phi << ", aLep_m: " << aLep_m << std::endl;

  // Kinematics in this event
  httInput.P4gen_t[0].SetPtEtaPhiM(pTop_pt, pTop_eta, pTop_phi, pTop_m);
  httInput.P4gen_t[1].SetPtEtaPhiM(aTop_pt, aTop_eta, aTop_phi, aTop_m);
  httInput.P4gen_d[1].SetPtEtaPhiM(pLep_pt, pLep_eta, pLep_phi, pLep_m);
  httInput.P4gen_d[0].SetPtEtaPhiM(aLep_pt, aLep_eta, aLep_phi, aLep_m);
 
  //Number2 event_weight = 1.;
  Number2 event_weight = weight_ggHtt(httInput);

  TLorentzVector tops = httInput.P4gen_t[0] + httInput.P4gen_t[1];
  
  std::cout << "ttbar mass: " << tops.M() << std::endl;
  //std::cout << "antitop mass: " << aTop_m << std::endl;
  //printf("Weight=%.3e\n",event_weight);
  return event_weight;

}

template <typename Number1, typename Number2>
auto calculate_weight_juan_paper()
{

  return [] (Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
             Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
             Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
             Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m)
           //            Number ini1_pt, Number ini1_eta, Number ini1_phi, Number ini1_m,
           //            Number ini2_pt, Number ini2_eta, Number ini2_phi, Number ini2_m)
  {
    return compute_weight<Number1, Number2>(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                          aTop_pt, aTop_eta, aTop_phi, aTop_m,
                          pLep_pt, pLep_eta, pLep_phi, pLep_m,
                          aLep_pt, aLep_eta, aLep_phi, aLep_m);
   //                       ini1_pt, ini1_eta, ini1_phi, ini1_m,
   //                       ini2_pt, ini2_eta, ini2_phi, ini2_m,
   //                       pdg1, pdg2);

  };
}



