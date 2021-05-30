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



TLorentzVector f_zmf_tt (TLorentzVector p4_particle, TLorentzVector p4_new_rest_frame) {
  p4_particle.Boost( -1. * p4_new_rest_frame.BoostVector() );
  return p4_particle;
}

enum class calc_weight_version{
  juan_paper, juan_code, juan_paper_different_M2, undefined
};

enum class higgs_type_t{
  scalar, pseudo_scalar, undefined
};

enum class res_int_t{
  resonance, interference, both, undefined
};


template <typename Number2>
struct constants{
      static constexpr Number2 m_t_ref = 172.5;
      static constexpr Number2 m_t_ref_sq = m_t_ref * m_t_ref;
      //Fermi-Konstante
      static constexpr Number2 G_F = 1.16637876e-5;
      static constexpr Number2 alpha_s = 0.12;
      //using namespace std::complex_literals;
      static constexpr std::complex i_cmplx = std::complex<Number2>(0., 1.L);
      // virtual method, so that constants is a pure vitual class, that is not instantiable 
      virtual ~constants() = 0;
};

template <typename Number2>
struct event_t{
    TLorentzVector vec_4_t;
    TLorentzVector vec_4_tbar;
    TLorentzVector vec_4_l;
    TLorentzVector vec_4_lbar;
    higgs_type_t higgs_type;
    Number2 higgs_mass;
    Number2 higgs_width;
    res_int_t res_int;
    calc_weight_version version;

    TLorentzVector vec_4_ttbar;
    Number2 sqrt_s;
    Number2 s;
    Number2 m_t;
    Number2 m_tbar;
    Number2 beta;
    Number2 beta_sq;
    Number2 beta_ref;
    Number2 beta_ref_sq;
    Number2 z;
    Number2 z_sq;
    TLorentzVector vec_4_l_trest;
    TLorentzVector vec_4_lbar_trest;
    Number2 ss_factor_1;
    Number2 ss_factor_beta;
    Number2 os_factor_juan_code;
    Number2 os_factor_juan_paper;
    Number2 y_top = 1;
    
    event_t(TLorentzVector vec_4_t, TLorentzVector vec_4_tbar, TLorentzVector vec_4_l, TLorentzVector vec_4_lbar, const higgs_type_t higgs_type, const Number2 mass, const Number2 width, const res_int_t res_int, const calc_weight_version version) 
    : vec_4_t(vec_4_t), vec_4_tbar(vec_4_tbar), vec_4_l(vec_4_l), vec_4_lbar(vec_4_lbar), higgs_type(higgs_type), higgs_mass(mass), higgs_width(width), res_int(res_int), version(version)

    {    
        TLorentzVector vec_4_ttbar = vec_4_t + vec_4_tbar;
        this->vec_4_ttbar = vec_4_ttbar;
        Number2 sqrt_s = vec_4_ttbar.M();
        this->sqrt_s = sqrt_s;
        Number2 s = sqrt_s * sqrt_s;
        std::cout << "s = " << s << std::endl;
        this->s = s;
        Number2 m_t = vec_4_t.M();
        std::cout << "m_t = " << m_t << std::endl;
        //Number2 m_t = 172.5;
        this->m_t = m_t;
        Number2 m_tbar = vec_4_tbar.M();
        std::cout << "m_tbar = " << m_tbar << std::endl;
        //Number2 m_tbar = 172.5;
        this->m_tbar = m_tbar;
        Number2 beta_sq = (1 - (m_t - m_tbar) * (m_t - m_tbar) / s) * (1 - (m_t + m_tbar) * (m_t + m_tbar) / s);
        std::cout << "beta_sq = " << beta_sq << std::endl;
        this->beta_sq = beta_sq;
        Number2 beta = 0;
        if (beta_sq >= 0){
            beta = sqrt(beta_sq);
        }
        else {
            std::cout << "Beta is smaller than 0!!! using 0";
        }
        this->beta = beta;
        //Number2 beta = sqrt(beta_sq);

        //Number2 beta3 = beta*beta_sq;
        //std::cout << "top mass: " << m_t << std::endl;
        Number2 mtrefOverE = 2 * constants<Number2>::m_t_ref / sqrt_s;
        std::cout << "mtrefOverE = " << mtrefOverE << std::endl;
        //std::cout << "square root of s: " << sqrt_s << std::endl;
        //std::cout << "mt over E: " << mtrefOverE << std::endl; 
        Number2 beta_ref_sq = 1 - mtrefOverE * mtrefOverE;
        
        //this->beta_ref_sq = beta_ref_sq;
        this->beta_ref_sq = beta_ref_sq;
        //std::cout << "beta squared ref: " << beta_ref_sq << std::endl;
        Number2 beta_ref = 0;
        if (beta_ref_sq > 0) {
            beta_ref = sqrt(beta_ref_sq);
        }
        //this->beta_ref = beta_ref;
        this->beta_ref = beta_ref;
        std::cout << "beta_ref = " << beta_ref << std::endl;

        static const TVector3 zBase(0., 0., 1.);
        TLorentzVector vec_4_t_hel = f_zmf_tt(vec_4_t, vec_4_ttbar);
        TLorentzVector vec_4_tbar_hel = f_zmf_tt(vec_4_tbar, vec_4_ttbar);
        TLorentzVector vec_4_l_hel = f_zmf_tt(vec_4_l, vec_4_ttbar);
        TLorentzVector vec_4_lbar_hel = f_zmf_tt(vec_4_lbar, vec_4_ttbar);
        //Number2 z = vec_4_t_hel.Vect().Unit().Dot( vec_4_ttbar.Vect().Unit() ); 
        Number2 z = vec_4_t_hel.Vect().Unit().Dot(zBase);

        //TVector3 zBoost(0.,0.,-vec_4_ttbar.Pz()/vec_4_ttbar.E()); vec_4_ttbar.Boost(zBoost);
        //TVector3 transverseBoost(-vec_4_ttbar.Px()/vec_4_ttbar.E(),-vec_4_ttbar.Py()/vec_4_ttbar.E(),0.); vec_4_ttbar.Boost(transverseBoost);
        //vec_4_t.Boost(zBoost); vec_4_t.Boost(transverseBoost);
        //vec_4_lbar.Boost(zBoost); vec_4_lbar.Boost(transverseBoost);
        //vec_4_tbar.Boost(zBoost); vec_4_tbar.Boost(transverseBoost);
        //vec_4_l.Boost(zBoost); vec_4_l.Boost(transverseBoost);
        
        //Number2 z = vec_4_t.CosTheta();
        this->z = z;
        std::cout << "z = " << z << std::endl;
        Number2 z_sq = z * z;
        this->z_sq = z_sq; 

        TLorentzVector vec_4_l_tbar_rest = f_zmf_tt(vec_4_l_hel, vec_4_tbar_hel);
        TLorentzVector vec_4_lbar_t_rest = f_zmf_tt(vec_4_lbar_hel, vec_4_t_hel);

        const Number2 spTP = std::sqrt(1. - z_sq);
        const Number2 sY = (z >= 0.) ? 1. : -1.;
        const TVector3 kBase = vec_4_t_hel.Vect().Unit();
        const TVector3 rBase = ( (sY / spTP) * (zBase - (z * kBase)) ).Unit();
        const TVector3 nBase = ( (sY / spTP) * zBase.Cross(kBase) ).Unit();
        Number2 b1k = vec_4_lbar_t_rest.Vect().Unit().Dot( kBase );
        Number2 b2k_like = vec_4_l_tbar_rest.Vect().Unit().Dot( kBase );
        Number2 b1r = vec_4_lbar_t_rest.Vect().Unit().Dot( rBase );
        Number2 b2r_like = vec_4_l_tbar_rest.Vect().Unit().Dot( rBase );
        Number2 b1n = vec_4_lbar_t_rest.Vect().Unit().Dot( nBase );
        Number2 b2n_like = vec_4_l_tbar_rest.Vect().Unit().Dot( nBase );

        Number2 c1 = b1k; 
        Number2 c2 = b2k_like;
        Number2 phi1 = atan2(b1n, b1r);
        Number2 phi2 = atan2(b2n_like, b2r_like); 

//        Number2 minusPhi1Plane = -vec_4_t_hel.Phi();
//        TVector3 zVect(0.,0.,1.);
//        vec_4_t_hel.Rotate(minusPhi1Plane,zVect);
//        vec_4_lbar_hel.Rotate(minusPhi1Plane,zVect);
//        vec_4_tbar_hel.Rotate(minusPhi1Plane,zVect);
//        vec_4_l_hel.Rotate(minusPhi1Plane,zVect);
//        Number2 t1Yangle = -atan2(vec_4_t_hel.Px(),vec_4_t_hel.Pz());
//        TVector3 yVect(0.,1.,0.);
//        vec_4_t_hel.Rotate(t1Yangle,yVect);
//        vec_4_lbar_hel.Rotate(t1Yangle,yVect);
//        vec_4_tbar_hel.Rotate(t1Yangle,yVect);
//        vec_4_l_hel.Rotate(t1Yangle,yVect);
//
//        TVector3 t1Boost(-vec_4_t_hel.Px()/vec_4_t_hel.E(),-vec_4_t_hel.Py()/vec_4_t_hel.E(),-vec_4_t_hel.Pz()/vec_4_t_hel.E()); 
//        vec_4_lbar_hel.Boost(t1Boost);
//        TLorentzVector vec_4_lbar_trest = vec_4_lbar_hel;
//        TVector3 t2Boost(-vec_4_tbar_hel.Px()/vec_4_tbar_hel.E(),-vec_4_tbar_hel.Py()/vec_4_tbar_hel.E(),-vec_4_tbar_hel.Pz()/vec_4_tbar_hel.E()); 
//        vec_4_l_hel.Boost(t2Boost);
//        TLorentzVector vec_4_l_trest = vec_4_l_hel;
//        this->vec_4_l_trest = vec_4_l_trest;
//        this->vec_4_lbar_trest = vec_4_lbar_trest;    
//        
//
//        Number2 c1 = vec_4_lbar_trest.CosTheta();
//        //std::cout << "c1: " << c1 << std::endl;
//        Number2 c2 = vec_4_l_trest.CosTheta();
        Number2 s1 = sqrt(1-c1*c1);
        std::cout << "s1 = " << s1 << std::endl;
        Number2 s2 = sqrt(1-c2*c2);
        std::cout << "s2 = " << s2 << std::endl;
//        Number2 phi1 = vec_4_lbar_trest.Phi();
//        Number2 phi2 = vec_4_l_trest.Phi();
     
        // os factor deleted one sqrt in last line 
        Number2 os_factor_juan_paper = (1 + z_sq + (1 - z_sq) * (1 - beta_ref_sq)) * (1. - c1*c2)
              + 2 * (1 - beta_ref_sq) * (1 - z_sq) * c1 * c2
              - beta_ref_sq * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
              - 2 * (1 - beta_ref_sq) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
              + 4 * sqrt(1 - beta_ref_sq) * z * (1 - z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
        this->os_factor_juan_paper = os_factor_juan_paper;
        std::cout << "os_factor_juan_paper = " << os_factor_juan_paper << std::endl;

        // old os factor
        Number2 os_factor_juan_code = (1 + z_sq + (1 - z_sq) * (1 - beta_ref_sq)) * (1. - c1 * c2) 
              + 2 * (1 - beta_ref_sq) * (1 - z_sq) * c1 * c2
              - beta_ref_sq * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
              - 2 * (1 - beta_ref_sq) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
              + 2 * sqrt(1 - beta_ref_sq) * z * sqrt(1 - z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
        this->os_factor_juan_code = os_factor_juan_code;
        std::cout << "os_factor_juan_code = " << os_factor_juan_code << std::endl;

        //std::cout << "os_factor: " << os_factor << std::endl;
        Number2 ss_factor_1 = 1 + c1 * c2 + s1 * s2 * cos(phi1 - phi2);
        std::cout << "ss_factor_1 = " << ss_factor_1 << std::endl;
        this->ss_factor_1 = ss_factor_1;
        //std::cout << "ss_factor_1: " << ss_factor_1 << std::endl;
        Number2 ss_factor_beta = beta_sq * (1 + c1 * c2 - s1 * s2 * cos(phi1 - phi2));
        this->ss_factor_beta = ss_factor_beta;
        std::cout << "ss_factor_beta = " << ss_factor_beta << std::endl;
    }

};


// Reweighting gg-> H/A -> ttbar events
// Input is an HTT_Input structure; output is the event weight
// template <typename Number2>
// Number2 weight_ggHtt(event_t<Number2> event);

// Simple code to decide whether the initial state should be gg or not
// pdgId1 is the the PDG id of first parton, along the +Z direction
// pdgId2 is the the PDG id of other parton, along the -Z direction
// pz_hard_radiation is the pz of the system made of the hard radiated partons recoiling 
// template <typename Number2>
// bool is_gg_initial_state(const int& pdgId1, const int& pdgId2, const Number2& pz_hard_radiation);

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

      Number2 MH2_REF = pow(httInput.MH,2);
      Number2 MA2_REF = pow(httInput.MA,2);
      Number2 ytop2 = pow(httInput.YTOP,2);

      TLorentzVector P4gen_t1 = httInput.P4gen_t[0]; 
      TLorentzVector P4gen_t2 = httInput.P4gen_t[1]; 
      TLorentzVector P4gen_d1 = httInput.P4gen_d[0]; 
      TLorentzVector P4gen_d2 = httInput.P4gen_d[1]; 

      const Number2 MT_REF = 172.5;
      const Number2 MT2_REF = MT_REF*MT_REF;
      const Number2 GF = 1.16637876e-5;
      const Number2 pi = TMath::Pi();
      const Number2 pi2 = pi*pi;
      const Number2 sqrt2 = sqrt(2.);

      TLorentzVector P4gen_tt = P4gen_t1+P4gen_t2;
      Number2 sqrts = P4gen_tt.M();
      Number2 sqrts2 = sqrts*sqrts;
      
      Number2 mt1 = httInput.P4gen_t[0].M();
      std::cout << "printf(\"m_t = " << mt1 << "\\n\");" << std::endl;
      Number2 mt2 = httInput.P4gen_t[1].M();
      std::cout << "printf(\"m_tbar = " << mt2 << "\\n\");" << std::endl;
      
      Number2 beta2 = (1-(mt1-mt2)*(mt1-mt2)/sqrts2)*(1-(mt1+mt2)*(mt1+mt2)/sqrts2);
      //Number2 beta = sqrt(beta2);
      //Number2 beta3 = beta*beta2;
      Number2 mtrefOverE = 2*MT_REF/sqrts;
      Number2 beta2Ref = 1-mtrefOverE*mtrefOverE;

      TVector3 zBoost(0.,0.,-P4gen_tt.Pz()/P4gen_tt.E()); P4gen_tt.Boost(zBoost);
      TVector3 transverseBoost(-P4gen_tt.Px()/P4gen_tt.E(),-P4gen_tt.Py()/P4gen_tt.E(),0.); P4gen_tt.Boost(transverseBoost);
      P4gen_t1.Boost(zBoost); P4gen_t1.Boost(transverseBoost);
      P4gen_d1.Boost(zBoost); P4gen_d1.Boost(transverseBoost);
      P4gen_t2.Boost(zBoost); P4gen_t2.Boost(transverseBoost);
      P4gen_d2.Boost(zBoost); P4gen_d2.Boost(transverseBoost);
      Number2 z1 = P4gen_t1.CosTheta();
      std::cout << "printf(\"z = " << z1 << "\\n\");" << std::endl;
      Number2 z2 = z1*z1;

      Number2 minusPhi1Plane = -P4gen_t1.Phi();
      TVector3 zVect(0.,0.,1.);
      P4gen_t1.Rotate(minusPhi1Plane,zVect);
      P4gen_d1.Rotate(minusPhi1Plane,zVect);
      P4gen_t2.Rotate(minusPhi1Plane,zVect);
      P4gen_d2.Rotate(minusPhi1Plane,zVect);
      Number2 t1Yangle = -atan2(P4gen_t1.Px(),P4gen_t1.Pz());
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

      Number2 common_factor_for_M2_gg_QCD = pi2/12*(7+9*beta2*z2)/pow(1-beta2*z2,2);
      std::cout << "printf(\"common_factor_for_M2_gg_QCD = " << common_factor_for_M2_gg_QCD << "\\n\");" << std::endl;
      Number2 common_factor_for_M2_gg_ss_nointerf = pi2/12*(5+9*beta2*z2)/pow(1-beta2*z2,2);
      std::cout << "printf(\"common_factor_for_M2_gg_ss_nointerf = " << common_factor_for_M2_gg_ss_nointerf << "\\n\");" << std::endl;
      Number2 common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta2*z2);
      std::cout << "printf(\"common_factor_for_M_gg_ss_interf = " << common_factor_for_M_gg_ss_interf << "\\n\");" << std::endl;

      using namespace std::literals::complex_literals;
      using namespace std::literals;
      const std::complex<Number2> i_cmplx = std::complex<Number2>(0., 1.L);

      std::complex<Number2> common_factor_for_M_H = 0.;
      Number2 mh_gh_partial = 0.;
      Number2 mh_gh = 0.;
      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
            if (httInput.MH>2*MT_REF) {
                  Number2 beta3_tdecay = pow(1-4*MT2_REF/MH2_REF,1.5);
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

            Number2 NB_real, NB_imag;
            Number2 auxh = 1.5*(1-beta2Ref)*ytop2;
            if (beta2Ref>0) {
                  Number2 betaRef = sqrt(beta2Ref);
                  NB_real = auxh*(1-beta2Ref/4*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2));
                  NB_imag = auxh*(beta2Ref/2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  NB_real = auxh*(1+beta2Ref*pow(asin(sqrts/2/MT_REF),2));
                  NB_imag = 0.;
            }
            std::complex<Number2> NB = NB_real + NB_imag*i_cmplx;
            std::complex<Number2> denomH = sqrts2 - MH2_REF + mh_gh*i_cmplx;

            common_factor_for_M_H = NB / denomH;
            common_factor_for_M_H *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
      }

      std::complex<Number2> common_factor_for_M_A = 0.;
      Number2 ma_ga_partial = 0.;
      Number2 ma_ga = 0.;
      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
            if (httInput.MA>2*MT_REF) {
                  Number2 beta_tdecay = pow(1-4*MT2_REF/MA2_REF,0.5);
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

            Number2 PB_real, PB_imag;
            Number2 auxa = -1.5*(1-beta2Ref)/4.*ytop2;
            if (beta2Ref>0) {
                  Number2 betaRef = sqrt(beta2Ref);
                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2);
                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
            } else {
                  PB_real = -4*auxa*pow(asin(sqrts/2/MT_REF),2);
                  PB_imag = 0.;
            }
            std::complex<Number2> PB = PB_real + PB_imag*i_cmplx;
            std::complex<Number2> denomA = sqrts2 - MA2_REF + ma_ga*i_cmplx;

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

      Number2 c1 = P4gen_d1_trest.CosTheta();
      Number2 c2 = P4gen_d2_trest.CosTheta();
      Number2 s1 = sqrt(1-c1*c1);
      Number2 s2 = sqrt(1-c2*c2);
      Number2 phi1 = P4gen_d1_trest.Phi();
      Number2 phi2 = P4gen_d2_trest.Phi();
      Number2 os_factor = (1+z2+(1-z2)*(1-beta2Ref))*(1.-c1*c2) + 2*(1-beta2Ref)*(1-z2)*c1*c2
            - beta2Ref*(1-z2)*s1*s2*cos(phi1+phi2)
            -2*(1-beta2Ref)*(1-z2)*s1*s2*cos(phi1)*cos(phi2)
            +2*sqrt(1-beta2Ref)*z1*sqrt(1-z2)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
      Number2 ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
      Number2 ss_factor_beta = beta2 * (1+c1*c2 - s1*s2*cos(phi1-phi2));

      Number2 M2_QCD = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
                     + common_factor_for_M2_gg_QCD*(1-beta2Ref)*(ss_factor_1+ss_factor_beta);

      Number2 M2 = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
                + common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
                + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
                + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;
      //Number2 M2 = common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
      //             + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
      //             + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;

      // Final weight
      Number2 weight = M2 / M2_QCD;

      std::cout << "printf(\"M2_QCD = " << M2_QCD << "\\n\");" << std::endl;
      std::cout << "printf(\"M2_BSM = " << M2 << "\\n\");" << std::endl;
      std::cout << "printf(\"weight = " << weight << "\\n\");" << std::endl;

      if (isnan(weight)) {
            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
            return 0.;
      }

      return weight;
}


// ---- N and P ----
template<typename Number2>
std::complex<Number2> calc_N(event_t<Number2> event){

     std::complex<Number2> N;

     if (event.beta_ref_sq>0) {
           Number2 beta_ref = sqrt(event.beta_ref_sq);
           N = (Number2) 3. / (Number2) 8. * (Number2) (1 - event.beta_ref_sq) * (Number2) pow(event.y_top, 2) * ( (Number2) 4. - (std::complex<Number2>) event.beta_ref_sq * (std::complex<Number2>) pow(log(((1 + beta_ref) / (1 - beta_ref))) - constants<Number2>::i_cmplx * (Number2) TMath::Pi(), 2) );
           //N = (Number2) 3. / (Number2) 8. * (Number2) (1 - event.beta_ref_sq) * (Number2) pow(event.y_top, 2) * ( (Number2) 4. - (std::complex<Number2>) event.beta_ref_sq * (std::complex<Number2>) std::norm(log(((1 + beta_ref) / (1 - beta_ref))) - constants<Number2>::i_cmplx * (Number2) TMath::Pi()) );
     } else {
           Number2 auxh = 1.5 * (1 - event.beta_ref_sq) * pow(event.y_top, 2);
           Number2 NB_real = auxh * (1 + event.beta_ref_sq * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2));
           Number2 NB_imag = 0.;
           N = NB_real + NB_imag * constants<Number2>::i_cmplx;
     }
     std::cout << "N=" << N << std::endl;
     return N;
}


template<typename Number2>
std::complex<Number2> calc_P(event_t<Number2> event){

      std::complex<Number2> P;
      if (event.beta_ref_sq>0) {
            Number2 beta_ref = sqrt(event.beta_ref_sq);
            P = - (Number2) pow(event.y_top, 2) * (Number2) 3. / (Number2) 8. * (Number2) (1 - event.beta_ref_sq) * pow( (Number2) log( (1 + beta_ref) / ( 1 - beta_ref)) - constants<Number2>::i_cmplx * (Number2) TMath::Pi(), 2);
            //P = - (Number2) pow(event.y_top, 2) * (Number2) 3. / (Number2) 8. * (Number2) (1 - event.beta_ref_sq) * std::norm( (Number2) log( (1 + beta_ref) / ( 1 - beta_ref)) - constants<Number2>::i_cmplx * (Number2) TMath::Pi());
      } else {
            Number2 auxa = -1.5 * (1 - event.beta_ref_sq) / 4. * pow(event.y_top, 2);
            Number2 PB_real = -4 * auxa * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2);
            Number2 PB_imag = 0.;
            P = PB_real + PB_imag * constants<Number2>::i_cmplx;
      }
      std::cout << "P=" << P << std::endl;
      return P;
}


// ---- methods without decay ----
template<typename Number2>
Number2 calc_resonance_scalar_without_decay(event_t<Number2> event){
     Number2 common_factor = pow(constants<Number2>::G_F, 2) * pow(constants<Number2>::m_t_ref, 2) * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3) );
     
     Number2 mh_gh = event.higgs_width * event.higgs_mass;
     std::complex<Number2> denomH = event.s - (Number2) pow(event.higgs_mass, 2) + mh_gh * constants<Number2>::i_cmplx;
     std::complex<Number2> N = calc_N(event);

     std::cout << "N / denomH = " << N / denomH << std::endl;

     //Number2 resonance = common_factor * pow(event.beta, 3) * std::norm( N / denomH ); 
     Number2 resonance = common_factor * pow(event.beta, 3) * pow(N.real(), 2) * std::norm( denomH ); 
     return resonance;   
}

template<typename Number2>
Number2 calc_resonance_pseudo_without_decay(event_t<Number2> event){
     Number2 common_factor = pow(constants<Number2>::G_F, 2) * pow(constants<Number2>::m_t_ref, 2) * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3) );
     
     Number2 ma_ga = event.higgs_width * event.higgs_mass;
     std::complex<Number2> denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;
     std::complex<Number2> P = calc_P(event);
     
     std::cout << "P / denomA = " << P / denomA << std::endl;

     //Number2 resonance = common_factor * event.beta * std::norm( P / denomA ); 
     Number2 resonance = common_factor * event.beta * pow(P.real(), 2) / std::norm( denomA ); 
     return resonance;   
}

template<typename Number2>
Number2 calc_QCD_without_decay(event_t<Number2> event){

     Number2 factor_1 = TMath::Pi() * event.beta / ( 96 * event.s);
     Number2 factor_2 = (7 + event.beta_sq * event.z_sq) / pow( 1 - event.beta_sq * event.z_sq, 2 ); 
     Number2 factor_3 = event.beta_sq * (1 - pow(event.z , 4)) + 4 * pow(constants<Number2>::m_t_ref, 2) / event.s * ( 1 + event.beta_sq * pow(event.z, 4) + 2 * event.beta_sq - 2 * event.beta_sq * event.z_sq);
     return factor_1 * factor_2 * factor_3;
}

template<typename Number2>
Number2 calc_QCD_without_decay_dicus(event_t<Number2> event){

     Number2 factor_1 = TMath::Pi() * event.beta / ( 12 * event.s);
     Number2 s_sq = pow(event.s, 2);
     Number2 p1_p3 = event.s / 4 * (1. - event.beta * event.z);
     Number2 p2_p3 = event.s / 4 * (1. + event.beta * event.z);
     Number2 factor_2 = s_sq / ( p1_p3 * p2_p3) - 9;
     Number2 factor_3 = pow(p1_p3, 2) / s_sq + pow(p2_p3, 2) / s_sq + constants<Number2>::m_t_ref_sq / event.s - pow(constants<Number2>::m_t_ref, 4) / ( 4 * p1_p3 * p2_p3); 
     return factor_1 * factor_2 * factor_3;
}


template<typename Number2>
Number2 calc_qcd_interference_scalar(event_t<Number2> event){

     if (event.higgs_type != higgs_type_t::scalar)
       return 0;

     //Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
     std::complex<Number2> NB;
     std::complex<Number2> denomH;

     std::complex<Number2> common_factor_for_M_H = 0.;
     Number2 mh_gh = event.higgs_width * event.higgs_mass;

     std::complex<Number2> NB_real, NB_imag;
     Number2 auxh = 1.5 * (1 - event.beta_ref_sq) * pow(event.y_top, 2);
     if (event.beta_ref_sq>0) {
           Number2 beta_ref = sqrt(event.beta_ref_sq);
           NB_real = auxh * (1 - event.beta_ref_sq / 4 * (pow( log((1+beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2) ));
           NB_imag = auxh * (event.beta_ref_sq / 2 * TMath::Pi() * (log((1+beta_ref)/(1-beta_ref))));
     } else {
           NB_real = auxh * (1 + event.beta_ref_sq * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2));
           NB_imag = 0.;
     }
     NB = NB_real + NB_imag * constants<Number2>::i_cmplx;
     denomH = event.s - (Number2) pow(event.higgs_mass, 2) + mh_gh * constants<Number2>::i_cmplx;

     common_factor_for_M_H = NB / denomH;
     common_factor_for_M_H *= pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6.* sqrt(3.) / 4. / TMath::Pi();
     
     Number2 qcd_interference_scalar = 0;
     Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
     Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
     switch (event.version){
        
          case calc_weight_version::juan_code:{
               //qcd_interference_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H)
               //                          * event.ss_factor_beta;
               qcd_interference_scalar = 4. * constants<Number2>::m_t_ref_sq / event.s * pow(TMath::Pi(), 2) / ( 6. * beta_z_factor_sq ) * event.ss_factor_beta; 
               break;}
          case calc_weight_version::juan_paper:{
              Number2 factor_interference = 4. / 48. * constants<Number2>::m_t_ref_sq * TMath::Pi() 
                                            / (pow(event.s, 2) * beta_z_factor_sq);
              //Number2 factor_interference = 4 / 3 * pow(constants<Number2>::m_t_ref, 2) * pow(TMath::Pi(), 3) 
              //                              / (pow(event.s, 2) * beta_z_factor_sq);
              //Number2 factor_in_norm_interference = sqrt(2.) * pow(event.s, 2) * constants<Number2>::G_F * beta_z_factor / (16 * pow(TMath::Pi(), 2));
              // return zero if higgs is not scalar 
              // if (event.higgs_type != higgs_type_t::scalar) {
              //     return 0;
              // } 
              qcd_interference_scalar = factor_interference * pow(event.beta, 3)
                                       // setting N to zero
                                       // * std::norm((Number2) 1. - factor_in_norm_interference * (NB / denomH))
                                       * event.ss_factor_beta;
              break;}
          default:
              break;
     }
      Number2 qcd_interference_scalar_code = 4. * constants<Number2>::m_t_ref_sq / event.s * pow(TMath::Pi(), 2) / ( 6. * beta_z_factor_sq ) * event.ss_factor_beta; 
      std::cout << "qcd interference scalar code: " << qcd_interference_scalar_code << std::endl;      

      Number2 factor_ich_juan = 8. * TMath::Pi() * event.s / pow(event.beta, 3);
      std::cout << "factor ich -> Juan: " << factor_ich_juan << std::endl;
      std::cout << "ich * factor: " << factor_ich_juan * qcd_interference_scalar << std::endl;
     return qcd_interference_scalar;
}

template <typename Number2>
Number2 calc_qcd_interference_pseudo(event_t<Number2> event){

      if (event.higgs_type != higgs_type_t::pseudo_scalar)
        return 0;

      //Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
      std::complex<Number2> PB = calc_P(event);
      std::complex<Number2> denomA;

      std::complex<Number2> common_factor_for_M_A = 0.;
      Number2 ma_ga = event.higgs_mass * event.higgs_width;

      //std::complex<Number2> PB_real, PB_imag;
      //Number2 auxa = -1.5 * (1 - event.beta_ref_sq)/ 4. * pow(event.y_top, 2);
      //if (event.beta_ref_sq>0) {
      //      Number2 beta_ref = sqrt(event.beta_ref_sq);
      //      PB_real = auxa*(pow(log((1 + beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2));
      //      PB_imag = auxa*(-2 * TMath::Pi() * (log((1 + beta_ref) / (1 - beta_ref))));
      //} else {
      //      PB_real = -4 * auxa * pow(asin(event.sqrt_s/ 2 / constants<Number2>::m_t_ref_sq), 2);
      //      PB_imag = 0.;
      //}
      //PB = PB_real + PB_imag * constants<Number2>::i_cmplx;
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;

      common_factor_for_M_A = PB / denomA;
      common_factor_for_M_A *= pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 qcd_interference_pseudo = 0;
 
      Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
      Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
      switch (event.version){
         case calc_weight_version::juan_code:{
             //qcd_interference_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
             //                             * event.ss_factor_1;
             //qcd_interference_pseudo = 4. * constants<Number2>::m_t_ref_sq / event.s * pow(TMath::Pi(), 2) / 36. / beta_z_factor_sq * event.ss_factor_1; 
             qcd_interference_pseudo = 4. * constants<Number2>::m_t_ref_sq / event.s * pow(TMath::Pi(), 2) / ( 6. * beta_z_factor_sq ) * event.ss_factor_1; 
             break;}
         case calc_weight_version::juan_paper:{
             Number2 factor_interference = 4. / 48. * constants<Number2>::m_t_ref_sq * TMath::Pi() / (pow(event.s, 2) * beta_z_factor_sq);
             //Number2 factor_interference = 4 / 3 * pow(constants<Number2>::m_t_ref, 2) * pow(TMath::Pi(), 3)  
             //                              / (pow(event.s, 2) * beta_z_factor_sq);
             //Number2 factor_in_norm_interference = sqrt(2.) * pow(event.s, 2) * constants<Number2>::G_F * beta_z_factor / (16 * pow(TMath::Pi(), 2));
             // return zero if higgs_type is not pseudo
             //if (event.higgs_type != higgs_type_t::pseudo_scalar){
             //    return 0;
             //}         
             qcd_interference_pseudo = factor_interference * event.beta
                                     // setting P to zero??
                                     //* std::norm((Number2) 1. - factor_in_norm_interference * (PB / denomA))
                                      * event.ss_factor_1;
             std::cout << "qcd interference pseudo paper: " <<  qcd_interference_pseudo << std::endl;
             break;}
         default:
             break;
      }

      //Number2 qcd_interference_pseudo_code = 4 * constants<Number2>::m_t_ref_sq / event.s * pow(TMath::Pi(), 2) / 36 / beta_z_factor_sq * event.ss_factor_1; 
      Number2 qcd_interference_pseudo_code = 4. * constants<Number2>::m_t_ref_sq / event.s * pow(TMath::Pi(), 2) / ( 6. * beta_z_factor_sq ) * event.ss_factor_1; 
      std::cout << "qcd interference pseudo code: " << qcd_interference_pseudo_code << std::endl;      

      Number2 factor_ich_juan = 8. * TMath::Pi() * event.s / event.beta;
      std::cout << "factor ich -> Juan: " << factor_ich_juan << std::endl;
      std::cout << "ich * factor: " << factor_ich_juan * qcd_interference_pseudo << std::endl;

      return qcd_interference_pseudo;
}

template <typename Number2>
Number2 calc_qcd_no_interference(event_t<Number2> event){
    Number2 common_factor_for_M2_gg_ss_nointerf = pow(TMath::Pi(), 2) / 12 * (5 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    Number2 additional_factor_for_nointerf = 0;
    switch (event.version){ 
        case calc_weight_version::juan_paper:
             //additional_factor_for_nointerf = 8 * pow(constants<Number2>::m_t_ref, 2) * TMath::Pi() * event.beta / pow(event.s, 2);
             additional_factor_for_nointerf = pow(constants<Number2>::m_t_ref, 2) * event.beta / (2 * TMath::Pi() * pow(event.s, 2));
            break;
        case calc_weight_version::juan_code: 
            additional_factor_for_nointerf = 1 - event.beta_ref_sq;
            break;
        default:
            break;
    }
    
    Number2 additional_factor_for_nointerf_paper = pow(constants<Number2>::m_t_ref, 2) * event.beta / (2 * TMath::Pi() * pow(event.s, 2));
    Number2 additional_factor_for_nointerf_code = 1 - event.beta_ref_sq;
    
    std::cout << "additional factor for no interf paper: " << additional_factor_for_nointerf_paper << std::endl;
    std::cout << "additional factor for no interf code: " << additional_factor_for_nointerf_code << std::endl;

    Number2 factor_ich_juan = 8 * TMath::Pi() * event.s / event.beta;
    std::cout << "factor ich -> Juan: " << factor_ich_juan << std::endl;
    std::cout << "ich * factor: " << factor_ich_juan * additional_factor_for_nointerf_paper << std::endl;

    Number2 no_interference_term = common_factor_for_M2_gg_ss_nointerf * additional_factor_for_nointerf * (event.ss_factor_1 + event.ss_factor_beta);
    return no_interference_term;
}


template <typename Number2>
Number2 calc_qcd_opp_gluon(event_t<Number2> event){
    Number2 common_factor_for_M2_gg_QCD = pow(TMath::Pi(), 2) / 12. * (7. + 9. * event.beta_sq * event.z_sq) / pow(1. - event.beta_sq * event.z_sq, 2);
    std::cout << "common_factor_for_M2_gg_QCD = " << common_factor_for_M2_gg_QCD << std::endl;
    Number2 opp_gluon_qcd_term = 0;
 
    Number2 opp_gluon_qcd_term_juan_paper_pre = pow(event.beta, 3) / (8. * TMath::Pi() * event.s ) * common_factor_for_M2_gg_QCD * (1 - event.z_sq);
    Number2 opp_gluon_qcd_term_juan_paper = opp_gluon_qcd_term_juan_paper_pre * event.os_factor_juan_paper;
    std::cout << "opp gluon qcd paper pre = " << opp_gluon_qcd_term_juan_paper_pre << std::endl;

    Number2 opp_gluon_qcd_term_juan_code_pre = common_factor_for_M2_gg_QCD * event.beta_sq * (1. - event.z_sq); 
    Number2 opp_gluon_qcd_term_juan_code = opp_gluon_qcd_term_juan_code_pre * event.os_factor_juan_code; 
    std::cout << "opp gluon qcd code pre = " << opp_gluon_qcd_term_juan_code_pre << std::endl;
    std::cout << "opp gluon qcd code = " << opp_gluon_qcd_term_juan_code << std::endl;

    Number2 factor_ich_juan = 8. * TMath::Pi() * event.s / event.beta;
    std::cout << "factor ich -> Juan: " << factor_ich_juan << std::endl;
    std::cout << "ich * factor: " << factor_ich_juan * opp_gluon_qcd_term_juan_paper_pre << std::endl;

    if (event.version == calc_weight_version::juan_paper){
       return opp_gluon_qcd_term_juan_paper;
    }
    if (event.version == calc_weight_version::juan_code){
       return opp_gluon_qcd_term_juan_code;
    }
    return opp_gluon_qcd_term;
}

template <typename Number2>
Number2 calc_M_2_qcd(event_t<Number2> event){
    Number2 M2_QCD = 0;
    Number2 common_factor_for_M2_gg_QCD = pow(TMath::Pi(), 2) / 12 * (7 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    //Number2 additional_factor_for_nointerf = 8 * pow(constants<Number2>::m_t_ref, 2) * TMath::Pi() * event.beta_ref / pow(event.s, 2);
    
    //if (event.version == calc_weight_version::juan_code){
        Number2 M2_QCD_juan = common_factor_for_M2_gg_QCD * event.beta_sq * (1 - event.z_sq) * event.os_factor_juan_code 
                 + common_factor_for_M2_gg_QCD * (1 - event.beta_ref_sq) * (event.ss_factor_1 + event.ss_factor_beta);
        std::cout << "M2_QCD code orig = " << M2_QCD_juan << std::endl;
    //}
    //if (event.version ==  calc_weight_version::juan_paper){
        // just to try it out, new M for QCD, looking at the other changes I made 
        // first part like the qcd part in BSM term and second part as it is, but (1-beta^2) becomes the same factor as in the no interference part
      
        std::cout << "beta_ref: " << event.beta_ref << std::endl;
        std::cout << "z_sq: " << event.z_sq << std::endl;
        std::cout << "s: " << event.s << std::endl;

        Number2 opp_gluon = calc_qcd_opp_gluon(event);
        Number2 same_gluon_no_interf = calc_qcd_no_interference(event);
        Number2 same_gluon_interf_pseudo = calc_qcd_interference_pseudo(event);
        Number2 same_gluon_interf_scalar = calc_qcd_interference_scalar(event);
 
        std::cout << "qcd opp gluon: " << opp_gluon << std::endl;
        std::cout << "qcd same gluon no interf: " << same_gluon_no_interf << std::endl;
        std::cout << "qcd same gluon interf: " << same_gluon_interf_pseudo << std::endl;
        std::cout << "qcd same gluon interf: " << same_gluon_interf_scalar << std::endl;

        M2_QCD = opp_gluon + same_gluon_no_interf + same_gluon_interf_pseudo + same_gluon_interf_scalar;
        std::cout << "M2_QCD = " << M2_QCD << std::endl;
        
        //M2_QCD = common_factor_for_M2_gg_QCD * pow(event.beta_ref, 3) / (8 * TMath::Pi() * event.s) * (1 - event.z_sq)
        //         * event.os_factor_juan_paper 
        //         + common_factor_for_M2_gg_QCD * additional_factor_for_nointerf * (event.ss_factor_1 + event.ss_factor_beta);
     //}   
    if (event.version == calc_weight_version::juan_code){
       return M2_QCD_juan;
    }
    if (event.version ==  calc_weight_version::juan_paper){
     return M2_QCD;
    } 
     return M2_QCD;
}


template <typename Number2>
Number2 calc_resonance_pseudo_new(event_t<Number2> event){

      //Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
      std::complex<Number2> PB = calc_P(event);
      std::complex<Number2> denomA;

      // std::complex<Number2> common_factor_for_M_A = 0.;
      Number2 ma_ga = event.higgs_mass * event.higgs_width;

     // std::complex<Number2> PB_real, PB_imag;
     // Number2 auxa = -1.5 * (1 - event.beta_ref_sq)/ 4. * pow(event.y_top, 2);
     // if (event.beta_ref_sq>0) {
     //       Number2 beta_ref = sqrt(event.beta_ref_sq);
     //       PB_real = auxa*(pow(log((1 + beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2));
     //       PB_imag = auxa*(-2 * TMath::Pi() * (log((1 + beta_ref) / (1 - beta_ref))));
     // } else {
     //       PB_real = -4 * auxa * pow(asin(event.sqrt_s/ 2 / constants<Number2>::m_t_ref_sq), 2);
     //       PB_imag = 0.;
     // }
      //PB = PB_real + PB_imag * constants<Number2>::i_cmplx;
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;

      std::cout << "denom A = " << denomA << std::endl;

      //common_factor_for_M_A = PB / denomA;
      Number2 common_factor_for_M_A_part = pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 resonance_pseudo = 0;
 
      switch (event.version){
         case calc_weight_version::juan_code:{
             //resonance_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
             Number2 mtrefOverE = 2 * constants<Number2>::m_t_ref / event.sqrt_s;
             // (pow(common_factor_for_M_gg_ss_interf * mtrefOverE, 2) + s
             resonance_pseudo = std::norm(PB / denomA) * pow(common_factor_for_M_A_part, 2) * event.ss_factor_1;
             std::cout << "mtrefOverE = " << mtrefOverE << std::endl;
             std::cout << "PB / denomA = " << (PB / denomA) << std::endl; 
             break;}
         case calc_weight_version::juan_paper:{
             Number2 factor_resonance = pow(constants<Number2>::G_F, 2) * constants<Number2>::m_t_ref_sq * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3));
             resonance_pseudo = factor_resonance * event.beta * std::norm(PB / denomA) * event.ss_factor_1;
             //resonance_pseudo = factor_resonance * event.beta * pow(PB.real(), 2) / std::norm(denomA) * event.ss_factor_1;
             break;}
         default:
             break;
      }
      Number2 factor_resonance = pow(constants<Number2>::G_F, 2) * constants<Number2>::m_t_ref_sq * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3));
      Number2 resonance_pseudo_paper = factor_resonance * event.beta * std::norm(PB / denomA) * event.ss_factor_1;
      std::cout << "resonance pseudo paper: " << resonance_pseudo_paper << std::endl;
      std::cout << "resonance pseudo paper * factor: " << resonance_pseudo_paper * 8. * TMath::Pi() * event.s / event.beta << std::endl;
      return resonance_pseudo;
}

template <typename Number2>
Number2 calc_interference_pseudo_new(event_t<Number2> event){

      Number2 common_factor_for_M_gg_ss_interf = (Number2) TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
      std::complex<Number2> PB;
      std::complex<Number2> denomA;

      //std::complex<Number2> common_factor_for_M_A = 0.;
      Number2 ma_ga = event.higgs_mass * event.higgs_width;

      std::complex<Number2> PB_real, PB_imag;
      Number2 auxa = -1.5 * (1 - event.beta_ref_sq)/ 4. * pow(event.y_top, 2);
      if (event.beta_ref_sq>0) {
            Number2 beta_ref = sqrt(event.beta_ref_sq);
            PB_real = auxa*(pow(log((1 + beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2));
            PB_imag = auxa*(-2 * TMath::Pi() * (log((1 + beta_ref) / (1 - beta_ref))));
      } else {
            PB_real = -4 * auxa * pow(asin(event.sqrt_s/ 2 / constants<Number2>::m_t_ref_sq), 2);
            PB_imag = 0.;
      }
      PB = PB_real + PB_imag * constants<Number2>::i_cmplx;
      std::cout << "P=" << PB << std::endl;
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;
      std::cout << "denomA=" << denomA << std::endl;
      std::cout << "denomA=" << denomA << std::endl;
      std::cout << "PB/denomA=" << PB / denomA << std::endl;
      

      //common_factor_for_M_A = PB / denomA;
      Number2 common_factor_for_M_A_part = pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 interference_pseudo = 0;
 
      switch (event.version){
         case calc_weight_version::juan_code:{
             // not fixed yet, still res and int or so?
             //interference_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
             interference_pseudo = (PB / denomA).real() * 2 * common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s * common_factor_for_M_A_part * event.ss_factor_1;
             break;}
         case calc_weight_version::juan_paper:{
             Number2 factor_interference = constants<Number2>::G_F * constants<Number2>::m_t_ref_sq / ( 48 * sqrt(2) * TMath::Pi());
             Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
             interference_pseudo = factor_interference * event.beta / beta_z_factor * (PB / denomA).real() * event.ss_factor_1;
             break;}
         default:
             break;
      }
      Number2 factor_interference = constants<Number2>::G_F * constants<Number2>::m_t_ref_sq / ( 48 * sqrt(2) * TMath::Pi());
      Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
      Number2 interference_pseudo_paper = factor_interference * event.beta / beta_z_factor * (PB / denomA).real() * event.ss_factor_1;
      std::cout << "interference pseudo paper: " << interference_pseudo_paper << std::endl;
      std::cout << "interference pseudo paper * factor: " << interference_pseudo_paper * 8. * TMath::Pi() * event.s / event.beta << std::endl;
      return interference_pseudo;
}

template<typename Number2>
Number2 calc_resonance_scalar_new(event_t<Number2> event){

     //Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
     std::complex<Number2> NB = calc_N(event);
     std::complex<Number2> denomH;

     //std::complex<Number2> common_factor_for_M_H = 0.;
     Number2 mh_gh = event.higgs_width * event.higgs_mass;

    // std::complex<Number2> NB_real, NB_imag;
    // Number2 auxh = 1.5 * (1 - event.beta_ref_sq) * pow(event.y_top, 2);
    // if (event.beta_ref_sq>0) {
    //       Number2 beta_ref = sqrt(event.beta_ref_sq);
    //       NB_real = auxh * (1 - event.beta_ref_sq / 4 * (pow( log((1+beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2) ));
    //       NB_imag = auxh * (event.beta_ref_sq / 2 * TMath::Pi() * (log((1+beta_ref)/(1-beta_ref))));
    // } else {
    //       NB_real = auxh * (1 + event.beta_ref_sq * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2));
    //       NB_imag = 0.;
    // }
    // NB = NB_real + NB_imag * constants<Number2>::i_cmplx;
     denomH = event.s - (Number2) pow(event.higgs_mass, 2) + mh_gh * constants<Number2>::i_cmplx;

     std::cout << "N = " << NB << std::endl;
     std::cout << "denomH = " << denomH << std::endl;

     //common_factor_for_M_H = NB / denomH;
     Number2 common_factor_for_M_H_part = pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6.* sqrt(3.) / 4. / TMath::Pi();
     
     Number2 resonance_scalar = 0;
     switch (event.version){
          // not fixed yet
          case calc_weight_version::juan_code:{
               //resonance_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H)
               //                          * event.ss_factor_beta;
               resonance_scalar = std::norm(NB / denomH) * pow(common_factor_for_M_H_part, 2) * event.ss_factor_beta;
               break;}
          case calc_weight_version::juan_paper:{
              Number2 factor_resonance = pow(constants<Number2>::G_F, 2) * constants<Number2>::m_t_ref_sq * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3)); 
              resonance_scalar = factor_resonance * pow(event.beta, 3) * std::norm(NB / denomH) * event.ss_factor_beta;
              //resonance_scalar = factor_resonance * event.beta * pow(NB.real(), 2) / std::norm(denomH) * event.ss_factor_beta;
              break;}
          default:
              break;
     }
      Number2 factor_resonance = pow(constants<Number2>::G_F, 2) * constants<Number2>::m_t_ref_sq * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3));
      Number2 resonance_scalar_paper = factor_resonance * pow(event.beta, 3) * std::norm(NB / denomH) * event.ss_factor_beta;
      std::cout << "resonance scalar paper: " << resonance_scalar_paper << std::endl;
      std::cout << "resonance scalar paper * factor: " << resonance_scalar_paper * 8. * TMath::Pi() * event.s / pow(event.beta, 3) << std::endl;
     return resonance_scalar;
}

template<typename Number2>
Number2 calc_interference_scalar_new(event_t<Number2> event){

     Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
     std::complex<Number2> NB;
     std::complex<Number2> denomH;

     std::complex<Number2> common_factor_for_M_H = 0.;
     Number2 mh_gh = event.higgs_width * event.higgs_mass;

     std::complex<Number2> NB_real, NB_imag;
     Number2 auxh = 1.5 * (1 - event.beta_ref_sq) * pow(event.y_top, 2);
     if (event.beta_ref_sq>0) {
           Number2 beta_ref = sqrt(event.beta_ref_sq);
           NB_real = auxh * (1 - event.beta_ref_sq / 4 * (pow( log((1+beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2) ));
           NB_imag = auxh * (event.beta_ref_sq / 2 * TMath::Pi() * (log((1+beta_ref)/(1-beta_ref))));
     } else {
           NB_real = auxh * (1 + event.beta_ref_sq * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2));
           NB_imag = 0.;
     }
     NB = NB_real + NB_imag * constants<Number2>::i_cmplx;
     denomH = event.s - (Number2) pow(event.higgs_mass, 2) + mh_gh * constants<Number2>::i_cmplx;

     common_factor_for_M_H = NB / denomH;
     Number2 common_factor_for_M_H_part = pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6.* sqrt(3.) / 4. / TMath::Pi();
     
     Number2 interference_scalar = 0;
     switch (event.version){
          // not fixed yet
          case calc_weight_version::juan_code:{
               //interference_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H) 
               //                          * event.ss_factor_beta;
               interference_scalar = (NB / denomH).real() * 2 * common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s * common_factor_for_M_H_part * event.ss_factor_beta;
               break;}
          case calc_weight_version::juan_paper:{
              Number2 factor_interference = constants<Number2>::G_F * constants<Number2>::m_t_ref_sq / (48 * sqrt(2.) * TMath::Pi()); 
              Number2 beta_z_factor = 1. - event.beta_sq * event.z_sq;
              interference_scalar = factor_interference * pow(event.beta, 3) / beta_z_factor * (NB / denomH).real() * event.ss_factor_beta;
              break;}
          default:
              break;
     }
      Number2 factor_interference = constants<Number2>::G_F * constants<Number2>::m_t_ref_sq / ( 48. * sqrt(2.) * TMath::Pi());
      Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
      Number2 interference_scalar_paper = factor_interference * pow(event.beta, 3) / beta_z_factor * (NB / denomH).real() * event.ss_factor_beta;
      std::cout << "interference scalar paper: " << interference_scalar_paper << std::endl;
      std::cout << "interference scalar paper * factor: " << interference_scalar_paper * 8. * TMath::Pi() * event.s / pow(event.beta, 3) << std::endl;
     return interference_scalar;
}

template <typename Number2>
Number2 M_2_bsm_ggHtt(event_t<Number2> event){


      Number2 resonance_scalar = 0;
      Number2 interference_scalar = 0;
      Number2 resonance_pseudo = 0;
      Number2 interference_pseudo = 0;
      if (event.higgs_type == higgs_type_t::scalar){
          if (event.res_int == res_int_t::resonance or event.res_int == res_int_t::both){
              resonance_scalar = calc_resonance_scalar_new(event);
              //resonance_scalar = calc_resonance_scalar_without_decay(event);
          }
          if (event.res_int == res_int_t::interference or event.res_int == res_int_t::both){
              interference_scalar = - calc_interference_scalar_new(event);
          }
      }
      else if (event.higgs_type == higgs_type_t::pseudo_scalar) {
          if (event.res_int == res_int_t::resonance or event.res_int == res_int_t::both){
              resonance_pseudo = calc_resonance_pseudo_new(event);
              //resonance_pseudo = calc_resonance_pseudo_without_decay(event);
          }
          if (event.res_int == res_int_t::interference or event.res_int == res_int_t::both){
              interference_pseudo = - calc_interference_pseudo_new(event);
          }
      }
             
//      Number2 qcd_opp_gluon = 0;
//      Number2 no_interference_term = 0;
//      Number2 interference_term_scalar = 0;
//      Number2 interference_term_pseudo = 0;
//      if (event.res_int == res_int_t::resonance or event.res_int == res_int_t::both){
//          qcd_opp_gluon = calc_qcd_opp_gluon(event);
//          no_interference_term = calc_no_interference_term(event);
//      }
//      if (event.res_int == res_int_t::interference or event.res_int == res_int_t::both){
//          if (event.higgs_type == higgs_type_t::scalar){
//              interference_term_scalar = calc_interference_scalar(event);
//          }
//          else if (event.higgs_type == higgs_type_t::pseudo_scalar) {
//              interference_term_pseudo = calc_interference_pseudo(event);
//          }
//      }
     
 
      std::cout << "resonance scalar: " << resonance_scalar << std::endl;
      std::cout << "interference scalar: " << interference_scalar << std::endl;
      std::cout << "resonance pseudo scalar: " << resonance_pseudo << std::endl;
      std::cout << "interference pseudo scalar: " << interference_pseudo << std::endl;  
//      std::cout << "qcd term: " << qcd_opp_gluon << std::endl;
//      std::cout << "no interference term: " << no_interference_term << std::endl;
//      std::cout << "interference term pseudo scalar: " << interference_term_pseudo << std::endl;
//      std::cout << "interference term scalar: " << interference_term_scalar << std::endl;  
      
      //std::cout << "M2_gg_QCD factor: " << common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor << std::endl;
      //std::cout << "M_gg_ss_nointerf factor: " << common_factor_for_M2_gg_ss_nointerf*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta) << std::endl;
      //std::cout << "M_gg_ss_interf M_A factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1 << std::endl;
      //std::cout << "M_gg_ss_interf M_H factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta << std::endl;

  
      // Number2 M_bsm = qcd_opp_gluon + no_interference_term + interference_term_scalar + interference_term_pseudo;
      // Number2 M_qcd = calc_M_qcd(event);

      //Number2 M_2_QCD = calc_M_2_qcd(event);
      //Number2 M_2_QCD = calc_QCD_without_decay(event);
      //Number2 M_2_QCD = calc_QCD_without_decay_dicus(event);
      Number2 M_2_bsm = resonance_scalar + interference_scalar + resonance_pseudo + interference_pseudo;
 
      //std::cout << "M_2_qcd: " << M_2_QCD << std::endl;  
      std::cout << "M_2_bsm: " << M_2_bsm << std::endl;  

      //Number2 weight = M_2_bsm / M_2_QCD;
      //Number2 weight = 1. / M_2_QCD;

      //std::cout << "weight: " << weight << std::endl;  
     
      if (isnan(M_2_bsm) or isinf(M_2_bsm)) {
            printf("What?? M2: %.3e; returning weight=0 !!\n", M_2_bsm);
            return 0.;
      }

 
      //if (isnan(weight) or isinf(weight)) {
      //      printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M_2_QCD, M_2_bsm);
      //      return 0.;
      //}
      
      return M_2_bsm;
}


template <typename Number2>
Number2 M_2_QCD_ggHtt(event_t<Number2> event){

      Number2 M_2_QCD = calc_M_2_qcd(event);
      
      if (isnan(M_2_QCD) or isinf(M_2_QCD)) {
            printf("What?? M2QCD: %.3e, returning weight=0 !!\n", M_2_QCD);
            return 0.;
      }
      
      return M_2_QCD;
}

template <typename Number1, typename Number2>
Number2 compute_weight(Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
                       Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
                       Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
                       Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m,
                       calc_weight_version version, higgs_type_t higgs_type, Number2 mass, Number2 width, res_int_t res_int)
{

  // Kinematics in this event
 TLorentzVector vec_4_t;
 vec_4_t.SetPtEtaPhiM(pTop_pt, pTop_eta, pTop_phi, pTop_m);
 TLorentzVector vec_4_tbar;
 vec_4_tbar.SetPtEtaPhiM(aTop_pt, aTop_eta, aTop_phi, aTop_m);
 TLorentzVector vec_4_l;
 vec_4_l.SetPtEtaPhiM(pLep_pt, pLep_eta, pLep_phi, pLep_m);
 TLorentzVector vec_4_lbar;
 vec_4_lbar.SetPtEtaPhiM(aLep_pt, aLep_eta, aLep_phi, aLep_m);

 //std::cout << "top x: " << vec_4_t.X() << ", y: " << vec_4_t.Y() << ", z: " << vec_4_t.Z() << ", E " << vec_4_t.E() << std::endl;
 //std::cout << "antitop x: " << vec_4_tbar.X() << ", y: " << vec_4_tbar.Y() << ", z: " << vec_4_tbar.Z() << ", E " << vec_4_tbar.E() << std::endl;
 //std::cout << "lepton x: " << vec_4_l.X() << ", y: " << vec_4_l.Y() << ", z: " << vec_4_l.Z() << ", E " << vec_4_l.E() << std::endl;
 //std::cout << "antilepton x: " << vec_4_lbar.X() << ", y: " << vec_4_lbar.Y() << ", z: " << vec_4_lbar.Z() << ", E " << vec_4_lbar.E() << std::endl;
  

  std::cout << std::endl << " === New Event === " << std::endl;

  event_t<Number2> event =  event_t(vec_4_t, vec_4_tbar, vec_4_l, vec_4_lbar, higgs_type, mass, width, res_int, version);

  if (res_int == res_int_t::resonance)
    std::cout << "res" << std::endl;
  if (version == calc_weight_version::juan_code)
    std::cout << "juan code" << std::endl;

  HTT_Input<Number2> httInput;
  if (higgs_type == higgs_type_t::pseudo_scalar)
    httInput.HIGGS_OPTION = 1;
  if (higgs_type == higgs_type_t::scalar)
    httInput.HIGGS_OPTION = 0;
  httInput.WIDTH_OPTION = 0;
  httInput.YTOP = 1.0;
  httInput.MH = mass;
  httInput.GH = width;
  httInput.MA = mass;
  httInput.GA = width;

  httInput.P4gen_t[0].SetPtEtaPhiM(pTop_pt, pTop_eta, pTop_phi, pTop_m);
  httInput.P4gen_t[1].SetPtEtaPhiM(aTop_pt, aTop_eta, aTop_phi, aTop_m);
  httInput.P4gen_d[1].SetPtEtaPhiM(pLep_pt, pLep_eta, pLep_phi, pLep_m);
  httInput.P4gen_d[0].SetPtEtaPhiM(aTop_pt, aTop_eta, aTop_phi, aTop_m);
 
  //Number2 weight_juan = weight_ggHtt(httInput);

  Number2 M_2_bsm = M_2_bsm_ggHtt(event);

  Number2 M_2_QCD = M_2_QCD_ggHtt(event);

  Number2 weight = M_2_bsm / M_2_QCD;

 std::cout << "mass = " << mass << std::endl;
 std::cout << "width = " << width << std::endl;
 std::cout << "higgs option = " << httInput.HIGGS_OPTION << std::endl;
  
 std::cout << "nachher: " << std::endl; 
 std::cout << "httInput.P4gen_t[0].SetPtEtaPhiM(" << pTop_pt << ", " << pTop_eta << ", " << pTop_phi << ", " << pTop_m << ");" << std::endl;  
 std::cout << "httInput.P4gen_t[1].SetPtEtaPhiM(" << aTop_pt << ", " << aTop_eta << ", " << aTop_phi << ", " << aTop_m << ");" << std::endl;  
 std::cout << "httInput.P4gen_d[1].SetPtEtaPhiM(" << pLep_pt << ", " << pLep_eta << ", " << pLep_phi << ", " << pLep_m << ");" << std::endl;  
 std::cout << "httInput.P4gen_d[0].SetPtEtaPhiM(" << aLep_pt << ", " << aLep_eta << ", " << aLep_phi << ", " << aLep_m << ");" << std::endl;  
 std::cout << "compare weight: " << weight << std::endl; 
 std::cout << "compare M2 QCD: " << M_2_QCD << std::endl; 
 std::cout << "compare M2 BSM: " << M_2_bsm << std::endl; 

 std::cout << "compare M2 QCD * factor: " << M_2_QCD  * 8. * event.s * TMath::Pi() / pow(event.beta, 3) << std::endl; 
 std::cout << "compare M2 BSM * factor: " << M_2_bsm  * 8. * event.s * TMath::Pi() / pow(event.beta, 3) << std::endl; 


 // 
 // if (isnan(weight) or isinf(weight)) {
 //      printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M_2_QCD, M_2_bsm);
 //      return 0.;
 // }
  
  return weight;

}




template <typename Number1, typename Number2>
auto calculate_weight(calc_weight_version version, higgs_type_t higgs_type, float mass, float width, res_int_t res_int)
{

  return [version, higgs_type, mass, width, res_int] (Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
             Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
             Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
             Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m)
  {

       return compute_weight<Number1, Number2>(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                          aTop_pt, aTop_eta, aTop_phi, aTop_m,
                          pLep_pt, pLep_eta, pLep_phi, pLep_m,
                          aLep_pt, aLep_eta, aLep_phi, aLep_m, 
                          version, higgs_type, mass, width, res_int);
  };
}



// ========= methods to return qcd and bsm separately =========
template <typename Number1, typename Number2>
Number2 compute_qcd(Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
                       Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
                       Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
                       Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m,
                       calc_weight_version version, higgs_type_t higgs_type, Number2 mass, Number2 width, res_int_t res_int)
{

  // Kinematics in this event
  TLorentzVector vec_4_t;
  vec_4_t.SetPtEtaPhiM(pTop_pt, pTop_eta, pTop_phi, pTop_m);
  TLorentzVector vec_4_tbar;
  vec_4_tbar.SetPtEtaPhiM(aTop_pt, aTop_eta, aTop_phi, aTop_m);
  TLorentzVector vec_4_l;
  vec_4_l.SetPtEtaPhiM(pLep_pt, pLep_eta, pLep_phi, pLep_m);
  TLorentzVector vec_4_lbar;
  vec_4_lbar.SetPtEtaPhiM(aLep_pt, aLep_eta, aLep_phi, aLep_m);

  event_t<Number2> event =  event_t(vec_4_t, vec_4_tbar, vec_4_l, vec_4_lbar, higgs_type, mass, width, res_int, version);
 
  Number2 M_2_QCD = M_2_QCD_ggHtt(event);
 
  return M_2_QCD;

}

template <typename Number1, typename Number2>
auto calculate_qcd(calc_weight_version version, higgs_type_t higgs_type, float mass, float width, res_int_t res_int)
{

  return [version, higgs_type, mass, width, res_int] (Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
             Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
             Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
             Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m)
  {

       return compute_qcd<Number1, Number2>(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                          aTop_pt, aTop_eta, aTop_phi, aTop_m,
                          pLep_pt, pLep_eta, pLep_phi, pLep_m,
                          aLep_pt, aLep_eta, aLep_phi, aLep_m, 
                          version, higgs_type, mass, width, res_int);
  };
}


template <typename Number1, typename Number2>
Number2 compute_bsm(Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
                       Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
                       Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
                       Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m,
                       calc_weight_version version, higgs_type_t higgs_type, Number2 mass, Number2 width, res_int_t res_int)
{

  // Kinematics in this event
  TLorentzVector vec_4_t;
  vec_4_t.SetPtEtaPhiM(pTop_pt, pTop_eta, pTop_phi, pTop_m);
  TLorentzVector vec_4_tbar;
  vec_4_tbar.SetPtEtaPhiM(aTop_pt, aTop_eta, aTop_phi, aTop_m);
  TLorentzVector vec_4_l;
  vec_4_l.SetPtEtaPhiM(pLep_pt, pLep_eta, pLep_phi, pLep_m);
  TLorentzVector vec_4_lbar;
  vec_4_lbar.SetPtEtaPhiM(aLep_pt, aLep_eta, aLep_phi, aLep_m);

  event_t<Number2> event =  event_t(vec_4_t, vec_4_tbar, vec_4_l, vec_4_lbar, higgs_type, mass, width, res_int, version);
 
  Number2 M_2_bsm = M_2_bsm_ggHtt(event);

  return M_2_bsm;

}
template <typename Number1, typename Number2>
auto calculate_bsm(calc_weight_version version, higgs_type_t higgs_type, float mass, float width, res_int_t res_int)
{

  return [version, higgs_type, mass, width, res_int] (Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
             Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
             Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
             Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m)
  {

       return compute_bsm<Number1, Number2>(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                          aTop_pt, aTop_eta, aTop_phi, aTop_m,
                          pLep_pt, pLep_eta, pLep_phi, pLep_m,
                          aLep_pt, aLep_eta, aLep_phi, aLep_m, 
                          version, higgs_type, mass, width, res_int);
  };
}
