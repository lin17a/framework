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



enum class calc_weight_version{
  juan_paper, juan_code, juan_paper_different_M2
};

enum class higgs_type_t{
  scalar, pseudo_scalar
};

enum class res_int_t{
  resonance, interference, both
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
    Number2 beta_sq_ref;
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
        this->s = s;
        Number2 m_t = vec_4_t.M();
        this->m_t = m_t;
        Number2 m_tbar = vec_4_tbar.M();
        this->m_tbar = m_tbar;
        Number2 beta_sq = (1 - (m_t - m_tbar) * (m_t - m_tbar) / s) * (1 - (m_t + m_tbar) * (m_t + m_tbar) / s);
        this->beta_sq = beta_sq;
        this->beta = sqrt(beta);
        //Number2 beta = sqrt(beta_sq);
        //Number2 beta3 = beta*beta_sq;
        //std::cout << "top mass: " << m_t << std::endl;
        Number2 mtrefOverE = 2 * constants<Number2>::m_t_ref / sqrt_s;
        //std::cout << "square root of s: " << sqrt_s << std::endl;
        //std::cout << "mt over E: " << mtrefOverE << std::endl; 
        Number2 beta_sq_ref = 1 - mtrefOverE * mtrefOverE;
        this->beta_sq_ref = beta_sq_ref;
        //std::cout << "beta squared ref: " << beta_sq_ref << std::endl;
        Number2 beta_ref = 0;
        if (beta_sq_ref > 0) {
            beta_ref = sqrt(beta_sq_ref);
        }
        this->beta_ref = beta_ref;

        TVector3 zBoost(0.,0.,-vec_4_ttbar.Pz()/vec_4_ttbar.E()); vec_4_ttbar.Boost(zBoost);
        TVector3 transverseBoost(-vec_4_ttbar.Px()/vec_4_ttbar.E(),-vec_4_ttbar.Py()/vec_4_ttbar.E(),0.); vec_4_ttbar.Boost(transverseBoost);
        vec_4_t.Boost(zBoost); vec_4_t.Boost(transverseBoost);
        vec_4_lbar.Boost(zBoost); vec_4_lbar.Boost(transverseBoost);
        vec_4_tbar.Boost(zBoost); vec_4_tbar.Boost(transverseBoost);
        vec_4_l.Boost(zBoost); vec_4_l.Boost(transverseBoost);
        
        Number2 z = vec_4_t.CosTheta();
        this->z = z;
        Number2 z_sq = z * z;
        this->z_sq = z_sq;  
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
        TLorentzVector vec_4_lbar_trest = vec_4_lbar;
        TVector3 t2Boost(-vec_4_tbar.Px()/vec_4_tbar.E(),-vec_4_tbar.Py()/vec_4_tbar.E(),-vec_4_tbar.Pz()/vec_4_tbar.E()); 
        vec_4_l.Boost(t2Boost);
        TLorentzVector vec_4_l_trest = vec_4_l;
        this->vec_4_l_trest = vec_4_l_trest;
        this->vec_4_lbar_trest = vec_4_lbar_trest;    
        

        Number2 c1 = vec_4_lbar_trest.CosTheta();
        //std::cout << "c1: " << c1 << std::endl;
        Number2 c2 = vec_4_l_trest.CosTheta();
        Number2 s1 = sqrt(1-c1*c1);
        Number2 s2 = sqrt(1-c2*c2);
        Number2 phi1 = vec_4_lbar_trest.Phi();
        Number2 phi2 = vec_4_l_trest.Phi();
     
        // os factor deleted one sqrt in last line 
        Number2 os_factor_juan_paper = (1 + z_sq + (1 - z_sq) * (1 - beta_sq_ref)) * (1. - c1*c2)
              + 2 * (1 - beta_sq_ref) * (1 - z_sq) * c1 * c2
              - beta_sq_ref * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
              - 2 * (1 - beta_sq_ref) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
              + 2 * sqrt(1 - beta_sq_ref) * z * (1 - z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
        this->os_factor_juan_paper = os_factor_juan_paper;

        // old os factor
        Number2 os_factor_juan_code = (1 + z_sq + (1 - z_sq) * (1 - beta_sq_ref)) * (1. - c1 * c2) 
              + 2 * (1 - beta_sq_ref) * (1 - z_sq) * c1 * c2
              - beta_sq_ref * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
              - 2 * (1 - beta_sq_ref) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
              + 2 * sqrt(1 - beta_sq_ref) * z * sqrt(1 - z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
        this->os_factor_juan_code = os_factor_juan_code;

        //std::cout << "os_factor: " << os_factor << std::endl;
        Number2 ss_factor_1 = 1 + c1 * c2 + s1 * s2 * cos(phi1 - phi2);
        this->ss_factor_1 = ss_factor_1;
        //std::cout << "ss_factor_1: " << ss_factor_1 << std::endl;
        Number2 ss_factor_beta = beta_sq * (1 + c1 * c2 - s1 * s2 * cos(phi1 - phi2));
        this->ss_factor_beta = ss_factor_beta;
    }

};


// Reweighting gg-> H/A -> ttbar events
// Input is an HTT_Input structure; output is the event weight
template <typename Number2>
Number2 weight_ggHtt(event_t<Number2> event);

// Simple code to decide whether the initial state should be gg or not
// pdgId1 is the the PDG id of first parton, along the +Z direction
// pdgId2 is the the PDG id of other parton, along the -Z direction
// pz_hard_radiation is the pz of the system made of the hard radiated partons recoiling 
template <typename Number2>
bool is_gg_initial_state(const int& pdgId1, const int& pdgId2, const Number2& pz_hard_radiation);


template<typename Number2>
Number2 calc_no_interference_scalar(event_t<Number2> event){

     Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
     std::complex<Number2> NB;
     std::complex<Number2> denomH;

     std::complex<Number2> common_factor_for_M_H = 0.;
     Number2 mh_gh_partial = 0.;
     Number2 mh_gh = 0.;
     if (event.higgs_mass > 2* constants<Number2>::m_t_ref_sq) {
           Number2 beta3_tdecay = pow(1-4* pow(constants<Number2>::m_t_ref, 2) / pow(event.higgs_mass, 2), 1.5);
           mh_gh_partial += 3 * constants<Number2>::G_F * constants<Number2>::m_t_ref_sq * pow(event.higgs_mass, 2) / (4 * TMath::Pi() * sqrt(2.)) * beta3_tdecay * pow(event.y_top, 2);
     }

     mh_gh = event.higgs_width * event.higgs_mass;

     std::complex<Number2> NB_real, NB_imag;
     Number2 auxh = 1.5 * (1 - event.beta_sq_ref) * pow(event.y_top, 2);
     if (event.beta_sq_ref>0) {
           Number2 beta_ref = sqrt(event.beta_sq_ref);
           NB_real = auxh * (1 - event.beta_sq_ref / 4 * (pow( log((1+beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2) ));
           NB_imag = auxh * (event.beta_sq_ref / 2 * TMath::Pi() * (log((1+beta_ref)/(1-beta_ref))));
     } else {
           NB_real = auxh * (1 + event.beta_sq_ref * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2));
           NB_imag = 0.;
     }
     NB = NB_real + NB_imag * constants<Number2>::i_cmplx;
     denomH = event.s - (Number2) pow(event.higgs_mass, 2) + mh_gh * constants<Number2>::i_cmplx;

     common_factor_for_M_H = NB / denomH;
     common_factor_for_M_H *= pow(event.s,1.5) * constants<Number2>::m_t_ref_sq * constants<Number2>::G_F / 6.* sqrt(3.) / 4. / TMath::Pi();
     
     Number2 no_interference_scalar = 0;
     if (event.version == calc_weight_version::juan_code){
        no_interference_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H) \
                                         * event.ss_factor_beta;
     }
     else if (event.version == calc_weight_version::juan_paper){
         Number2 beta_z_factor = 1 - event.beta_sq_ref * event.z_sq;
         Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
         Number2 factor_interference = 4 / 3 * pow(constants<Number2>::m_t_ref, 2) * pow(TMath::Pi(), 3) * pow(constants<Number2>::alpha_s, 2) / (pow(event.s, 2) * beta_z_factor_sq);
         Number2 factor_in_norm_interference = sqrt(2.) * pow(event.s, 2) * constants<Number2>::G_F * beta_z_factor / (16 * pow(TMath::Pi(), 2));
         
         no_interference_scalar = factor_interference * pow(event.beta, 3) \
                                  * pow(std::norm((Number2) 1. - factor_in_norm_interference * (NB / denomH)), 2) \
                                  * event.ss_factor_beta;
     }
     return no_interference_scalar;
}

template <typename Number2>
Number2 calc_no_interference_pseudo(event_t<Number2> event){

      Number2 common_factor_for_M_gg_ss_interf = TMath::Pi() / sqrt(6.) / (1 - event.beta_sq * event.z_sq);
      std::complex<Number2> PB;
      std::complex<Number2> denomA;

      std::complex<Number2> common_factor_for_M_A = 0.;
      Number2 ma_ga_partial = 0.;
      Number2 ma_ga = 0.;
      if (event.higgs_mass > 2 * constants<Number2>::m_t_ref_sq) {
            Number2 beta_tdecay = pow(1 - 4 * constants<Number2>::m_t_ref_sq / pow(event.higgs_mass, 2) ,0.5);
            ma_ga_partial += 3 * constants<Number2>::G_F * constants<Number2>::m_t_ref_sq * pow(event.higgs_mass, 2) / ( 4 * TMath::Pi() * sqrt(2.)) * beta_tdecay * pow(event.y_top, 2);
      }
      ma_ga = event.higgs_mass * event.higgs_width;

      std::complex<Number2> PB_real, PB_imag;
      Number2 auxa = -1.5 * (1 - event.beta_sq_ref)/ 4. * pow(event.y_top, 2);
      if (event.beta_sq_ref>0) {
            Number2 beta_ref = sqrt(event.beta_sq_ref);
            PB_real = auxa*(pow(log((1 + beta_ref) / (1-beta_ref)), 2) - pow(TMath::Pi(), 2));
            PB_imag = auxa*(-2 * TMath::Pi() * (log((1 + beta_ref) / (1 - beta_ref))));
      } else {
            PB_real = -4 * auxa * pow(asin(event.sqrt_s/ 2 / constants<Number2>::m_t_ref_sq), 2);
            PB_imag = 0.;
      }
      PB = PB_real + PB_imag * constants<Number2>::i_cmplx;
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;

      common_factor_for_M_A = PB / denomA;
      common_factor_for_M_A *= pow(event.s,1.5) * constants<Number2>::m_t_ref_sq * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 no_interference_pseudo = 0;
      if (event.version == calc_weight_version::juan_code){
         no_interference_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
                                          * event.ss_factor_1;
      }
      else if (event.version == calc_weight_version::juan_paper){
          Number2 beta_z_factor = 1 - event.beta_sq_ref * event.z_sq;
          Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
          Number2 factor_interference = 4 / 3 * pow(constants<Number2>::m_t_ref, 2) * pow(TMath::Pi(), 3) * pow(constants<Number2>::alpha_s, 2) / (pow(event.s, 2) * beta_z_factor_sq);
          Number2 factor_in_norm_interference = sqrt(2.) * pow(event.s, 2) * constants<Number2>::G_F * beta_z_factor / (16 * pow(TMath::Pi(), 2));
          
          no_interference_pseudo = factor_interference * event.beta 
                                   * pow(std::norm((Number2) 1. - factor_in_norm_interference * (PB / denomA)), 2) 
                                   * event.ss_factor_1;
      }
      return no_interference_pseudo;
}

template <typename Number2>
Number2 calc_no_interference_term(event_t<Number2> event){
    const Number2 alpha_s = 0.12;
    Number2 common_factor_for_M2_gg_ss_nointerf = pow(TMath::Pi(), 2) / 12 * (5 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    Number2 additional_factor_for_nointerf = 0;
    Number2 m_t_ref = 172.5;
    if (event.version == calc_weight_version::juan_paper){
        additional_factor_for_nointerf = 8 * pow(m_t_ref, 2) * TMath::Pi() * alpha_s * event.beta / pow(event.s, 2);
    }
    else{
        additional_factor_for_nointerf = 1 - pow(event.beta_ref, 2);
    }
    Number2 no_interference_term = common_factor_for_M2_gg_ss_nointerf * additional_factor_for_nointerf * (event.ss_factor_1 + event.ss_factor_beta);
    return no_interference_term;
}


template <typename Number2>
Number2 calc_qcd_opp_gluon(event_t<Number2> event){
    const Number2 alpha_s = 0.12; 
    Number2 common_factor_for_M2_gg_QCD = pow(TMath::Pi(), 2) / 12 * (7 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    Number2 qcd_term = 0;
    if (event.version == calc_weight_version::juan_paper) {
        qcd_term = pow(alpha_s, 2) * pow(event.beta, 3) / (8 * event.s * TMath::Pi()) * common_factor_for_M2_gg_QCD * (1 - event.z_sq) * event.os_factor_juan_paper;
    }
    else {
        qcd_term = common_factor_for_M2_gg_QCD * pow(event.beta, 2) * (1 - event.z_sq) * event.os_factor_juan_code; 
    }
    return qcd_term;
}

template <typename Number2>
Number2 calc_M_qcd(event_t<Number2> event){
    Number2 M2_QCD = 0;
    Number2 common_factor_for_M2_gg_QCD = pow(TMath::Pi(), 2) / 12 * (7 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    Number2 additional_factor_for_nointerf = 8 * pow(constants<Number2>::m_t_ref, 2) * TMath::Pi() * constants<Number2>::alpha_s * event.beta / pow(event.s, 2);
    if (event.version == calc_weight_version::juan_code){
        M2_QCD = common_factor_for_M2_gg_QCD * pow(event.beta, 2) * (1 - event.z_sq) * event.os_factor_juan_code \
                 + common_factor_for_M2_gg_QCD * (1 - event.beta_sq_ref) * (event.ss_factor_1 + event.ss_factor_beta);
    }
    else if (event.version == calc_weight_version::juan_paper){
        // just to try it out, new M for QCD, looking at the other changes I made 
        // first part like the qcd part in BSM term and second part as it is, but (1-beta^2) becomes the same factor as in the no interference part
        M2_QCD = common_factor_for_M2_gg_QCD * pow(constants<Number2>::alpha_s, 2) * pow(event.beta, 3) / (8 * TMath::Pi() * event.s) * (1 - event.z_sq) \
                 * event.os_factor_juan_paper \
                 + common_factor_for_M2_gg_QCD * additional_factor_for_nointerf * (event.ss_factor_1 + event.ss_factor_beta);
    }
    return M2_QCD;
}

template <typename Number2>
Number2 weight_ggHtt(event_t<Number2> event){


      Number2 qcd_opp_gluon = 0;
      Number2 interference_term = 0;
      Number2 no_interference_term_scalar = 0;
      Number2 no_interference_term_pseudo = 0;
      if (event.res_int == res_int_t::resonance or event.res_int == res_int_t::both){
          qcd_opp_gluon = calc_qcd_opp_gluon(event);
          if (event.higgs_type == higgs_type_t::scalar){
              no_interference_term_scalar = calc_no_interference_scalar(event);
          }
          else if (event.higgs_type == higgs_type_t::pseudo_scalar) {
              no_interference_term_pseudo = calc_no_interference_pseudo(event);
          }
      }
      if (event.res_int == res_int_t::interference or event.res_int == res_int_t::both){
          interference_term = calc_no_interference_term(event);
      }
             
      
      std::cout << "qcd term: " << qcd_opp_gluon << std::endl;
      std::cout << "no interference term: " << interference_term << std::endl;
      std::cout << "interference term pseudo scalar: " << no_interference_term_pseudo << std::endl;
      std::cout << "interference term scalar: " << no_interference_term_scalar << std::endl;  
      
      //std::cout << "M2_gg_QCD factor: " << common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor << std::endl;
      //std::cout << "M_gg_ss_nointerf factor: " << common_factor_for_M2_gg_ss_nointerf*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta) << std::endl;
      //std::cout << "M_gg_ss_interf M_A factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1 << std::endl;
      //std::cout << "M_gg_ss_interf M_H factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta << std::endl;

  
      Number2 M_bsm = qcd_opp_gluon + interference_term + no_interference_term_scalar + no_interference_term_pseudo;
      Number2 M_qcd = calc_M_qcd(event);
 
      std::cout << "M_qcd: " << M_qcd << std::endl;  

      Number2 weight = M_bsm / M_qcd;
      
      if (isnan(weight)) {
            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M_qcd, M_bsm);
            return 0.;
      }
      return weight;
}



//template <typename Number2>
//Number2 weight_ggHtt_juan_paper(const HTT_Input<Number2>& httInput, res_int_t res_int) {
//      // Check that the chosen options are valid
//      if (httInput.HIGGS_OPTION>2 || httInput.WIDTH_OPTION>3) {
//            printf("What? HIGGS_OPTION=%d, WIDTH_OPTION=%d; returning weight=1 !!\n", httInput.HIGGS_OPTION, httInput.WIDTH_OPTION);
//            return 1.;
//      }
//
//      // define all important variables
//      Number2 m_H_ref_sq = pow(httInput.MH,2);
//      Number2 m_A_ref_sq = pow(httInput.MA,2);
//      // coupling constant
//      Number2 ytop_sq = pow(httInput.YTOP,2);
//
//      TLorentzVector vec_4_t = httInput.P4gen_t[0]; 
//      TLorentzVector vec_4_tbar = httInput.P4gen_t[1]; 
//      TLorentzVector vec_4_lbar = httInput.P4gen_d[0]; 
//      TLorentzVector vec_4_l = httInput.P4gen_d[1]; 
//
//      const Number2 m_t_ref = 172.5;
//      const Number2 m_t_ref_sq = m_t_ref*m_t_ref;
//      //Fermi-Konstante
//      const Number2 G_F = 1.16637876e-5;
//      const Number2 alpha_s = 0.12;
//      const Number2 pi = TMath::Pi();
//      const Number2 pi_sq = pi*pi;
//      const Number2 sqrt_2 = sqrt(2.);
//
//      TLorentzVector vec_4_ttbar = vec_4_t+vec_4_tbar;
//      Number2 sqrt_s = vec_4_ttbar.M();
//      Number2 s = sqrt_s*sqrt_s;
//      Number2 m_t = httInput.P4gen_t[0].M();
//      Number2 m_tbar = httInput.P4gen_t[1].M();
//      Number2 beta_sq = (1-(m_t-m_tbar)*(m_t-m_tbar)/s)*(1-(m_t+m_tbar)*(m_t+m_tbar)/s);
//      Number2 beta = sqrt(beta_sq);
//      //Number2 beta3 = beta*beta_sq;
//      std::cout << "top mass: " << m_t << std::endl;
//      Number2 mtrefOverE = 2 * m_t_ref / sqrt_s;
//      std::cout << "square root of s: " << sqrt_s << std::endl;
//      std::cout << "mt over E: " << mtrefOverE << std::endl; 
//      Number2 beta_sq_ref = 1 - mtrefOverE * mtrefOverE;
//      std::cout << "beta squared ref: " << beta_sq_ref << std::endl;
//      Number2 beta_ref = 0;
//      if (beta_sq_ref > 0) {
//          beta_ref = sqrt(beta_sq_ref);
//      }
//
//      TVector3 zBoost(0.,0.,-vec_4_ttbar.Pz()/vec_4_ttbar.E()); vec_4_ttbar.Boost(zBoost);
//      TVector3 transverseBoost(-vec_4_ttbar.Px()/vec_4_ttbar.E(),-vec_4_ttbar.Py()/vec_4_ttbar.E(),0.); vec_4_ttbar.Boost(transverseBoost);
//      vec_4_t.Boost(zBoost); vec_4_t.Boost(transverseBoost);
//      vec_4_lbar.Boost(zBoost); vec_4_lbar.Boost(transverseBoost);
//      vec_4_tbar.Boost(zBoost); vec_4_tbar.Boost(transverseBoost);
//      vec_4_l.Boost(zBoost); vec_4_l.Boost(transverseBoost);
//      
//      Number2 z = vec_4_t.CosTheta();
//      Number2 z_sq = z*z;
//
//      Number2 minusPhi1Plane = -vec_4_t.Phi();
//      TVector3 zVect(0.,0.,1.);
//      vec_4_t.Rotate(minusPhi1Plane,zVect);
//      vec_4_lbar.Rotate(minusPhi1Plane,zVect);
//      vec_4_tbar.Rotate(minusPhi1Plane,zVect);
//      vec_4_l.Rotate(minusPhi1Plane,zVect);
//      Number2 t1Yangle = -atan2(vec_4_t.Px(),vec_4_t.Pz());
//      TVector3 yVect(0.,1.,0.);
//      vec_4_t.Rotate(t1Yangle,yVect);
//      vec_4_lbar.Rotate(t1Yangle,yVect);
//      vec_4_tbar.Rotate(t1Yangle,yVect);
//      vec_4_l.Rotate(t1Yangle,yVect);
//
//      TVector3 t1Boost(-vec_4_t.Px()/vec_4_t.E(),-vec_4_t.Py()/vec_4_t.E(),-vec_4_t.Pz()/vec_4_t.E()); 
//      vec_4_lbar.Boost(t1Boost);
//      TLorentzVector P4gen_d1_trest = vec_4_lbar;
//      TVector3 t2Boost(-vec_4_tbar.Px()/vec_4_tbar.E(),-vec_4_tbar.Py()/vec_4_tbar.E(),-vec_4_tbar.Pz()/vec_4_tbar.E()); 
//      vec_4_l.Boost(t2Boost);
//      TLorentzVector vec_4_l_trest = vec_4_l;
//
//
//      Number2 common_factor_for_M2_gg_QCD = pi_sq/12*(7+9*beta_sq*z_sq)/pow(1-beta_sq*z_sq,2);
//      Number2 common_factor_for_M2_gg_ss_nointerf = pi_sq/12*(5+9*beta_sq*z_sq)/pow(1-beta_sq*z_sq,2);
//      Number2 common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta_sq*z_sq);
//
//
//      using namespace std::complex_literals;
//      //const std::complex<Number2> i_cmplx = 1i;
//      //const std::complex<Number2>{0.0, static_cast<double>(d)};
//      const std::complex i_cmplx = std::complex<Number2>(0., 1.L);
//      //const std::complex<Number2> i_cmplx = static_cast<Number2>(1);
//
//
//      std::complex<Number2> NB;
//      std::complex<Number2> denomH;
//
//      std::complex<Number2> common_factor_for_M_H = 0.;
//      Number2 mh_gh_partial = 0.;
//      Number2 mh_gh = 0.;
//      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
//            if (httInput.MH>2*m_t_ref_sq) {
//                  Number2 beta3_tdecay = pow(1-4*m_t_ref_sq/m_H_ref_sq,1.5);
//                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
//                        mh_gh_partial += 3*G_F*m_t_ref_sq*s/(4*pi*sqrt_2)*beta3_tdecay*ytop_sq;
//                  } else {
//                        mh_gh_partial += 3*G_F*m_t_ref_sq*m_H_ref_sq/(4*pi*sqrt_2)*beta3_tdecay*ytop_sq;
//                  }
//            }
//
//            if (httInput.WIDTH_OPTION==1) {
//                  mh_gh = httInput.GH*s/httInput.MH;
//            } else if (httInput.WIDTH_OPTION==0) {
//                  mh_gh = httInput.GH*httInput.MH;
//            } else {
//                  mh_gh = mh_gh_partial;
//            }
//
//            Number2 NB_real, NB_imag;
//            Number2 auxh = 1.5 * (1-beta_sq_ref) * ytop_sq;
//            if (beta_sq_ref>0) {
//                  Number2 betaRef = sqrt(beta_sq_ref);
//                  NB_real = auxh * (1 - beta_sq_ref / 4 * (pow( log((1+betaRef) / (1-betaRef)), 2) - pi_sq ));
//                  NB_imag = auxh * (beta_sq_ref / 2 * pi * (log((1+betaRef)/(1-betaRef))));
//            } else {
//                  NB_real = auxh * (1 + beta_sq_ref * pow(asin(sqrt_s/2/m_t_ref_sq), 2));
//                  NB_imag = 0.;
//            }
//            NB = NB_real + NB_imag * i_cmplx;
//            denomH = s - m_H_ref_sq + mh_gh * i_cmplx;
//
//            common_factor_for_M_H = NB / denomH;
//            common_factor_for_M_H *= pow(s,1.5)*m_t_ref_sq*G_F/6.*sqrt(3.)/4./pi;
//      }
//
//      std::complex<Number2> PB;
//      std::complex<Number2> denomA;
//
//      std::complex<Number2> common_factor_for_M_A = 0.;
//      Number2 ma_ga_partial = 0.;
//      Number2 ma_ga = 0.;
//      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
//            if (httInput.MA>2*m_t_ref_sq) {
//                  Number2 beta_tdecay = pow(1-4*m_t_ref_sq/m_A_ref_sq,0.5);
//                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
//                        ma_ga_partial += 3*G_F*m_t_ref_sq*s/(4*pi*sqrt_2)*beta_tdecay*ytop_sq;
//                  } else {
//                        ma_ga_partial += 3*G_F*m_t_ref_sq*m_A_ref_sq/(4*pi*sqrt_2)*beta_tdecay*ytop_sq;
//                  }
//            }
//            if (httInput.WIDTH_OPTION==1) {
//                        ma_ga = httInput.GA*s/httInput.MA;
//            } else if (httInput.WIDTH_OPTION==0) {
//                        ma_ga = httInput.GA*httInput.MA;
//            } else {
//                        ma_ga = ma_ga_partial;
//            }
//
//            Number2 PB_real, PB_imag;
//            Number2 auxa = -1.5*(1-beta_sq_ref)/4.*ytop_sq;
//            if (beta_sq_ref>0) {
//                  Number2 betaRef = sqrt(beta_sq_ref);
//                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi_sq);
//                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
//            } else {
//                  PB_real = -4*auxa*pow(asin(sqrt_s/2/m_t_ref_sq),2);
//                  PB_imag = 0.;
//            }
//            PB = PB_real + PB_imag * i_cmplx;
//            denomA = s - m_A_ref_sq + ma_ga * i_cmplx;
//
//            common_factor_for_M_A = PB / denomA;
//            common_factor_for_M_A *= pow(s,1.5)*m_t_ref_sq*G_F/6.*sqrt(3.)/4./pi;
//      }
//
//      // Printing global input information (once)
//      static bool printOnce = true;
//      if (printOnce) {
//            printOnce = false;
//            printf("\nHTT_INPUT >>>> (printed only once)\n");
//            printf("\tYTOP=%.3f times the SM Higgs-ttbar coupling\n", httInput.YTOP);
//            if (httInput.HIGGS_OPTION==0) {
//                  printf("\tOnly SCALAR CASE, MH=%.3f GeV\n", httInput.MH);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed width: GH=%.3f GeV\n", httInput.GH);
//                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning width: GH(MH)=%.3f GeV\n", httInput.GH);
//                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/s*httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed width: GH(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", mh_gh/s*httInput.MH);
//                  }
//            } else if (httInput.HIGGS_OPTION==1) {
//                  printf("\tOnly PSEUDO-SCALAR CASE, MA=%.3f GeV\n", httInput.MA);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed width: GA=%.3f GeV\n", httInput.GA);
//                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning width: GA(MA)=%.3f GeV\n", httInput.GA);
//                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/s*httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed width: GA(top+bottom)=%.3f GeV\n", ma_ga/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", ma_ga/s*httInput.MA);
//                  }
//            } else {
//                  printf("\tSCALAR + PSEUDO-SCALAR CASE, MH=%.3f GeV, MA=%.3f GeV\n", httInput.MH, httInput.MA);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed widths: GH=%.3f GeV, GA=%.3f GeV\n", httInput.GH, httInput.GA);
//                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH, ma_ga_partial/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning widths: GH(MH)=%.3f GeV, GA(MA)=%.3f GeV\n", httInput.GH, httInput.GA);
//                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/s*httInput.MH, ma_ga_partial/s*httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH, ma_ga/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/s*httInput.MH, ma_ga/s*httInput.MA);
//                  }
//            }
//            printf("\n");
//      }
//
//      Number2 c1 = P4gen_d1_trest.CosTheta();
//      //std::cout << "c1: " << c1 << std::endl;
//      Number2 c2 = vec_4_l_trest.CosTheta();
//      Number2 s1 = sqrt(1-c1*c1);
//      Number2 s2 = sqrt(1-c2*c2);
//      Number2 phi1 = P4gen_d1_trest.Phi();
//      Number2 phi2 = vec_4_l_trest.Phi();
//     
//      // os factor deleted one sqrt in last line 
//      Number2 os_factor = (1 + z_sq + (1 - z_sq) * (1 - beta_sq_ref)) * (1. - c1*c2)
//            + 2 * (1 - beta_sq_ref) * (1 - z_sq) * c1 * c2
//            - beta_sq_ref * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
//            - 2 * (1 - beta_sq_ref) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
//            + 2 * sqrt(1 - beta_sq_ref) * z * (1-z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
//
//
//      // old os factor
//      //Number2 os_factor = (1+z_sq+(1-z_sq)*(1-beta_sq_ref))*(1.-c1*c2) + 2*(1-beta_sq_ref)*(1-z_sq)*c1*c2
//      //      - beta_sq_ref*(1-z_sq)*s1*s2*cos(phi1+phi2)
//      //      -2*(1-beta_sq_ref)*(1-z_sq)*s1*s2*cos(phi1)*cos(phi2)
//      //      +2*sqrt(1-beta_sq_ref)*z*sqrt(1-z_sq)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
//      //std::cout << "os_factor: " << os_factor << std::endl;
//      Number2 ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
//      //std::cout << "ss_factor_1: " << ss_factor_1 << std::endl;
//      Number2 ss_factor_beta = beta_sq * (1+c1*c2 - s1*s2*cos(phi1-phi2));
//      
//      Number2 additional_factor_for_nointerf = 8 * pow(m_t, 2) * pi * alpha_s * beta / pow(s, 2);
//      
//      //std::cout << "ss_factor_beta" << ss_factor_beta << std::endl;
//      //Number2 M2_QCD = common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor
//      //              + common_factor_for_M2_gg_QCD*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta);
//      
//      // just to try it out, new M for QCD, looking at the other changes I made 
//      // first part like the qcd part in BSM term and second part as it is, but (1-beta^2) becomes the same factor as in the no interference part
//      Number2 M2_QCD = common_factor_for_M2_gg_QCD * pow(alpha_s, 2) * pow(beta, 3) / (8 * pi * s) * (1 - z_sq) * os_factor
//                    + common_factor_for_M2_gg_QCD * additional_factor_for_nointerf * (ss_factor_1 + ss_factor_beta);
//
//      Number2 beta_z_factor = 1 - beta_sq_ref * z_sq;
//      Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
//      Number2 factor_interference = 4 / 3 * pow(m_t, 2) * pow(pi, 3) * pow(alpha_s, 2) / (pow(s, 2) * beta_z_factor_sq);
//      Number2 factor_in_norm_interference = sqrt_2 * pow(s, 2) * G_F * beta_z_factor / (16 * pi_sq);
//     
//      Number2 factor_interference_pseudo = 0;
//      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
//          factor_interference_pseudo = pow(std::norm(std::complex<Number2>(1., 0.) - ( factor_in_norm_interference * PB / denomA )), 2);
//      }
//      std::cout << "N = " << NB << std::endl;
//      std::cout << "denom H: " << denomH << std::endl;
//      Number2 factor_interference_scalar = 0; 
//      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
//          factor_interference_scalar = pow(std::norm(std::complex<Number2>(1., 0.) - ( factor_in_norm_interference * NB / denomH )), 2);
//      }
//      // terms to add for the BSM cross section
//      Number2 qcd_term = 0;
//      if (event.res_int == res_int::resonance or event.res_int == res_int::both) {
//          qcd_term = pow(alpha_s, 2) * pow(beta, 3) / (8 * s * pi) * common_factor_for_M2_gg_QCD * (1 - z_sq) * os_factor;
//      }
//
//      Number2 no_interference_term = 0;
//      if (res_int == res_int::resonance or res_int == res_int::both) {
//          no_interference_term = common_factor_for_M2_gg_ss_nointerf * additional_factor_for_nointerf * (ss_factor_1 + ss_factor_beta);
//      }
// 
//      // if beta_sq_ref < 0, calculate the old way, because there's no decription of what to do in that case in the paper
//      Number2 interference_term_pseudo;
//      if (beta_sq_ref>0 and (res_int == res_int::interference or res_int == res_int::both)) {
//            interference_term_pseudo = factor_interference * beta_ref * factor_interference_pseudo * ss_factor_1;
//      }
//      else if (res_int == res_int::interference or res_int == res_int::both) {
//            std::cout << "oh no! beta_sq_ref < 0, I'll use the old definition." << std::endl;
//            interference_term_pseudo = + std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_A) * ss_factor_1;
//      }
//      else {
//            interference_term_pseudo = 0;
//      }
//
//      Number2 interference_term_scalar;
//      if (beta_sq_ref>0 and (res_int == res_int::interference or res_int == res_int::both)) {
//            interference_term_scalar = factor_interference * pow(beta_ref, 3) * factor_interference_scalar * ss_factor_beta;
//      }
//      else if (res_int == res_int::interference or res_int == res_int::both) {
//            std::cout << "oh no! beta_sq_ref < 0, I'll use the old definition." << std::endl;
//            interference_term_scalar = std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_H) * ss_factor_beta;
//      }
//      else {
//            interference_term_scalar = 0;
//      }
//
//
//      std::cout << "qcd term: " << qcd_term << std::endl;
//      std::cout << "no interference term: " << no_interference_term << std::endl;
//      std::cout << "interference term pseudo scalar: " << interference_term_pseudo << std::endl;
//      std::cout << "interference term scalar: " << interference_term_scalar << std::endl;  
//
//      // BSM cross section
//      Number2 M2 = qcd_term + no_interference_term + interference_term_pseudo + interference_term_scalar;
//      
//      std::cout << "M2_gg_QCD factor: " << common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor << std::endl;
//      std::cout << "M_gg_ss_nointerf factor: " << common_factor_for_M2_gg_ss_nointerf*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta) << std::endl;
//      std::cout << "M_gg_ss_interf M_A factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1 << std::endl;
//      std::cout << "M_gg_ss_interf M_H factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta << std::endl;
//
//      // Final weight
//      Number2 weight = M2 / M2_QCD;
//      if (isnan(weight)) {
//            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
//            return 0.;
//      }
//
//     return weight;
//     // return weight - 1;
//}
//
//
//
//template <typename Number2>
//Number2 weight_ggHtt_juan_code(const HTT_Input<Number2>& httInput, res_int_t res_int) {
//      // Check that the chosen options are valid
//      if (httInput.HIGGS_OPTION>2 || httInput.WIDTH_OPTION>3) {
//            printf("What? HIGGS_OPTION=%d, WIDTH_OPTION=%d; returning weight=1 !!\n", httInput.HIGGS_OPTION, httInput.WIDTH_OPTION);
//            return 1.;
//      }
//
//      Number2 MH2_REF = pow(httInput.MH,2);
//      Number2 MA2_REF = pow(httInput.MA,2);
//      Number2 ytop2 = pow(httInput.YTOP,2);
//
//      TLorentzVector P4gen_t1 = httInput.P4gen_t[0]; 
//      TLorentzVector P4gen_t2 = httInput.P4gen_t[1]; 
//      TLorentzVector P4gen_d1 = httInput.P4gen_d[0]; 
//      TLorentzVector P4gen_d2 = httInput.P4gen_d[1]; 
//
//      const Number2 MT_REF = 172.5;
//      const Number2 MT2_REF = MT_REF*MT_REF;
//      //Fermi-Konstante
//      const Number2 GF = 1.16637876e-5;
//      const Number2 pi = TMath::Pi();
//      const Number2 pi2 = pi*pi;
//      const Number2 sqrt2 = sqrt(2.);
//
//      TLorentzVector P4gen_tt = P4gen_t1+P4gen_t2;
//      Number2 sqrts = P4gen_tt.M();
//      Number2 sqrts2 = sqrts*sqrts;
//      Number2 mt1 = httInput.P4gen_t[0].M();
//      Number2 mt2 = httInput.P4gen_t[1].M();
//      Number2 beta2 = (1-(mt1-mt2)*(mt1-mt2)/sqrts2)*(1-(mt1+mt2)*(mt1+mt2)/sqrts2);
//      //Number2 beta = sqrt(beta2);
//      //Number2 beta3 = beta*beta2;
//      Number2 mtrefOverE = 2*MT_REF/sqrts;
//      Number2 beta2Ref = 1-mtrefOverE*mtrefOverE;
//
//      TVector3 zBoost(0.,0.,-P4gen_tt.Pz()/P4gen_tt.E()); P4gen_tt.Boost(zBoost);
//      TVector3 transverseBoost(-P4gen_tt.Px()/P4gen_tt.E(),-P4gen_tt.Py()/P4gen_tt.E(),0.); P4gen_tt.Boost(transverseBoost);
//      P4gen_t1.Boost(zBoost); P4gen_t1.Boost(transverseBoost);
//      P4gen_d1.Boost(zBoost); P4gen_d1.Boost(transverseBoost);
//      P4gen_t2.Boost(zBoost); P4gen_t2.Boost(transverseBoost);
//      P4gen_d2.Boost(zBoost); P4gen_d2.Boost(transverseBoost);
//      Number2 z1 = P4gen_t1.CosTheta();
//      Number2 z2 = z1*z1;
//
//      Number2 minusPhi1Plane = -P4gen_t1.Phi();
//      TVector3 zVect(0.,0.,1.);
//      P4gen_t1.Rotate(minusPhi1Plane,zVect);
//      P4gen_d1.Rotate(minusPhi1Plane,zVect);
//      P4gen_t2.Rotate(minusPhi1Plane,zVect);
//      P4gen_d2.Rotate(minusPhi1Plane,zVect);
//      Number2 t1Yangle = -atan2(P4gen_t1.Px(),P4gen_t1.Pz());
//      TVector3 yVect(0.,1.,0.);
//      P4gen_t1.Rotate(t1Yangle,yVect);
//      P4gen_d1.Rotate(t1Yangle,yVect);
//      P4gen_t2.Rotate(t1Yangle,yVect);
//      P4gen_d2.Rotate(t1Yangle,yVect);
//
//      TVector3 t1Boost(-P4gen_t1.Px()/P4gen_t1.E(),-P4gen_t1.Py()/P4gen_t1.E(),-P4gen_t1.Pz()/P4gen_t1.E()); 
//      P4gen_d1.Boost(t1Boost);
//      TLorentzVector P4gen_d1_trest = P4gen_d1;
//      TVector3 t2Boost(-P4gen_t2.Px()/P4gen_t2.E(),-P4gen_t2.Py()/P4gen_t2.E(),-P4gen_t2.Pz()/P4gen_t2.E()); 
//      P4gen_d2.Boost(t2Boost);
//      TLorentzVector P4gen_d2_trest = P4gen_d2;
//
//
//      Number2 common_factor_for_M2_gg_QCD = 0;
//      if (res_int == res_int::resonance or res_int == res_int::both ){
//           common_factor_for_M2_gg_QCD = pi2/12*(7+9*beta2*z2)/pow(1-beta2*z2,2); 
//      }
//      Number2 common_factor_for_M2_gg_ss_nointerf = 0;
//      if (res_int == res_int::resonance or res_int == res_int::both ){
//          common_factor_for_M2_gg_ss_nointerf = pi2/12*(5+9*beta2*z2)/pow(1-beta2*z2,2);
//      }
//      Number2 common_factor_for_M_gg_ss_interf = 0;
//      if (res_int == res_int::interference or res_int == res_int::both ){
//          common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta2*z2);
//      }
//
//      using namespace std::complex_literals;
//      //const std::complex<Number2> i_cmplx = 1i;
//      //const std::complex<Number2>{0.0, static_cast<double>(d)};
//      const std::complex i_cmplx = std::complex<Number2>(0., 1.L);
//      //const std::complex<Number2> i_cmplx = static_cast<Number2>(1);
//
//      std::complex<Number2> common_factor_for_M_H = 0.;
//      Number2 mh_gh_partial = 0.;
//      Number2 mh_gh = 0.;
//      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
//            if (httInput.MH>2*MT_REF) {
//                  Number2 beta3_tdecay = pow(1-4*MT2_REF/MH2_REF,1.5);
//                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
//                        mh_gh_partial += 3*GF*MT2_REF*sqrts2/(4*pi*sqrt2)*beta3_tdecay*ytop2;
//                  } else {
//                        mh_gh_partial += 3*GF*MT2_REF*MH2_REF/(4*pi*sqrt2)*beta3_tdecay*ytop2;
//                  }
//            }
//
//            if (httInput.WIDTH_OPTION==1) {
//                  mh_gh = httInput.GH*sqrts2/httInput.MH;
//            } else if (httInput.WIDTH_OPTION==0) {
//                  mh_gh = httInput.GH*httInput.MH;
//            } else {
//                  mh_gh = mh_gh_partial;
//            }
//
//            Number2 NB_real, NB_imag;
//            Number2 auxh = 1.5 * (1-beta2Ref) * ytop2;
//            if (beta2Ref>0) {
//                  Number2 betaRef = sqrt(beta2Ref);
//                  NB_real = auxh * (1 - beta2Ref / 4 * (pow( log((1+betaRef) / (1-betaRef)), 2) - pi2 ));
//                  NB_imag = auxh*(beta2Ref/2*pi*(log((1+betaRef)/(1-betaRef))));
//            } else {
//                  NB_real = auxh * (1 + beta2Ref * pow(asin(sqrts/2/MT_REF), 2));
//                  NB_imag = 0.;
//            }
//            std::complex<Number2> NB = NB_real + NB_imag * i_cmplx;
//            std::complex<Number2> denomH = sqrts2 - MH2_REF + mh_gh * i_cmplx;
//
//            common_factor_for_M_H = NB / denomH;
//            common_factor_for_M_H *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
//      }
//
//      std::complex<Number2> common_factor_for_M_A = 0.;
//      Number2 ma_ga_partial = 0.;
//      Number2 ma_ga = 0.;
//      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
//            if (httInput.MA>2*MT_REF) {
//                  Number2 beta_tdecay = pow(1-4*MT2_REF/MA2_REF,0.5);
//                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
//                        ma_ga_partial += 3*GF*MT2_REF*sqrts2/(4*pi*sqrt2)*beta_tdecay*ytop2;
//                  } else {
//                        ma_ga_partial += 3*GF*MT2_REF*MA2_REF/(4*pi*sqrt2)*beta_tdecay*ytop2;
//                  }
//            }
//            if (httInput.WIDTH_OPTION==1) {
//                        ma_ga = httInput.GA*sqrts2/httInput.MA;
//            } else if (httInput.WIDTH_OPTION==0) {
//                        ma_ga = httInput.GA*httInput.MA;
//            } else {
//                        ma_ga = ma_ga_partial;
//            }
//
//            Number2 PB_real, PB_imag;
//            Number2 auxa = -1.5*(1-beta2Ref)/4.*ytop2;
//            if (beta2Ref>0) {
//                  Number2 betaRef = sqrt(beta2Ref);
//                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi2);
//                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
//            } else {
//                  PB_real = -4*auxa*pow(asin(sqrts/2/MT_REF),2);
//                  PB_imag = 0.;
//            }
//            std::complex<Number2> PB = PB_real + PB_imag * i_cmplx;
//            std::complex<Number2> denomA = sqrts2 - MA2_REF + ma_ga * i_cmplx;
//
//            common_factor_for_M_A = PB / denomA;
//            common_factor_for_M_A *= pow(sqrts2,1.5)*MT_REF*GF/6.*sqrt(3.)/4./pi;
//      }
//
//      // Printing global input information (once)
//      static bool printOnce = true;
//      if (printOnce) {
//            printOnce = false;
//            printf("\nHTT_INPUT >>>> (printed only once)\n");
//            printf("\tYTOP=%.3f times the SM Higgs-ttbar coupling\n", httInput.YTOP);
//            if (httInput.HIGGS_OPTION==0) {
//                  printf("\tOnly SCALAR CASE, MH=%.3f GeV\n", httInput.MH);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed width: GH=%.3f GeV\n", httInput.GH);
//                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning width: GH(MH)=%.3f GeV\n", httInput.GH);
//                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/sqrts2*httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed width: GH(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", mh_gh/sqrts2*httInput.MH);
//                  }
//            } else if (httInput.HIGGS_OPTION==1) {
//                  printf("\tOnly PSEUDO-SCALAR CASE, MA=%.3f GeV\n", httInput.MA);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed width: GA=%.3f GeV\n", httInput.GA);
//                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning width: GA(MA)=%.3f GeV\n", httInput.GA);
//                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/sqrts2*httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed width: GA(top+bottom)=%.3f GeV\n", ma_ga/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", ma_ga/sqrts2*httInput.MA);
//                  }
//            } else {
//                  printf("\tSCALAR + PSEUDO-SCALAR CASE, MH=%.3f GeV, MA=%.3f GeV\n", httInput.MH, httInput.MA);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed widths: GH=%.3f GeV, GA=%.3f GeV\n", httInput.GH, httInput.GA);
//                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH, ma_ga_partial/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning widths: GH(MH)=%.3f GeV, GA(MA)=%.3f GeV\n", httInput.GH, httInput.GA);
//                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/sqrts2*httInput.MH, ma_ga_partial/sqrts2*httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH, ma_ga/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/sqrts2*httInput.MH, ma_ga/sqrts2*httInput.MA);
//                  }
//            }
//            printf("\n");
//      }
//
//      Number2 c1 = P4gen_d1_trest.CosTheta();
//      //std::cout << "c1: " << c1 << std::endl;
//      Number2 c2 = P4gen_d2_trest.CosTheta();
//      Number2 s1 = sqrt(1-c1*c1);
//      Number2 s2 = sqrt(1-c2*c2);
//      Number2 phi1 = P4gen_d1_trest.Phi();
//      Number2 phi2 = P4gen_d2_trest.Phi();
//      Number2 os_factor = (1+z2+(1-z2)*(1-beta2Ref))*(1.-c1*c2) + 2*(1-beta2Ref)*(1-z2)*c1*c2
//            - beta2Ref*(1-z2)*s1*s2*cos(phi1+phi2)
//            -2*(1-beta2Ref)*(1-z2)*s1*s2*cos(phi1)*cos(phi2)
//            +2*sqrt(1-beta2Ref)*z1*sqrt(1-z2)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
//      //std::cout << "os_factor: " << os_factor << std::endl;
//      Number2 ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
//      //std::cout << "ss_factor_1: " << ss_factor_1 << std::endl;
//      Number2 ss_factor_beta = beta2 * (1+c1*c2 - s1*s2*cos(phi1-phi2));
//      //std::cout << "ss_factor_beta" << ss_factor_beta << std::endl;
//      Number2 M2_QCD = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
//                    + common_factor_for_M2_gg_QCD*(1-beta2Ref)*(ss_factor_1+ss_factor_beta);
//
//      //Number2 M2 = common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor
//      //          + common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
//      //          + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
//      //          + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;i
//      Number2 M2 = common_factor_for_M2_gg_QCD * beta2 * (1-z2) * os_factor 
//                  + common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta)
//                  + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1
//                  + 2*common_factor_for_M_gg_ss_interf*mtrefOverE*common_factor_for_M_A.real()*ss_factor_1
//                  + std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta;
//      std::cout << "M2_gg_QCD factor: " << common_factor_for_M2_gg_QCD*beta2*(1-z2)*os_factor << std::endl;
//      std::cout << "M_gg_ss_nointerf factor: " << common_factor_for_M2_gg_ss_nointerf*(1-beta2Ref)*(ss_factor_1+ss_factor_beta) << std::endl;
//      std::cout << "M_gg_ss_interf M_A factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1 << std::endl;
//      std::cout << "M_gg_ss_interf M_H factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta << std::endl;
//
//      // Final weight
//      Number2 weight = M2 / M2_QCD;
//      if (isnan(weight)) {
//            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
//            return 0.;
//      }
//
//     // return weight;
//     return weight - 1;
//}
//
//
//
//template <typename Number2>
//Number2 weight_ggHtt_juan_paper_different_M2(const HTT_Input<Number2>& httInput, res_int res_int) {
///*
//struct HTT_Input {
//       unsigned int HIGGS_OPTION; // 0: SCALAR only; 1: PSEUDOSCALAR only; 2: BOTH
//       double MH; // mass of scalar Higgs
//       double MA; // mass of pseudo-scalar Higgs
//       double GH; // width of scalar Higgs (only used if WIDTH_OPTION<2)
//       double GA; // width of pseudoscalar Higgs (only used if WIDTH_OPTION<2)
//       double YTOP; // times the SM Httbar coupling (1/tan(beta) in 2HDM aligned scenarios); also used to determine widths for WIDTH_OPTION>=2; note that when WIDTH_OPTION<2 this parameter is interpreted as a scale factor w.r.t the effective ggH (ggA) SM Yukawa coupling, and we assume that only the H->ttbar (A->ttbar) decay channel has sizable contributions to the total H (A) width
//       unsigned int WIDTH_OPTION; // 0: fixed widths taken from GH and GA; 1: sqrt(s) running widths taken from GH and GA; 2: determine FIXED widths from YTOP; 3: RUNNING widths determined from YTOP
//       TLorentzVector P4gen_t[2]; // index 0: top; index 1: antitop
//       TLorentzVector P4gen_d[2]; // index 0: antidown-type fermion from W+ decay; index 1: down-type fermion from W- decay
//};
//*/
//      // Check that the chosen options are valid
//      if (httInput.HIGGS_OPTION>2 || httInput.WIDTH_OPTION>3) {
//            printf("What? HIGGS_OPTION=%d, WIDTH_OPTION=%d; returning weight=1 !!\n", httInput.HIGGS_OPTION, httInput.WIDTH_OPTION);
//            return 1.;
//      }
//
//      // define all important variables
//      Number2 m_H_ref_sq = pow(httInput.MH,2);
//      Number2 m_A_ref_sq = pow(httInput.MA,2);
//      // coupling constant
//      Number2 ytop_sq = pow(httInput.YTOP,2);
//
//      TLorentzVector vec_4_t = httInput.P4gen_t[0]; 
//      TLorentzVector vec_4_tbar = httInput.P4gen_t[1]; 
//      TLorentzVector vec_4_lbar = httInput.P4gen_d[0]; 
//      TLorentzVector vec_4_l = httInput.P4gen_d[1]; 
//
//      const Number2 m_t_ref = 172.5;
//      const Number2 m_t_ref_sq = m_t_ref*m_t_ref;
//      //Fermi-Konstante
//      const Number2 G_F = 1.16637876e-5;
//      const Number2 alpha_s = 0.12;
//      const Number2 pi = TMath::Pi();
//      const Number2 pi_sq = pi*pi;
//      const Number2 sqrt_2 = sqrt(2.);
//
//      TLorentzVector vec_4_ttbar = vec_4_t+vec_4_tbar;
//      Number2 sqrt_s = vec_4_ttbar.M();
//      Number2 s = sqrt_s*sqrt_s;
//      Number2 m_t = httInput.P4gen_t[0].M();
//      Number2 m_tbar = httInput.P4gen_t[1].M();
//      Number2 beta_sq = (1-(m_t-m_tbar)*(m_t-m_tbar)/s)*(1-(m_t+m_tbar)*(m_t+m_tbar)/s);
//      Number2 beta = sqrt(beta_sq);
//      //Number2 beta3 = beta*beta_sq;
//      std::cout << "top mass: " << m_t << std::endl;
//      Number2 mtrefOverE = 2 * m_t_ref / sqrt_s;
//      std::cout << "square root of s: " << sqrt_s << std::endl;
//      std::cout << "mt over E: " << mtrefOverE << std::endl; 
//      Number2 beta_sq_ref = 1 - mtrefOverE * mtrefOverE;
//      std::cout << "beta squared ref: " << beta_sq_ref << std::endl;
//      Number2 beta_ref = 0;
//      if (beta_sq_ref > 0) {
//          beta_ref = sqrt(beta_sq_ref);
//      }
//
//      TVector3 zBoost(0.,0.,-vec_4_ttbar.Pz()/vec_4_ttbar.E()); vec_4_ttbar.Boost(zBoost);
//      TVector3 transverseBoost(-vec_4_ttbar.Px()/vec_4_ttbar.E(),-vec_4_ttbar.Py()/vec_4_ttbar.E(),0.); vec_4_ttbar.Boost(transverseBoost);
//      vec_4_t.Boost(zBoost); vec_4_t.Boost(transverseBoost);
//      vec_4_lbar.Boost(zBoost); vec_4_lbar.Boost(transverseBoost);
//      vec_4_tbar.Boost(zBoost); vec_4_tbar.Boost(transverseBoost);
//      vec_4_l.Boost(zBoost); vec_4_l.Boost(transverseBoost);
//      
//      Number2 z = vec_4_t.CosTheta();
//      Number2 z_sq = z*z;
//
//      Number2 minusPhi1Plane = -vec_4_t.Phi();
//      TVector3 zVect(0.,0.,1.);
//      vec_4_t.Rotate(minusPhi1Plane,zVect);
//      vec_4_lbar.Rotate(minusPhi1Plane,zVect);
//      vec_4_tbar.Rotate(minusPhi1Plane,zVect);
//      vec_4_l.Rotate(minusPhi1Plane,zVect);
//      Number2 t1Yangle = -atan2(vec_4_t.Px(),vec_4_t.Pz());
//      TVector3 yVect(0.,1.,0.);
//      vec_4_t.Rotate(t1Yangle,yVect);
//      vec_4_lbar.Rotate(t1Yangle,yVect);
//      vec_4_tbar.Rotate(t1Yangle,yVect);
//      vec_4_l.Rotate(t1Yangle,yVect);
//
//      TVector3 t1Boost(-vec_4_t.Px()/vec_4_t.E(),-vec_4_t.Py()/vec_4_t.E(),-vec_4_t.Pz()/vec_4_t.E()); 
//      vec_4_lbar.Boost(t1Boost);
//      TLorentzVector P4gen_d1_trest = vec_4_lbar;
//      TVector3 t2Boost(-vec_4_tbar.Px()/vec_4_tbar.E(),-vec_4_tbar.Py()/vec_4_tbar.E(),-vec_4_tbar.Pz()/vec_4_tbar.E()); 
//      vec_4_l.Boost(t2Boost);
//      TLorentzVector vec_4_l_trest = vec_4_l;
//
//
//      Number2 common_factor_for_M2_gg_QCD = pi_sq/12*(7+9*beta_sq*z_sq)/pow(1-beta_sq*z_sq,2);
//      Number2 common_factor_for_M2_gg_ss_nointerf = pi_sq/12*(5+9*beta_sq*z_sq)/pow(1-beta_sq*z_sq,2);
//      Number2 common_factor_for_M_gg_ss_interf = pi/sqrt(6.)/(1-beta_sq*z_sq);
//
//
//      using namespace std::complex_literals;
//      //const std::complex<Number2> i_cmplx = 1i;
//      //const std::complex<Number2>{0.0, static_cast<double>(d)};
//      const std::complex i_cmplx = std::complex<Number2>(0., 1.L);
//      //const std::complex<Number2> i_cmplx = static_cast<Number2>(1);
//
//
//      std::complex<Number2> NB;
//      std::complex<Number2> denomH;
//
//      std::complex<Number2> common_factor_for_M_H = 0.;
//      Number2 mh_gh_partial = 0.;
//      Number2 mh_gh = 0.;
//      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
//            if (httInput.MH>2*m_t_ref_sq) {
//                  Number2 beta3_tdecay = pow(1-4*m_t_ref_sq/m_H_ref_sq,1.5);
//                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
//                        mh_gh_partial += 3*G_F*m_t_ref_sq*s/(4*pi*sqrt_2)*beta3_tdecay*ytop_sq;
//                  } else {
//                        mh_gh_partial += 3*G_F*m_t_ref_sq*m_H_ref_sq/(4*pi*sqrt_2)*beta3_tdecay*ytop_sq;
//                  }
//            }
//
//            if (httInput.WIDTH_OPTION==1) {
//                  mh_gh = httInput.GH*s/httInput.MH;
//            } else if (httInput.WIDTH_OPTION==0) {
//                  mh_gh = httInput.GH*httInput.MH;
//            } else {
//                  mh_gh = mh_gh_partial;
//            }
//
//            Number2 NB_real, NB_imag;
//            Number2 auxh = 1.5 * (1-beta_sq_ref) * ytop_sq;
//            if (beta_sq_ref>0) {
//                  Number2 betaRef = sqrt(beta_sq_ref);
//                  NB_real = auxh * (1 - beta_sq_ref / 4 * (pow( log((1+betaRef) / (1-betaRef)), 2) - pi_sq ));
//                  NB_imag = auxh * (beta_sq_ref / 2 * pi * (log((1+betaRef)/(1-betaRef))));
//            } else {
//                  NB_real = auxh * (1 + beta_sq_ref * pow(asin(sqrt_s/2/m_t_ref_sq), 2));
//                  NB_imag = 0.;
//            }
//            NB = NB_real + NB_imag * i_cmplx;
//            denomH = s - m_H_ref_sq + mh_gh * i_cmplx;
//
//            common_factor_for_M_H = NB / denomH;
//            common_factor_for_M_H *= pow(s,1.5)*m_t_ref_sq*G_F/6.*sqrt(3.)/4./pi;
//      }
//
//      std::complex<Number2> PB;
//      std::complex<Number2> denomA;
//
//      std::complex<Number2> common_factor_for_M_A = 0.;
//      Number2 ma_ga_partial = 0.;
//      Number2 ma_ga = 0.;
//      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
//            if (httInput.MA>2*m_t_ref_sq) {
//                  Number2 beta_tdecay = pow(1-4*m_t_ref_sq/m_A_ref_sq,0.5);
//                  if (httInput.WIDTH_OPTION==1 || httInput.WIDTH_OPTION==3) {
//                        ma_ga_partial += 3*G_F*m_t_ref_sq*s/(4*pi*sqrt_2)*beta_tdecay*ytop_sq;
//                  } else {
//                        ma_ga_partial += 3*G_F*m_t_ref_sq*m_A_ref_sq/(4*pi*sqrt_2)*beta_tdecay*ytop_sq;
//                  }
//            }
//            if (httInput.WIDTH_OPTION==1) {
//                        ma_ga = httInput.GA*s/httInput.MA;
//            } else if (httInput.WIDTH_OPTION==0) {
//                        ma_ga = httInput.GA*httInput.MA;
//            } else {
//                        ma_ga = ma_ga_partial;
//            }
//
//            Number2 PB_real, PB_imag;
//            Number2 auxa = -1.5*(1-beta_sq_ref)/4.*ytop_sq;
//            if (beta_sq_ref>0) {
//                  Number2 betaRef = sqrt(beta_sq_ref);
//                  PB_real = auxa*(pow(log((1+betaRef)/(1-betaRef)),2)-pi_sq);
//                  PB_imag = auxa*(-2*pi*(log((1+betaRef)/(1-betaRef))));
//            } else {
//                  PB_real = -4*auxa*pow(asin(sqrt_s/2/m_t_ref_sq),2);
//                  PB_imag = 0.;
//            }
//            PB = PB_real + PB_imag * i_cmplx;
//            denomA = s - m_A_ref_sq + ma_ga * i_cmplx;
//
//            common_factor_for_M_A = PB / denomA;
//            common_factor_for_M_A *= pow(s,1.5)*m_t_ref_sq*G_F/6.*sqrt(3.)/4./pi;
//      }
//
//      // Printing global input information (once)
//      static bool printOnce = true;
//      if (printOnce) {
//            printOnce = false;
//            printf("\nHTT_INPUT >>>> (printed only once)\n");
//            printf("\tYTOP=%.3f times the SM Higgs-ttbar coupling\n", httInput.YTOP);
//            if (httInput.HIGGS_OPTION==0) {
//                  printf("\tOnly SCALAR CASE, MH=%.3f GeV\n", httInput.MH);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed width: GH=%.3f GeV\n", httInput.GH);
//                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning width: GH(MH)=%.3f GeV\n", httInput.GH);
//                        printf("\tSM-like width would be (ytop=%.3f): GH=%.3f GeV\n", httInput.YTOP, mh_gh_partial/s*httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed width: GH(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", mh_gh/s*httInput.MH);
//                  }
//            } else if (httInput.HIGGS_OPTION==1) {
//                  printf("\tOnly PSEUDO-SCALAR CASE, MA=%.3f GeV\n", httInput.MA);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed width: GA=%.3f GeV\n", httInput.GA);
//                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning width: GA(MA)=%.3f GeV\n", httInput.GA);
//                        printf("\tSM-like width would be (ytop=%.3f): GA=%.3f GeV\n", httInput.YTOP, ma_ga_partial/s*httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed width: GA(top+bottom)=%.3f GeV\n", ma_ga/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning width: GH(top+bottom)=%.3f GeV\n", ma_ga/s*httInput.MA);
//                  }
//            } else {
//                  printf("\tSCALAR + PSEUDO-SCALAR CASE, MH=%.3f GeV, MA=%.3f GeV\n", httInput.MH, httInput.MA);
//                  if (httInput.WIDTH_OPTION==0) {
//                        printf("\tNon-running fixed widths: GH=%.3f GeV, GA=%.3f GeV\n", httInput.GH, httInput.GA);
//                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/httInput.MH, ma_ga_partial/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==1) {
//                        printf("\tRunning widths: GH(MH)=%.3f GeV, GA(MA)=%.3f GeV\n", httInput.GH, httInput.GA);
//                        printf("\tSM-like widths would be (ytop=%.3f): GH=%.3f GeV, GA=%.3f GeV\n", httInput.YTOP, mh_gh_partial/s*httInput.MH, ma_ga_partial/s*httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==2) {
//                        printf("\tNon-running fixed widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/httInput.MH, ma_ga/httInput.MA);
//                  } else if (httInput.WIDTH_OPTION==3) {
//                        printf("\tRunning widths: GH(top+bottom)=%.3f GeV, GA(top+bottom)=%.3f GeV\n", mh_gh/s*httInput.MH, ma_ga/s*httInput.MA);
//                  }
//            }
//            printf("\n");
//      }
//
//      Number2 c1 = P4gen_d1_trest.CosTheta();
//      //std::cout << "c1: " << c1 << std::endl;
//      Number2 c2 = vec_4_l_trest.CosTheta();
//      Number2 s1 = sqrt(1-c1*c1);
//      Number2 s2 = sqrt(1-c2*c2);
//      Number2 phi1 = P4gen_d1_trest.Phi();
//      Number2 phi2 = vec_4_l_trest.Phi();
//     
//      // os factor deleted one sqrt in last line 
//      Number2 os_factor = (1 + z_sq + (1 - z_sq) * (1 - beta_sq_ref)) * (1. - c1*c2)
//            + 2 * (1 - beta_sq_ref) * (1 - z_sq) * c1 * c2
//            - beta_sq_ref * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
//            - 2 * (1 - beta_sq_ref) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
//            + 2 * sqrt(1 - beta_sq_ref) * z * (1-z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
//
//
//      // old os factor
//      //Number2 os_factor = (1+z_sq+(1-z_sq)*(1-beta_sq_ref))*(1.-c1*c2) + 2*(1-beta_sq_ref)*(1-z_sq)*c1*c2
//      //      - beta_sq_ref*(1-z_sq)*s1*s2*cos(phi1+phi2)
//      //      -2*(1-beta_sq_ref)*(1-z_sq)*s1*s2*cos(phi1)*cos(phi2)
//      //      +2*sqrt(1-beta_sq_ref)*z*sqrt(1-z_sq)*(c1*s2*cos(phi2)+c2*s1*cos(phi1));
//      //std::cout << "os_factor: " << os_factor << std::endl;
//      Number2 ss_factor_1 = 1+c1*c2 + s1*s2*cos(phi1-phi2);
//      //std::cout << "ss_factor_1: " << ss_factor_1 << std::endl;
//      Number2 ss_factor_beta = beta_sq * (1+c1*c2 - s1*s2*cos(phi1-phi2));
//      
//      Number2 additional_factor_for_nointerf = 8 * pow(m_t, 2) * pi * alpha_s * beta / pow(s, 2);
//      
//      //std::cout << "ss_factor_beta" << ss_factor_beta << std::endl;
//      //Number2 M2_QCD = common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor
//      //              + common_factor_for_M2_gg_QCD*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta);
//      
//      // just to try it out, new M for QCD, looking at the other changes I made 
//      // first part like the qcd part in BSM term and second part as it is, but (1-beta^2) becomes the same factor as in the no interference part
//      Number2 M2_QCD = common_factor_for_M2_gg_QCD * pow(alpha_s, 2) * pow(beta, 3) / (8 * pi * s) * (1 - z_sq) * os_factor
//                    + common_factor_for_M2_gg_QCD * additional_factor_for_nointerf * (ss_factor_1 + ss_factor_beta);
//
//
//      Number2 beta_z_factor = 1 - beta_sq_ref * z_sq;
//      Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
//      Number2 factor_interference = 4 / 3 * pow(m_t, 2) * pow(pi, 3) * pow(alpha_s, 2) / (pow(s, 2) * beta_z_factor_sq);
//      Number2 factor_in_norm_interference = sqrt_2 * pow(s, 2) * G_F * beta_z_factor / (16 * pi_sq);
//     
//      Number2 factor_interference_pseudo = 0;
//      if (httInput.HIGGS_OPTION==1 || httInput.HIGGS_OPTION==2) {
//          factor_interference_pseudo = pow(std::norm(std::complex<Number2>(1., 0.) - ( factor_in_norm_interference * PB / denomA )), 2);
//      }
//      std::cout << "N = " << NB << std::endl;
//      std::cout << "denom H: " << denomH << std::endl;
//      Number2 factor_interference_scalar = 0; 
//      if (httInput.HIGGS_OPTION==0 || httInput.HIGGS_OPTION==2) {
//          factor_interference_scalar = pow(std::norm(std::complex<Number2>(1., 0.) - ( factor_in_norm_interference * NB / denomH )), 2);
//      }
//      // terms to add for the BSM cross section
//      Number2 qcd_term = 0; 
//      if (res_int == res_int::interference or res_int == res_int::both) {
//           qcd_term = pow(alpha_s, 2) * pow(beta, 3) / (8 * s * pi) * common_factor_for_M2_gg_QCD * (1 - z_sq) * os_factor;
//      }
//
//      Number2 no_interference_term = 0;
//      if (res_int == res_int::resonance or res_int == res_int::both) {
//          no_interference_term = common_factor_for_M2_gg_ss_nointerf * additional_factor_for_nointerf * (ss_factor_1 + ss_factor_beta);
//      }
// 
//      // if beta_sq_ref < 0, calculate the old way, because there's no decription of what to do in that case in the paper
//      Number2 interference_term_pseudo;
//      if (beta_sq_ref>0 and (res_int == res_int::interference or res_int == res_int::both)) {
//            interference_term_pseudo = factor_interference * beta_ref * factor_interference_pseudo * ss_factor_1;
//      }
//      else if (res_int == res_int::interference or res_int == res_int::both) {
//            std::cout << "oh no! beta_sq_ref < 0, I'll use the old definition." << std::endl;
//            interference_term_pseudo = + std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_A) * ss_factor_1;
//      }
//      else {
//            interference_term_pseudo = 0;
//      }
//
//      Number2 interference_term_scalar;
//      if (beta_sq_ref>0 and (res_int == res_int::interference or res_int == res_int::both)) {
//            interference_term_scalar = factor_interference * pow(beta_ref, 3) * factor_interference_scalar * ss_factor_beta;
//      }
//      else if (res_int == res_int::interference or res_int == res_int::both) {
//            std::cout << "oh no! beta_sq_ref < 0, I'll use the old definition." << std::endl;
//            interference_term_scalar = std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_H) * ss_factor_beta;
//      }
//      else {
//            interference_term_scalar = 0;
//      }
//
//      std::cout << "qcd term: " << qcd_term << std::endl;
//      std::cout << "no interference term: " << no_interference_term << std::endl;
//      std::cout << "interference term pseudo scalar: " << interference_term_pseudo << std::endl;
//      std::cout << "interference term scalar: " << interference_term_scalar << std::endl;  
//
//      // BSM cross section
//      //Number2 M2 = qcd_term + no_interference_term + interference_term_pseudo + interference_term_scalar;
//      
//      Number2 M2 = std::norm(common_factor_for_M_gg_ss_interf * mtrefOverE - common_factor_for_M_A) * ss_factor_1;
//      
//
//      std::cout << "M2_gg_QCD factor: " << common_factor_for_M2_gg_QCD*beta_sq*(1-z_sq)*os_factor << std::endl;
//      std::cout << "M_gg_ss_nointerf factor: " << common_factor_for_M2_gg_ss_nointerf*(1-beta_sq_ref)*(ss_factor_1+ss_factor_beta) << std::endl;
//      std::cout << "M_gg_ss_interf M_A factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_A)*ss_factor_1 << std::endl;
//      std::cout << "M_gg_ss_interf M_H factor: " << std::norm(common_factor_for_M_gg_ss_interf*mtrefOverE-common_factor_for_M_H)*ss_factor_beta << std::endl;
//
//      // Final weight
//      Number2 weight = M2 / M2_QCD;
//      if (isnan(weight)) {
//            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M2_QCD, M2);
//            return 0.;
//      }
//
//     return weight;
//     //return weight - 1;
//}




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
                       Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m,
                       calc_weight_version version, higgs_type_t higgs_type, Number2 mass, Number2 width, res_int_t res_int)
{
  // what's that?
  //double pz_hard_radiation = 0.;

 // bool is_gg = is_gg_initial_state(pdg1, pdg2, pz_hard_radiation); 
	


  //HTT_Input<Number2> httInput;
  //if (higgs_type == higgs_type::scalar){
  //  httInput.HIGGS_OPTION = 0;
  //}
  //else if (higgs_type == higgs_type::pseudo_scalar){
  //  httInput.HIGGS_OPTION = 1;
  //}
  //else {
  //  std::cout << "no valid higgs type" << std::endl;
  //}

  //httInput.WIDTH_OPTION = 0;
  //httInput.YTOP = 1.0;
  //httInput.MH = mass;
  //httInput.GH = width;
  //httInput.MA = mass;
  //httInput.GA = width;

  // std::cout << "pTop_pt: " << pTop_pt << ", pTop_eta: " << pTop_eta << ", pTop_phi: " << pTop_phi << ", pTop_m: " << pTop_m << std::endl;
  // std::cout << "aTop_pt: " << aTop_pt << ", aTop_eta: " << aTop_eta << ", aTop_phi: " << aTop_phi << ", aTop_m: " << aTop_m << std::endl;
  // std::cout << "pLep_pt: " << pLep_pt << ", pLep_eta: " << pLep_eta << ", pLep_phi: " << pLep_phi << ", pLep_m: " << pLep_m << std::endl;
  // std::cout << "aLep_pt: " << aLep_pt << ", aLep_eta: " << aLep_eta << ", aLep_phi: " << aLep_phi << ", aLep_m: " << aLep_m << std::endl;

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
 
  Number2 event_weight = 1;
  
  event_weight = weight_ggHtt(event);

  //if (version==calc_weight_version::juan_paper) {
  //  event_weight = weight_ggHtt_juan_paper(httInput, res_int);
  //}
  //else if (version==calc_weight_version::juan_code){
  //  event_weight = weight_ggHtt_juan_code(httInput, res_int);
  //}
  //else if (version==calc_weight_version::juan_paper_different_M2){
  //  event_weight = weight_ggHtt_juan_paper_different_M2(httInput, res_int);
  //}
  //else {
  //  std::cout << "no valid version for weight calculation given!" << std::endl;
  //}

  //TLorentzVector tops = httInput.P4gen_t[0] + httInput.P4gen_t[1];
  
  return event_weight;

}




template <typename Number1, typename Number2>
auto calculate_weight(calc_weight_version version, higgs_type_t higgs_type, float mass, float width, res_int_t res_int)
{

  return [version, higgs_type, mass, width, res_int] (Number1 pTop_pt, Number1 pTop_eta, Number1 pTop_phi, Number1 pTop_m,
             Number1 aTop_pt, Number1 aTop_eta, Number1 aTop_phi, Number1 aTop_m,
             Number1 pLep_pt, Number1 pLep_eta, Number1 pLep_phi, Number1 pLep_m,
             Number1 aLep_pt, Number1 aLep_eta, Number1 aLep_phi, Number1 aLep_m)
           //            Number ini1_pt, Number ini1_eta, Number ini1_phi, Number ini1_m,
           //            Number ini2_pt, Number ini2_eta, Number ini2_phi, Number ini2_m)
  {

       return compute_weight<Number1, Number2>(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                          aTop_pt, aTop_eta, aTop_phi, aTop_m,
                          pLep_pt, pLep_eta, pLep_phi, pLep_m,
                          aLep_pt, aLep_eta, aLep_phi, aLep_m, 
                          version, higgs_type, mass, width, res_int);
  };
}
//}



