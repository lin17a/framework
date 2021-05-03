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
        this->s = s;
        Number2 m_t = vec_4_t.M();
        this->m_t = m_t;
        Number2 m_tbar = vec_4_tbar.M();
        this->m_tbar = m_tbar;
        Number2 beta_sq = (1 - (m_t - m_tbar) * (m_t - m_tbar) / s) * (1 - (m_t + m_tbar) * (m_t + m_tbar) / s);
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
        //std::cout << "square root of s: " << sqrt_s << std::endl;
        //std::cout << "mt over E: " << mtrefOverE << std::endl; 
        Number2 beta_ref_sq = 1 - mtrefOverE * mtrefOverE;
        this->beta_ref_sq = beta_ref_sq;
        //std::cout << "beta squared ref: " << beta_ref_sq << std::endl;
        Number2 beta_ref = 0;
        if (beta_ref_sq > 0) {
            beta_ref = sqrt(beta_ref_sq);
        }
        this->beta_ref = beta_ref;


        static const TVector3 zBase(0., 0., 1.);
        TLorentzVector vec_4_t_hel = f_zmf_tt(vec_4_t, vec_4_ttbar);
        TLorentzVector vec_4_tbar_hel = f_zmf_tt(vec_4_tbar, vec_4_ttbar);
        TLorentzVector vec_4_l_hel = f_zmf_tt(vec_4_l, vec_4_ttbar);
        TLorentzVector vec_4_lbar_hel = f_zmf_tt(vec_4_lbar, vec_4_ttbar);
        //Number2 z = vec_4_t_hel.Vect().Unit().Dot( vec_4_ttbar.Vect().Unit() ); 
        Number2 z = vec_4_t_hel.Vect().Unit().Dot(zBase);

//        TVector3 zBoost(0.,0.,-vec_4_ttbar.Pz()/vec_4_ttbar.E()); vec_4_ttbar.Boost(zBoost);
//        TVector3 transverseBoost(-vec_4_ttbar.Px()/vec_4_ttbar.E(),-vec_4_ttbar.Py()/vec_4_ttbar.E(),0.); vec_4_ttbar.Boost(transverseBoost);
//        vec_4_t.Boost(zBoost); vec_4_t.Boost(transverseBoost);
//        vec_4_lbar.Boost(zBoost); vec_4_lbar.Boost(transverseBoost);
//        vec_4_tbar.Boost(zBoost); vec_4_tbar.Boost(transverseBoost);
//        vec_4_l.Boost(zBoost); vec_4_l.Boost(transverseBoost);
//        
//        Number2 z = vec_4_t.CosTheta();
        this->z = z;
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
        Number2 s2 = sqrt(1-c2*c2);
//        Number2 phi1 = vec_4_lbar_trest.Phi();
//        Number2 phi2 = vec_4_l_trest.Phi();
     
        // os factor deleted one sqrt in last line 
        Number2 os_factor_juan_paper = (1 + z_sq + (1 - z_sq) * (1 - beta_ref_sq)) * (1. - c1*c2)
              + 2 * (1 - beta_ref_sq) * (1 - z_sq) * c1 * c2
              - beta_ref_sq * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
              - 2 * (1 - beta_ref_sq) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
              + 4 * sqrt(1 - beta_ref_sq) * z * (1 - z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
        this->os_factor_juan_paper = os_factor_juan_paper;

        // old os factor
        Number2 os_factor_juan_code = (1 + z_sq + (1 - z_sq) * (1 - beta_ref_sq)) * (1. - c1 * c2) 
              + 2 * (1 - beta_ref_sq) * (1 - z_sq) * c1 * c2
              - beta_ref_sq * (1 - z_sq) * s1 * s2 * cos(phi1 + phi2)
              - 2 * (1 - beta_ref_sq) * (1 - z_sq) * s1 * s2 * cos(phi1) * cos(phi2)
              + 2 * sqrt(1 - beta_ref_sq) * z * sqrt(1 - z_sq) * (c1 * s2 * cos(phi2) + c2 * s1 * cos(phi1));
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
// template <typename Number2>
// Number2 weight_ggHtt(event_t<Number2> event);

// Simple code to decide whether the initial state should be gg or not
// pdgId1 is the the PDG id of first parton, along the +Z direction
// pdgId2 is the the PDG id of other parton, along the -Z direction
// pz_hard_radiation is the pz of the system made of the hard radiated partons recoiling 
// template <typename Number2>
// bool is_gg_initial_state(const int& pdgId1, const int& pdgId2, const Number2& pz_hard_radiation);

// ---- N and P ----
template<typename Number2>
std::complex<Number2> calc_N(event_t<Number2> event){

     std::complex<Number2> N;

     if (event.beta_ref_sq>0) {
           Number2 beta_ref = sqrt(event.beta_ref_sq);
           N = (Number2) 3. / (Number2) 8. * (Number2) (1 - event.beta_ref_sq) * (Number2) pow(event.y_top, 2) * ( (Number2) 4. - (std::complex<Number2>) event.beta_ref_sq * (std::complex<Number2>) pow(log(((1 + beta_ref) / (1 - beta_ref))) - constants<Number2>::i_cmplx * (Number2) TMath::Pi(), 2) );
     } else {
           Number2 auxh = 1.5 * (1 - event.beta_ref_sq) * pow(event.y_top, 2);
           Number2 NB_real = auxh * (1 + event.beta_ref_sq * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2));
           Number2 NB_imag = 0.;
           N = NB_real + NB_imag * constants<Number2>::i_cmplx;
     }
     return N;
}


template<typename Number2>
std::complex<Number2> calc_P(event_t<Number2> event){

      std::complex<Number2> P;
      if (event.beta_ref_sq>0) {
            Number2 beta_ref = sqrt(event.beta_ref_sq);
            P = - (Number2) pow(event.y_top, 2) * (Number2) 3. / (Number2) 8. * (Number2) (1 - event.beta_ref_sq) * pow( (Number2) log( (1 + beta_ref) / ( 1 - beta_ref)) - constants<Number2>::i_cmplx * (Number2) TMath::Pi(), 2);
      } else {
            Number2 auxa = -1.5 * (1 - event.beta_ref_sq) / 4. * pow(event.y_top, 2);
            Number2 PB_real = -4 * auxa * pow(asin(event.sqrt_s / 2 / constants<Number2>::m_t_ref_sq), 2);
            Number2 PB_imag = 0.;
            P = PB_real + PB_imag * constants<Number2>::i_cmplx;
      }
      return P;
}


// ---- methods without decay ----
template<typename Number2>
Number2 calc_resonance_scalar_without_decay(event_t<Number2> event){
     Number2 common_factor = pow(constants<Number2>::G_F, 2) * pow(constants<Number2>::m_t_ref, 2) * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3) );
     
     Number2 mh_gh = event.higgs_width * event.higgs_mass;
     std::complex<Number2> denomH = event.s - (Number2) pow(event.higgs_mass, 2) + mh_gh * constants<Number2>::i_cmplx;
     std::complex<Number2> N = calc_N(event);

     Number2 resonance = common_factor * pow(event.beta, 3) * std::norm( N / denomH ); 
     return resonance;   
}

template<typename Number2>
Number2 calc_resonance_pseudo_without_decay(event_t<Number2> event){
     Number2 common_factor = pow(constants<Number2>::G_F, 2) * pow(constants<Number2>::m_t_ref, 2) * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3) );
     
     Number2 ma_ga = event.higgs_width * event.higgs_mass;
     std::complex<Number2> denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;
     std::complex<Number2> P = calc_P(event);

     Number2 resonance = common_factor * event.beta * std::norm( P / denomA ); 
     return resonance;   
}


template<typename Number2>
Number2 calc_QCD_without_decay(event_t<Number2> event){

     Number2 factor_1 = TMath::Pi() * event.beta / ( 96 * event.s);
     Number2 factor_2 = (7 + event.beta_sq * event.z_sq) / pow( 1 - event.beta_sq * event.z_sq, 2 ); 
     Number2 factor_3 = event.beta_sq * (1 - pow(event.z , 4)) + 4 * pow(constants<Number2>::m_t_ref, 2) / event.s * ( 1 + event.beta_sq * pow(event.z, 4) + 2 * event.beta_sq - 2 * event.beta_sq * event.z_sq);
     return factor_1 * factor_2 * factor_3;
}


// ---- old methods ----
template<typename Number2>
Number2 calc_qcd_interference_scalar(event_t<Number2> event){

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
     switch (event.version){
        
          case calc_weight_version::juan_code:{
               qcd_interference_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H)
                                         * event.ss_factor_beta;
               break;}
          case calc_weight_version::juan_paper:{
              Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
              Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
              Number2 factor_interference = 4 / 3 * pow(constants<Number2>::m_t_ref, 2) * pow(TMath::Pi(), 3) 
                                            / (pow(event.s, 2) * beta_z_factor_sq);
              //Number2 factor_in_norm_interference = sqrt(2.) * pow(event.s, 2) * constants<Number2>::G_F * beta_z_factor / (16 * pow(TMath::Pi(), 2));
              // return zero if higgs is not scalar 
              // if (event.higgs_type != higgs_type_t::scalar) {
              //     return 0;
              // } 
              qcd_interference_scalar = factor_interference * event.beta
                                       // setting N to zero
                                       // * std::norm((Number2) 1. - factor_in_norm_interference * (NB / denomH))
                                       * event.ss_factor_beta;
              break;}
          default:
              break;
     }
     return qcd_interference_scalar;
}

template <typename Number2>
Number2 calc_qcd_interference_pseudo(event_t<Number2> event){

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
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;

      common_factor_for_M_A = PB / denomA;
      common_factor_for_M_A *= pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 qcd_interference_pseudo = 0;
 
      switch (event.version){
         case calc_weight_version::juan_code:{
             qcd_interference_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
                                          * event.ss_factor_1;
             break;}
         case calc_weight_version::juan_paper:{
             Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
             Number2 beta_z_factor_sq = beta_z_factor * beta_z_factor;
             Number2 factor_interference = 4 / 3 * pow(constants<Number2>::m_t_ref, 2) * pow(TMath::Pi(), 3)  
                                           / (pow(event.s, 2) * beta_z_factor_sq);
             //Number2 factor_in_norm_interference = sqrt(2.) * pow(event.s, 2) * constants<Number2>::G_F * beta_z_factor / (16 * pow(TMath::Pi(), 2));
             // return zero if higgs_type is not pseudo
             //if (event.higgs_type != higgs_type_t::pseudo_scalar){
             //    return 0;
             //}         
             qcd_interference_pseudo = factor_interference * event.beta
 
                                     // setting P to zero??
                                     //  * std::norm((Number2) 1. - factor_in_norm_interference * (PB / denomA))
                                      * event.ss_factor_1;
             break;}
         default:
             break;
      }
      return qcd_interference_pseudo;
}

template <typename Number2>
Number2 calc_qcd_no_interference(event_t<Number2> event){
    Number2 common_factor_for_M2_gg_ss_nointerf = pow(TMath::Pi(), 2) / 12 * (5 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    Number2 additional_factor_for_nointerf = 0;
    switch (event.version){ 
        case calc_weight_version::juan_paper:
            additional_factor_for_nointerf = 8 * pow(constants<Number2>::m_t_ref, 2) * TMath::Pi() * event.beta / pow(event.s, 2);
            break;
        case calc_weight_version::juan_code: 
            additional_factor_for_nointerf = 1 - event.beta_sq;
            break;
        default:
            break;
    }
    Number2 no_interference_term = common_factor_for_M2_gg_ss_nointerf * additional_factor_for_nointerf * (event.ss_factor_1 + event.ss_factor_beta);
    return no_interference_term;
}


template <typename Number2>
Number2 calc_qcd_opp_gluon(event_t<Number2> event){
    Number2 common_factor_for_M2_gg_QCD = pow(TMath::Pi(), 2) / 12 * (7 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    Number2 opp_gluon_qcd_term = 0;
 
    if (event.version == calc_weight_version::juan_paper){
       opp_gluon_qcd_term = 2 * pow(event.beta, 3) * TMath::Pi() / event.s * common_factor_for_M2_gg_QCD * (1 - event.z_sq) 
                  * event.os_factor_juan_paper;
    }
    if (event.version == calc_weight_version::juan_code){
            opp_gluon_qcd_term = common_factor_for_M2_gg_QCD * event.beta_sq * (1 - event.z_sq) * event.os_factor_juan_code; 
    }
    return opp_gluon_qcd_term;
}

template <typename Number2>
Number2 calc_M_2_qcd(event_t<Number2> event){
    Number2 M2_QCD = 0;
    Number2 common_factor_for_M2_gg_QCD = pow(TMath::Pi(), 2) / 12 * (7 + 9 * event.beta_sq * event.z_sq) / pow(1 - event.beta_sq * event.z_sq, 2);
    //Number2 additional_factor_for_nointerf = 8 * pow(constants<Number2>::m_t_ref, 2) * TMath::Pi() * event.beta_ref / pow(event.s, 2);
    
    if (event.version == calc_weight_version::juan_code){
        M2_QCD = common_factor_for_M2_gg_QCD * event.beta_sq * (1 - event.z_sq) * event.os_factor_juan_code \
                 + common_factor_for_M2_gg_QCD * (1 - event.beta_ref_sq) * (event.ss_factor_1 + event.ss_factor_beta);
    }
    if (event.version ==  calc_weight_version::juan_paper){
        // just to try it out, new M for QCD, looking at the other changes I made 
        // first part like the qcd part in BSM term and second part as it is, but (1-beta^2) becomes the same factor as in the no interference part
      
        std::cout << "beta_ref: " << event.beta_ref << std::endl;
        std::cout << "z_sq: " << event.z_sq << std::endl;
        std::cout << "s: " << event.s << std::endl;

        Number2 opp_gluon = calc_qcd_opp_gluon(event);
        Number2 same_gluon_no_interf = calc_qcd_no_interference(event);
        Number2 same_gluon_interf_pseudo = calc_qcd_interference_pseudo(event);
        Number2 same_gluon_interf_scalar = calc_qcd_interference_scalar(event);
 
        M2_QCD = opp_gluon + same_gluon_no_interf + same_gluon_interf_pseudo + same_gluon_interf_scalar;
        
        //M2_QCD = common_factor_for_M2_gg_QCD * pow(event.beta_ref, 3) / (8 * TMath::Pi() * event.s) * (1 - event.z_sq)
        //         * event.os_factor_juan_paper 
        //         + common_factor_for_M2_gg_QCD * additional_factor_for_nointerf * (event.ss_factor_1 + event.ss_factor_beta);
     }   
     return M2_QCD;
}


template <typename Number2>
Number2 calc_resonance_pseudo_new(event_t<Number2> event){

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
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;

      common_factor_for_M_A = PB / denomA;
      common_factor_for_M_A *= pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 resonance_pseudo = 0;
 
      switch (event.version){
         case calc_weight_version::juan_code:{
             // not fixed yet, still res and int or so?
             resonance_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
                                          * event.ss_factor_1;
             break;}
         case calc_weight_version::juan_paper:{
             Number2 factor_resonance = pow(constants<Number2>::G_F, 2) * constants<Number2>::m_t_ref_sq * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3));
             resonance_pseudo = factor_resonance * event.beta * std::norm(PB / denomA) * event.ss_factor_1;
             break;}
         default:
             break;
      }
      return resonance_pseudo;
}

template <typename Number2>
Number2 calc_interference_pseudo_new(event_t<Number2> event){

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
      denomA = event.s - (Number2) pow(event.higgs_mass, 2) + ma_ga * constants<Number2>::i_cmplx;

      common_factor_for_M_A = PB / denomA;
      common_factor_for_M_A *= pow(event.s,1.5) * constants<Number2>::m_t_ref * constants<Number2>::G_F / 6. * sqrt(3.) / 4. / TMath::Pi();
      
      Number2 interference_pseudo = 0;
 
      switch (event.version){
         case calc_weight_version::juan_code:{
             // not fixed yet, still res and int or so?
             interference_pseudo = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_A)
                                          * event.ss_factor_1;
             break;}
         case calc_weight_version::juan_paper:{
             Number2 factor_interference = constants<Number2>::G_F * constants<Number2>::m_t_ref_sq / ( 48 * sqrt(2) * TMath::Pi());
             Number2 beta_z_factor = 1 - event.beta_sq * event.z_sq;
             interference_pseudo = factor_interference * event.beta / beta_z_factor * (PB / denomA).real() * event.ss_factor_1;
             break;}
         default:
             break;
      }
      return interference_pseudo;
}

template<typename Number2>
Number2 calc_resonance_scalar_new(event_t<Number2> event){

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
     
     Number2 resonance_scalar = 0;
     switch (event.version){
          // not fixed yet
          case calc_weight_version::juan_code:{
               resonance_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H) \
                                         * event.ss_factor_beta;
               break;}
          case calc_weight_version::juan_paper:{
              Number2 factor_resonance = pow(constants<Number2>::G_F, 2) * constants<Number2>::m_t_ref_sq * pow(event.s, 2) / ( 1536 * pow(TMath::Pi(), 3)); 
              resonance_scalar = factor_resonance * event.beta * std::norm(NB / denomH) * event.ss_factor_beta;
              break;}
          default:
              break;
     }
     return resonance_scalar;
}

template<typename Number2>
Number2 calc_interference_scalar_new(event_t<Number2> event){

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
     
     Number2 interference_scalar = 0;
     switch (event.version){
          // not fixed yet
          case calc_weight_version::juan_code:{
               interference_scalar = std::norm(common_factor_for_M_gg_ss_interf * 2 * constants<Number2>::m_t_ref / event.sqrt_s - common_factor_for_M_H) \
                                         * event.ss_factor_beta;
               break;}
          case calc_weight_version::juan_paper:{
              Number2 factor_interference = constants<Number2>::G_F * constants<Number2>::m_t_ref_sq / (48 * sqrt(2.) * TMath::Pi()); 
              Number2 beta_z_factor = 1. - event.beta_sq * event.z_sq;
              interference_scalar = factor_interference * event.beta / beta_z_factor * (NB / denomH).real() * event.ss_factor_beta;
              break;}
          default:
              break;
     }
     return interference_scalar;
}

template <typename Number2>
Number2 weight_ggHtt(event_t<Number2> event){


      Number2 resonance_scalar = 0;
      Number2 interference_scalar = 0;
      Number2 resonance_pseudo = 0;
      Number2 interference_pseudo = 0;
      if (event.higgs_type == higgs_type_t::scalar){
          if (event.res_int == res_int_t::resonance or event.res_int == res_int_t::both){
              //resonance_scalar = calc_resonance_scalar_new(event);
              resonance_scalar = calc_resonance_scalar_without_decay(event);
          }
          if (event.res_int == res_int_t::interference or event.res_int == res_int_t::both){
              interference_scalar = calc_interference_scalar_new(event);
          }
      }
      else if (event.higgs_type == higgs_type_t::pseudo_scalar) {
          if (event.res_int == res_int_t::resonance or event.res_int == res_int_t::both){
              //resonance_pseudo = calc_resonance_pseudo_new(event);
              resonance_pseudo = calc_resonance_pseudo_without_decay(event);
          }
          if (event.res_int == res_int_t::interference or event.res_int == res_int_t::both){
              interference_pseudo = calc_interference_pseudo_new(event);
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
      Number2 M_2_QCD = calc_QCD_without_decay(event);
      Number2 M_2_bsm = resonance_scalar + interference_scalar + resonance_pseudo + interference_pseudo;
 
      std::cout << "M_2_qcd: " << M_2_QCD << std::endl;  
      std::cout << "M_2_bsm: " << M_2_bsm << std::endl;  

      Number2 weight = M_2_bsm / M_2_QCD;
      //Number2 weight = 1. / M_2_QCD;

      std::cout << "weight: " << weight << std::endl;  
      
      if (isnan(weight) or isinf(weight)) {
            printf("What?? weight: %.3e, M2QCD: %.3e, M2: %.3e; returning weight=0 !!\n", weight, M_2_QCD, M_2_bsm);
            return 0.;
      }
      return weight;
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

  event_t<Number2> event =  event_t(vec_4_t, vec_4_tbar, vec_4_l, vec_4_lbar, higgs_type, mass, width, res_int, version);
 
  Number2 event_weight = 1;
  
  event_weight = weight_ggHtt(event);

  
  return event_weight;

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



