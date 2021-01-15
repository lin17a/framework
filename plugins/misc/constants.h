// -*- C++ -*-
// author: afiq anuar
// short: some constants that might be needed in analysis
// note: masses are in GeV, most are taken from PDG 2020

#ifndef FWK_CONSTANTS_H
#define FWK_CONSTANTS_H

namespace Framework {
  namespace constants {
    template <typename T = float>
    constexpr T sqrt2      = T(1.4142135623730950488L);

    template <typename T = float>
    constexpr T ln2        = T(6.9314718055994530941L);

    template <typename T = float>
    constexpr T pi         = T(3.1415926535897932384L);

    template <typename T = float>
    constexpr T e          = T(2.7182818284590452353L);

    template <typename T = float>
    constexpr T m_proton   = T(0.93827208816L);

    template <typename T = float>
    constexpr T m_neutron  = T(0.93956542052L);

    template <typename T = float>
    constexpr T m_electron = T(0.00051099895L);

    template <typename T = float>
    constexpr T m_muon     = T(0.1056583745L);

    template <typename T = float>
    constexpr T m_tau      = T(1.77686L);

    // quarks are inherently problematic, so use with care
    template <typename T = float>
    constexpr T m_bottom   = T(4.18L);

    template <typename T = float>
    constexpr T m_top      = T(172.76L);

    template <typename T = float>
    constexpr T m_w        = T(80.379L);

    template <typename T = float>
    constexpr T m_z        = T(91.1876L);

    template <typename T = float>
    constexpr T m_h        = T(125.10L);

    template <typename T = float>
    constexpr T br_wlv     = T(0.1086L);

    template <typename T = float>
    constexpr T br_wqq     = T(0.6741L);

    template <typename T = float>
    constexpr T br_zll     = T(0.033658L);

    template <typename T = float>
    constexpr T br_zqq     = T(0.69911L);
  }
}

#endif
