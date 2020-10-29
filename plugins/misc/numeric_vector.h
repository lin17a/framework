#ifndef FWK_NUMERIC_VECTOR_H
#define FWK_NUMERIC_VECTOR_H

// -*- C++ -*-
// author: afiq anuar
// short: a listing of free functions that contains common numerical/vector operations

#include <numeric>
#include <cmath>
#include "TLorentzVector.h"

/// poor man's std::hypot (tested with g++/clang++, fine for double: 1e-158 < i < 1e10)
template <typename Number = float>
Number quadratic_sum(Number num1, Number num2)
{
  return std::sqrt((num1 * num1) + (num2 * num2));
}



template <typename Number = float>
Number absolute_difference(Number num1, Number num2)
{
  return std::abs(num1 - num2);
}



/// actually abs(dPhi) i.e. defined within [0, pi]
template <typename Number = float>
Number dPhi(Number phi1, Number phi2) 
{
  return std::acos( std::cos(phi1 - phi2) );
}



template <typename Number = float>
Number dR(Number eta1, Number phi1, Number eta2, Number phi2) 
{
  return quadratic_sum(eta1 - eta2, dPhi(phi1, phi2));
}



template <typename Number = float>
Number angle(Number eta1, Number phi1, Number eta2, Number phi2) 
{
  Number xdot = std::cos(phi1) * std::cos(phi2);
  Number ydot = std::sin(phi1) * std::sin(phi2);
  Number zdot = std::sinh(eta1) * std::sinh(eta2);

  return (xdot + ydot + zdot) / (std::cosh(eta1) * std::cosh(eta2));
}



/// helper for multivector_system_impl, doing the actual daughter assignment
template <size_t ...N, typename Number>
void assign_vector(std::array<TLorentzVector, sizeof...(N)> &p4, std::array<Number, 4 * sizeof...(N)> &arg, std::index_sequence<N...>)
{
  (std::get<N>(p4).SetPtEtaPhiM(std::get<4 * N>(arg), std::get<(4 * N) + 1>(arg), std::get<(4 * N) + 2>(arg), std::get<(4 * N) + 3>(arg)), ...);
}



/// provide a reference to sum of N 4-momenta
/// with arg checking to minimize the reassignment i.e. reuse the p4 as much as possible
template <size_t ...N, typename ...Numbers>
const TLorentzVector& multivector_system_impl(std::index_sequence<N...>, Numbers ...numbers)
{
  static_assert(sizeof...(N) > 0 and sizeof...(numbers) == sizeof...(N) and sizeof...(N) % 4 == 0, 
                "ERROR: multivector_system: arguments are interpreted as pt, eta, phi and mass of "
                "the vectors comprising the system. As such, the number of arguments must be 4N!! "
                "Provide arguments as pt1, eta1, phi1, mass1, pt2, eta2, phi2, mass2...");

  using Number = typename std::tuple_element<0, std::tuple<Numbers...>>::type;

  static TLorentzVector psum;
  static std::array<Number, sizeof...(N)> arg;
  if (((std::get<N>(arg) == numbers) and ...))
    return psum;

  ((std::get<N>(arg) = numbers), ...);

  static std::array<TLorentzVector, sizeof...(N) / 4> p4s;
  assign_vector(p4s, arg, std::make_index_sequence<sizeof...(N) / 4>{});
  psum = p4s[0];
  for (int ip4 = 1; ip4 < sizeof...(N) / 4; ++ip4)
    psum += p4s[ip4];

  return psum;
}



/// whenever such a p4 sum is desired, without having to deal with details
/// for use mainly in the Aggregate attributes that require index masking
template <typename ...Numbers>
const TLorentzVector& multivector_system(Numbers ...numbers)
{
  return multivector_system_impl(std::make_index_sequence<sizeof...(numbers)>{}, numbers...);
}



/// implementation of function returning invariant mass of a p4
template <typename Number, typename ...Numbers>
Number invariant_mass_impl(Number number, Numbers ...numbers) 
{ 
  const auto &p4 = multivector_system(number, numbers...);
  return p4.M();
}



/// helper that returns a function pointer to invariant_mass_impl
template <typename Number, size_t ...N>
auto invariant_mass_helper(std::index_sequence<N...>) -> Number(*)(typename std::tuple_element<N, std::array<Number, sizeof...(N)>>::type...)
{
  return invariant_mass_impl<typename std::tuple_element<N, std::array<Number, sizeof...(N)>>::type...>;
}



/// a function returning a function that calculates the invariant mass of system consisting of N 4-momenta 
template <size_t N = 2, typename Number = float>
auto invariant_mass() -> decltype(invariant_mass_helper<Number>(std::make_index_sequence<4 * N>{}))
{
  return invariant_mass_helper<Number>(std::make_index_sequence<4 * N>{});
}



/// works similarly as invariant_mass_* above, but returns the transverse momentum instead
template <typename Number, typename ...Numbers>
Number system_pt_impl(Number number, Numbers ...numbers) 
{ 
  const auto &p4 = multivector_system(number, numbers...);
  return p4.Pt();
}



template <typename Number, size_t ...N>
auto system_pt_helper(std::index_sequence<N...>) -> Number(*)(typename std::tuple_element<N, std::array<Number, sizeof...(N)>>::type...)
{
  return system_pt_impl<typename std::tuple_element<N, std::array<Number, sizeof...(N)>>::type...>;
}



template <size_t N = 2, typename Number = float>
auto system_pt() -> decltype(system_pt_helper<Number>(std::make_index_sequence<4 * N>{}))
{
  return system_pt_helper<Number>(std::make_index_sequence<4 * N>{});
}



/// works similarly as invariant_mass_* above, but returns the rapidity instead
template <typename Number, typename ...Numbers>
Number system_rapidity_impl(Number number, Numbers ...numbers) 
{ 
  const auto &p4 = multivector_system(number, numbers...);
  return p4.Rapidity();
}



template <typename Number, size_t ...N>
auto system_rapidity_helper(std::index_sequence<N...>) -> Number(*)(typename std::tuple_element<N, std::array<Number, sizeof...(N)>>::type...)
{
  return system_rapidity_impl<typename std::tuple_element<N, std::array<Number, sizeof...(N)>>::type...>;
}



template <size_t N = 2, typename Number = float>
auto system_rapidity() -> decltype(system_rapidity_helper<Number>(std::make_index_sequence<4 * N>{}))
{
  return system_rapidity_helper<Number>(std::make_index_sequence<4 * N>{});
}

#endif
