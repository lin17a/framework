// -*- C++ -*-
// author: afiq anuar
// short: a listing of free convenience functions for use in the framework

#ifndef FWK_FUNCTION_UTIL_H
#define FWK_FUNCTION_UTIL_H

#include "../src/Heap.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// a function that performs a simple copy of the argument
template <typename T = float>
T identity(T t)
{
  return t;
}



// fillers for the histogram class
template <typename ...Groups>
auto filler_count(const Groups &...groups)
{
  static_assert(sizeof...(groups) > 0 and sizeof...(groups) < 4, "ERROR: filler_count: currently only 1D - 3D histograms are supported!!");

  using Hist = typename std::conditional_t<sizeof...(groups) != 1, typename std::conditional_t<sizeof...(groups) != 2, TH3, TH2>, TH1>;

  return [&groups...] (Hist *hist, const double &weight) {
    hist->Fill(groups.n_elements()..., weight);
  };
}



template <typename Group, typename ...Attributes>
auto filler_first_of(const Group &group, Attributes &&...attrs)
{
  static_assert(sizeof...(attrs) > 0 and sizeof...(attrs) < 4, "ERROR: filler_first_of: currently only 1D - 3D histograms are supported!!");

  using Hist = typename std::conditional_t<sizeof...(attrs) != 1, typename std::conditional_t<sizeof...(attrs) != 2, TH3, TH2>, TH1>;

  auto iattrs = std::make_tuple( group.inquire(attrs)... );
  auto filler = [&group] (auto &&...idxs) {
    return [&group, idxs...] (Hist *hist, const double &weight) {
      std::visit([&hist, &weight, &indices = group.ref_to_indices()] (const auto &...vec) {
          if (indices) 
            hist->Fill(vec[indices[0]]..., weight); 
        }, group(idxs)...);
    };
  };

  return std::apply(filler, iattrs);
}



template <typename Group, typename ...Attributes>
auto filler_all_of(const Group &group, Attributes &&...attrs)
{
  static_assert(sizeof...(attrs) > 0 and sizeof...(attrs) < 4, "ERROR: filler_all_of: currently only 1D - 3D histograms are supported!!");

  using Hist = typename std::conditional_t<sizeof...(attrs) != 1, typename std::conditional_t<sizeof...(attrs) != 2, TH3, TH2>, TH1>;

  auto iattrs = std::make_tuple( group.inquire(attrs)... );
  auto filler = [&group] (auto &&...idxs) {
    return [&group, idxs...] (Hist *hist, const double &weight) {
      std::visit([&hist, &weight, &indices = group.ref_to_indices()] (const auto &...vec) {
          for (int iE = 0; iE < indices.size(); ++iE)
            hist->Fill(vec[indices[iE]]..., weight);
        }, group(idxs)...);
    };
  };

  return std::apply(filler, iattrs);
}



// a or b or c or ... in function form
// call any_of<N> to get a function that takes N bools and return the OR of them all
template <typename ...Bools>
bool any_of_impl(Bools ...bools)
{
  return (bools or ...);
}



template <size_t ...N>
auto any_of_helper(std::index_sequence<N...>) -> bool(*)(typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...)
{
  return any_of_impl<typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...>;
}



template <size_t N = 1>
auto any_of() -> decltype(any_of_helper(std::make_index_sequence<N>{}))
{
  return any_of_helper(std::make_index_sequence<N>{});
}



// a and b and c and ... in function form
// call all_of<N> to get a function that takes N bools and return the AND of them all
template <typename ...Bools>
bool all_of_impl(Bools ...bools)
{
  return (bools and ...);
}



template <size_t ...N>
auto all_of_helper(std::index_sequence<N...>) -> bool(*)(typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...)
{
  return all_of_impl<typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...>;
}



template <size_t N = 1>
auto all_of() -> decltype(all_of_helper(std::make_index_sequence<N>{}))
{
  return all_of_helper(std::make_index_sequence<N>{});
}



// implementation of the apply_to<N, F>, refer to that for more info
template <typename Number, size_t ...I, typename Function, size_t ...N>
auto apply_to_impl(std::index_sequence<I...>, Function function, std::index_sequence<N...>)
{
  using Traits = function_traits<decltype(function)>;
  auto f_apply = [function] (typename std::tuple_element_t<N, std::array<Number, sizeof...(N)>> ...args) -> typename Traits::result_type { 
    return function( std::get<sizeof...(N) - Traits::arity + I>(std::array<Number, sizeof...(N)>{args...})... );
  };

  return f_apply;
}



// helper of apply_to<N, F>, refer to that for more info
template <typename Number, typename Function, size_t ...N>
auto apply_to_helper(Function function, std::index_sequence<N...>)
{
  using Traits = function_traits<decltype(function)>;
  return apply_to_impl<typename Traits::result_type>(std::make_index_sequence<Traits::arity>{}, function, std::make_index_sequence<sizeof...(N)>{});
}



// a function that outputs a function that applies another function on its last F arguments
// let ff be a function taking F arguments
// then a call of apply_to<N>(ff) returns a function that takes (N + 1) * F arguments
// whose result is the same as calling ff on the last F arguments, ignoring the first NF arguments
// there is no reason why would apply_to<0> not work; it's just prevented as no use case is envisioned for it
// being that this is written mainly to make index masking a touch easier
template <size_t N, typename Function>
auto apply_to(Function function)
{
  static_assert(N > 0, "ERROR: apply_to is not callable with N = 0, as in this case no masking is necessary!!");
  using Traits = function_traits<decltype(function)>;
  return apply_to_helper<typename Traits::result_type>(function, std::make_index_sequence<(N + 1) * Traits::arity>{});
}

#endif
