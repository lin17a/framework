// -*- C++ -*-
// author: afiq anuar
// short: a listing of free convenience functions pertaining to containers

#ifndef FWK_CONTAINER_UTIL_H
#define FWK_CONTAINER_UTIL_H

#include <vector>
#include <type_traits>

template <typename K, typename V, typename C, typename std::enable_if_t<std::is_convertible_v<C, K>>* = nullptr>
int index_with_key(const std::vector<std::pair<K, V>> &vec, const C &key)
{
  for (int iE = 0; iE < vec.size(); ++iE) {
    auto &alias = vec[iE].first;
    if (alias == key)
      return iE;
  }

  return -1;
}



/// fill an interval from the last element in vector to a number in some steps
/// step serves also as epsilon (should be ok for intended usage of this function)
template <typename Number, typename std::enable_if_t<std::is_floating_point_v<Number>>* = nullptr> 
void fill_interval(std::vector<Number> &vec, Number to, Number step)
{
  // do nothing in various cases the function won't make sense
  const Number eps = std::abs(step / 50.);
  if (step == 0. or vec.empty() or (step > 0. and to - vec.back() < step) or (step < 0. and to - vec.back() > step))
    return;

  while (std::abs(to - vec.back()) >= eps)
    vec.emplace_back(vec.back() + step);
}



/// also also for integral types where there's none of the epsilon shenanigans
template <typename Number, typename std::enable_if_t<std::is_integral_v<Number>>* = nullptr> 
void fill_interval(std::vector<Number> &vec, Number to, Number step)
{
  // do nothing in various cases the function won't make sense
  if (step == 0 or vec.empty() or (step > 0 and to - vec.back() < step) or (step < 0 and to - vec.back() > step))
    return;

  while (std::abs(to - vec.back()) >= std::abs(step))
    vec.emplace_back(vec.back() + step);
}



/// method that makes such an interval vector relying on the above
template <typename Number> 
std::vector<Number> make_interval(Number from, Number to, Number step)
{
  std::vector<Number> v_interval = {from};
  fill_interval(v_interval, to, step);
  return (v_interval.size() > 1) ? v_interval : {};
}

#endif
