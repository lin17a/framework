#ifndef FWK_HEAP_H
#define FWK_HEAP_H

// -*- C++ -*-
// author: afiq anuar
// short: collection of tools that are needed by group and friends

#include <string>

#include <memory>
#include <utility>

#include <vector>
#include <array>
#include <tuple>
#include <variant>

#include <algorithm>
#include <numeric>

#include <type_traits>
#include <functional>
#include <stdexcept>

// credit https://stackoverflow.com/questions/34099597/check-if-a-type-is-passed-in-variadic-template-parameter-pack
template<typename T, typename ...Ts>
constexpr bool contained_in = std::disjunction_v<std::is_same<T, Ts>...>;



// credit https://stackoverflow.com/questions/18986560/check-variadic-templates-parameters-for-uniqueness
template <typename T> 
struct Base {};

template <typename ...Ts>
struct Types : Base<Ts>...
{
  template <typename T>
  constexpr auto operator+(Base<T>)
  {
    if constexpr (std::is_base_of_v<Base<T>, Types>)
      return Types{};
    else
      return Types<Ts..., T>{};
  }

  constexpr size_t size() const
  {
    return sizeof...(Ts);
  }
};

template <typename ...Ts>
constexpr bool unique_types = ( (Types<>{} + ... + Base<Ts>{}).size() == sizeof...(Ts) );



// credit https://stackoverflow.com/questions/42580997/check-if-one-set-of-types-is-a-subset-of-the-other/
template <typename T, typename U>
constexpr bool is_subset_of = false;

template <template <typename, typename...> typename T, typename ...Ts, template <typename, typename...> typename U, typename ...Us>
constexpr bool is_subset_of<T<Ts...>, U<Us...>> 
= (contained_in<Ts, Us...> and ...);

template <typename T, typename U>
constexpr bool mutual_overlap = false;

template <template <typename, typename...> typename T, typename ...Ts, template <typename, typename...> typename U, typename ...Us>
constexpr bool mutual_overlap<T<Ts...>, U<Us...>> 
= (is_subset_of<T<Ts...>, U<Us...>> or is_subset_of<U<Us...>, T<Ts...>>);



// for std::visit
// credit https://en.cppreference.com/w/cpp/utility/variant/visit
template<typename ...Ts> struct overload : Ts... { using Ts::operator()...; };
template<typename ...Ts> overload(Ts...) -> overload<Ts...>;
template <typename> inline constexpr bool always_false_v = false;


// credit https://stackoverflow.com/questions/7943525/is-it-possible-to-figure-out-the-parameter-type-and-return-type-of-a-lambda
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{};

template <typename Ret, typename ...Args>
struct function_traits<Ret(*)(Args...)> : public function_traits<Ret(Args...)>
{};

// specialization for normal lambda
template <typename Cls, typename Ret, typename ...Args>
struct function_traits<Ret(Cls::*)(Args...) const>
{
  static constexpr int arity = sizeof...(Args);
  using result_type = Ret;
  using tuple_arg_types = std::tuple<Args...>;
  using tuple_arg_bare_types = std::tuple<std::remove_cv_t<std::remove_reference_t<Args>>...>;

  template <size_t i>
  using arg = typename std::tuple_element<i, tuple_arg_types>::type;

  template <size_t i>
  using bare_arg = typename std::tuple_element<i, tuple_arg_bare_types>::type;
};

// specialization for mutable lambda
template <typename Cls, typename Ret, typename ...Args>
struct function_traits<Ret(Cls::*)(Args...)>
{
  static constexpr int arity = sizeof...(Args);
  using result_type = Ret;
  using tuple_arg_types = std::tuple<Args...>;
  using tuple_arg_bare_types = std::tuple<std::remove_cv_t<std::remove_reference_t<Args>>...>;

  template <size_t i>
  using arg = typename std::tuple_element<i, tuple_arg_types>::type;

  template <size_t i>
  using bare_arg = typename std::tuple_element<i, tuple_arg_bare_types>::type;
};

// specialization for regular functions
template <typename Ret, typename ...Args>
struct function_traits<Ret(Args...)>
{
  static constexpr int arity = sizeof...(Args);
  using result_type = Ret;
  using tuple_arg_types = std::tuple<Args...>;
  using tuple_arg_bare_types = std::tuple<std::remove_cv_t<std::remove_reference_t<Args>>...>;

  template <size_t i>
  using arg = typename std::tuple_element<i, tuple_arg_types>::type;

  template <size_t i>
  using bare_arg = typename std::tuple_element<i, tuple_arg_bare_types>::type;
};



// make a tuple to references to Group data per function arg types
template <typename Tuple, typename Traits, std::size_t ...Is>
auto tuple_of_ref(const Tuple &tuple, Traits, std::index_sequence<Is...>)
{
  return std::make_tuple( std::ref(std::get<std::vector<typename Traits::template bare_arg<Is>>>( std::get<Is>(tuple) ))... );
}



// credit https://stackoverflow.com/questions/11322095/how-to-make-a-function-that-zips-two-tuples-in-c11-stl
// given tuples vec = {vec1, vec2, ... vecN} and idx = {i1, i2, ... iN}, where vec is anything that has the subscript operator []
// zip_index returns another tuple of references tup = {vec1[i1], vec2[i2], ... vecN[iN]}
// which can be passed to functions with std::apply
// _nn_ in name refers to "mapping n to n"
template <template <typename ...> typename T1, typename ...T1s, template <typename, size_t> typename A2, typename T2, std::size_t ...N>
constexpr auto zip_nn_helper(const T1<T1s...> &t1, const A2<T2, sizeof...(N)> &a2, std::index_sequence<N...>) 
-> decltype(std::forward_as_tuple( std::get<N>(t1)[ std::get<N>(a2) ]... ))
{
  return std::forward_as_tuple( std::get<N>(t1)[ std::get<N>(a2) ]... );
}

template <template <typename ...> typename T1, typename ...T1s, template <typename, size_t> typename A2, typename T2, std::size_t N>
constexpr auto zip_nn(const T1<T1s...> &t1, const A2<T2, N> &a2) 
-> decltype( zip_nn_helper(t1, a2, std::make_index_sequence<N>{}) ) 
{
  static_assert(sizeof...(T1s) == N, "ERROR: the tuple sizes must be the compatible!");
  return zip_nn_helper(t1, a2, std::make_index_sequence<N>{});
}

// a version implementing tup = {vec[i1], vec[i2], ... vec[N]}
// FIXME funny that _nn_ works in Aggregate:add_attribute swapping the forward_as_tuple to make_tuple
// FIXME but the same is not true for _1n_ in Group::transform_attribute
// FIXME make_tuple( std::ref(...)... ) works, but forward_as_tuple is chosen as it's more homogenous
template <typename T1, template <typename, size_t> typename A2, typename T2, std::size_t ...N>
constexpr auto zip_1n_helper(T1 &t1, const A2<T2, sizeof...(N)> &a2, std::index_sequence<N...>) 
-> decltype(std::forward_as_tuple( t1[ std::get<N>(a2) ]... ))
{
  return std::forward_as_tuple( t1[ std::get<N>(a2) ]... );
}

template <typename T1, template <typename, size_t> typename A2, typename T2, std::size_t N>
constexpr auto zip_1n(T1 &t1, const A2<T2, N> &a2) 
-> decltype( zip_1n_helper(t1, a2, std::make_index_sequence<N>{}) ) 
{
  return zip_1n_helper(t1, a2, std::make_index_sequence<N>{});
}

#endif
