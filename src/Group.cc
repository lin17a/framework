// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename ...Ts>
Framework::Group<Ts...>::Group(const std::string &name_, int counter_) :
name(name_),
counter(counter_),
selected(counter_),
v_index(this)
{
  if (counter > 0) {
    for (int iC = 0; iC < counter; ++iC)
      v_index.emplace_back(iC);
  }
}



template <typename ...Ts>
int Framework::Group<Ts...>::n_elements() const
{
  return selected;
}



template <typename ...Ts>
const int& Framework::Group<Ts...>::ref_to_n_elements() const
{
  return selected;
}



template <typename ...Ts>
int& Framework::Group<Ts...>::mref_to_n_elements()
{
  return selected;
}



template <typename ...Ts>
int Framework::Group<Ts...>::n_attributes() const
{
  return v_attr.size();
}



template <typename ...Ts>
bool Framework::Group<Ts...>::has_attribute(const std::string &name) const
{
  return inquire(name) != -1;
}



template <typename ...Ts>
void Framework::Group<Ts...>::reserve(int attr)
{
  v_attr.reserve(attr);
  v_data.reserve(attr);
}



template <typename ...Ts>
template <typename Function, typename ...Attributes>
bool Framework::Group<Ts...>::transform_attribute(const std::string &attr, Function function, Attributes &&...attrs)
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::transform_attribute requires some attributes to be provided!!");

  using Traits = function_traits<decltype(function)>;
  static_assert(contained_in<typename Traits::result_type, Ts...>, 
                "ERROR: Group::transform_attribute: the function return type is not among the types expected by the Group!!");
  static_assert(mutual_overlap<typename Traits::tuple_arg_bare_types, Group<Ts...>>, 
                "ERROR: Group::transform_attribute: the function argument types do not match the types expected by the Group!!");

  if (has_attribute(attr))
    return false;

  auto iA = (has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Group::transform_attribute: some of the requested attributes are not within the group!!" );

  // note on the structure
  // the three stage function execution is needed because
  // f_loop runs on actual data and produces the result
  // f_apply matches the attribute indices and forward the relevant refs to f_loop
  // all the functions need to be copied instead of referred 
  // due to scoping and/or lambda vs function pointer support

  auto f_loop = [function, this] (auto &vec, const auto &...vecs) -> void {
    for (int iE = 0; iE < this->counter; ++iE)
      vec[iE] = function(vecs[iE]...);
  };

  const std::array<int, sizeof...(attrs)> iattrs = {inquire(attrs)...};
  retype_per_function(zip_1n(v_data, iattrs), Traits{}, std::make_index_sequence<Traits::arity>{});

  auto f_apply = [f_loop, this, iattr = v_data.size(), iattrs] () -> void {
    auto refs = std::tuple_cat(std::make_tuple(std::ref( std::get<std::vector<typename Traits::result_type>>(v_data[iattr]) )), 
                               tuple_of_ref( zip_1n(v_data, iattrs), Traits{}, std::make_index_sequence<Traits::arity>{}) );

    std::apply(f_loop, refs);
  };

  v_attr.emplace_back(std::make_pair(attr, std::function<void()>(f_apply)));
  v_data.emplace_back(std::vector<typename Traits::result_type>());

  return true;
}



template <typename ...Ts>
std::vector<std::string> Framework::Group<Ts...>::attributes() const
{
  std::vector<std::string> v_attr_name;
  for (const auto &attr : v_attr)
    v_attr_name.emplace_back(attr.first);
  return v_attr_name;
}



template <typename ...Ts>
const std::vector<std::variant<std::vector<Ts>...>>& Framework::Group<Ts...>::data() const
{
  return v_data;
}



template <typename ...Ts>
const std::variant<std::vector<Ts>...>& Framework::Group<Ts...>::operator()(const std::string &name) const
{
  auto iA = inquire(name);
  if (iA == -1)
    throw std::invalid_argument( "ERROR: Group::get: requested attribute " + name + " is not within the group!!" );

  return v_data[iA];
}



template <typename ...Ts>
const std::variant<std::vector<Ts>...>& Framework::Group<Ts...>::operator()(int iattr) const
{
  return v_data[iattr];
}



template <typename ...Ts>
std::variant<std::vector<Ts>...>& Framework::Group<Ts...>::mref_to_attribute(const std::string &name)
{
  return const_cast<std::variant<std::vector<Ts>...>&>( (*const_cast<const Framework::Group<Ts...>*>(this))(name) );
}



template <typename ...Ts>
template <typename T>
const std::vector<T>& Framework::Group<Ts...>::get(const std::string &name) const
{
  static_assert(contained_in<T, Ts...>, "ERROR: Group::get: called with a type not among the types of by the Group!!");
  return std::get<std::vector<T>>((*this)(name));
}



template <typename ...Ts>
template <typename T>
const T& Framework::Group<Ts...>::get(const std::string &name, int index) const
{
  if (index >= selected)
    throw std::invalid_argument( "ERROR: Group::get: an invalid element index is requested!!" );

  return this->get<T>(name)[v_index[index]];
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::indices() const
{
  return v_index;
}



template <typename ...Ts>
const typename Framework::Group<Ts...>::idxs& Framework::Group<Ts...>::ref_to_indices() const
{
  return v_index;
}



template <typename ...Ts>
int Framework::Group<Ts...>::index(int order) const
{
  int index = (order > -1 and order < selected) ? v_index[order] : -1;
  return index;
}



template <typename ...Ts>
void Framework::Group<Ts...>::update_indices(const typename Framework::Group<Ts...>::idxs &v_idx)
{
  v_index.clear();
  for (auto idx : v_idx)
    v_index.emplace_back(idx);
  selected = v_index.size();
}



template <typename ...Ts>
template <typename Function, typename ...Attributes>
void Framework::Group<Ts...>::iterate(Function function, const typename Framework::Group<Ts...>::idxs &v_idx, Attributes &&...attrs) const
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::iterate makes no sense without specifying attributes!!");

  auto iA = (has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Group::iterate: some of the requested attributes are not within the group!!" );

  if (v_idx.ref != this)
    throw std::invalid_argument( "ERROR: Group::iterate: requesting iteration over index set of another Group is nonsensical!!" );

  std::visit([&function, &v_idx, this] (const auto &...vec) {
      for (auto idx : v_idx)
        function(vec[idx]...);
    }, v_data[inquire(attrs)]...);
}



template <typename ...Ts>
template <typename Compare, typename ...Attributes>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                       Compare compare, Attributes &&...attrs) const
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::filter makes no sense without specifying attributes!!");

  auto iA = (has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Group::filter some of the requested attributes are not within the group!!" );

  if (v_idx.ref == this and v_idx.empty())
    return v_idx;

  return filter_helper( (v_idx.ref != this) ? v_index : v_idx, compare, inquire(attrs)... );
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_less(const std::string &name, Number value, 
                                                                            const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data < value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_less_equal(const std::string &name, Number value, 
                                                                                  const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data <= value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_greater(const std::string &name, Number value, 
                                                                               const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data > value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_greater_equal(const std::string &name, Number value, 
                                                                                     const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data >= value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_equal(const std::string &name, Number value, 
                                                                             const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data == value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_not(const std::string &name, Number value, 
                                                                           const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data != value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_absolute_equal(const std::string &name, Number value, 
                                                                                      const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return std::abs(data) == value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_absolute_not(const std::string &name, Number value, 
                                                                                    const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return std::abs(data) != value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_bit_and(const std::string &name, Number value, 
                                                                               const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {
      if constexpr(std::is_integral_v<std::remove_cv_t<std::remove_reference_t<decltype(data)>>> and 
                   std::is_integral_v<std::remove_cv_t<std::remove_reference_t<Number>>>)
                    return (data & value);
      else
        return false;
    }, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_in(const std::string &name, Number min, Number max, 
                                                                          const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&min, &max] (auto &data) {return (data > min and data < max);}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_out(const std::string &name, Number min, Number max, 
                                                                           const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&min, &max] (auto &data) {return (data < min and data > max);}, name);
}


template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_3values(const std::string &name, Number value1, Number value2, Number value3) const
{
  return filter([&value1, &value2, &value3] (auto &data) {return ((data == value1 or data == value2) or data == value3);}, name);
}


template <typename ...Ts>
template <typename Compare>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                     Compare compare, const std::string &name) const
{
  auto iA = inquire(name);
  if (iA == -1)
    throw std::invalid_argument( "ERROR: Group::sort: some of the requested attributes are not within the group!!" );

  if (v_idx.ref == this and v_idx.size() < 2)
    return v_idx;

  return sort_helper( (v_idx.ref != this) ? v_index : v_idx, compare, iA );
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_ascending(const std::string &name, 
                                                                               const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (p1.second < p2.second); }, name);
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_descending(const std::string &name, 
                                                                                const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (p1.second > p2.second); }, name);
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_absolute_ascending(const std::string &name, 
                                                                                        const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (std::abs(p1.second) < std::abs(p2.second)); }, name);
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_absolute_descending(const std::string &name, 
                                                                                         const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (std::abs(p1.second) > std::abs(p2.second)); }, name);
}



template <typename ...Ts>
int Framework::Group<Ts...>::inquire(const std::string &name) const
{
  for (int iA = 0; iA < v_attr.size(); ++iA) {
    auto &[alias, _] = v_attr[iA];
    (void) _;

    if (alias == name)
      return iA;
  }

  return -1;
}



template <typename ...Ts>
void Framework::Group<Ts...>::reorder()
{
  for (int iS = 0; iS < selected; ++iS) {
    if (iS != v_index[iS]) {
      for (auto &dat : v_data)
        std::visit([iS, iI = v_index[iS]] (auto &vec) {std::swap(vec[iS], vec[iI]);}, dat);

      v_index[iS] = iS;
    }
  }
}



template <typename ...Ts>
void Framework::Group<Ts...>::initialize(int init)
{
  v_index.reserve(init);

  for (auto &dat : v_data)
    std::visit([init] (auto &vec) {vec.reserve(init); vec.clear();}, dat);
}



template <typename ...Ts>
template <typename Number>
void Framework::Group<Ts...>::retype(std::variant<std::vector<Ts>...> &dat)
{
  if constexpr (contained_in<Number, Ts...>) {
      if (std::get_if<std::vector<Number>>(&dat) == nullptr) {
        dat = std::vector<Number>();
        std::visit([init = v_index.capacity()] (auto &vec) {vec.reserve(init); vec.clear();}, dat);
      }
    }
}



template <typename ...Ts>
template <typename Tuple, typename Traits, std::size_t ...Is>
void Framework::Group<Ts...>::retype_per_function(const Tuple &tuple, Traits, std::index_sequence<Is...>)
{
  (retype<typename Traits::template bare_arg<Is>>(std::get<Is>(tuple)), ...);
}



template <typename ...Ts>
template <typename Compare, typename ...Attributes>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_helper(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                              Compare &compare, Attributes &&...attrs) const
{
  typename Framework::Group<Ts...>::idxs v_pass(this);
  std::visit([&v_pass, &v_idx, &compare] (const auto &...vec) {
      for (auto index : v_idx) {
        if (compare(vec[index]...))
          v_pass.emplace_back(index);
      }
    }, v_data[attrs]...);

  return v_pass;
}



template <typename ...Ts>
template <typename Compare>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_helper(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                            Compare &compare, int attr) const
{
  typename Framework::Group<Ts...>::idxs v_sort(this);
  std::visit([&v_sort, &v_idx, &compare] (const auto &vec) {
      using VT = typename std::decay_t<decltype(vec)>::value_type;
      std::vector<std::pair<int, VT>> v_zip;
      for (auto index : v_idx)
        v_zip.emplace_back(index, vec[index]);

      std::sort(std::begin(v_zip), std::end(v_zip), compare);

      for (auto &zip : v_zip)
        v_sort.emplace_back(zip.first);
    }, v_data[attr]);

  return v_sort;
}
