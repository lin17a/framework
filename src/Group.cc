// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename ...Ts>
Framework::Group<Ts...>::Group(const std::string &name_, int counter_) :
name(name_),
counter(counter_),
selected(counter_)
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
std::vector<int> Framework::Group<Ts...>::indices() const
{
  return v_index;
}



template <typename ...Ts>
const std::vector<int>& Framework::Group<Ts...>::ref_to_indices() const
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
void Framework::Group<Ts...>::update_indices(const std::vector<int> &v_idx)
{
  v_index.clear();
  for (auto idx : v_idx)
    v_index.emplace_back(idx);
  selected = v_index.size();
}



template <typename ...Ts>
template <typename Function, typename ...Attributes>
void Framework::Group<Ts...>::iterate(Function function, int begin, int end, Attributes &&...attrs) const
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::iterate makes no sense without specifying attributes!!");

  begin = (begin < 0 or begin > selected) ? 0 : begin;
  end = (end < begin + 1 or end > selected) ? selected : end;

  auto iA = (has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Group::iterate: some of the requested attributes are not within the group!!" );

  std::visit([&function, &begin, &end, this] (const auto &...vec) {
      for (int iE = begin; iE < end; ++iE)
        function(vec[this->v_index[iE]]...);
    }, v_data[inquire(attrs)]...);
}



template <typename ...Ts>
template <typename Compare, typename ...Attributes>
std::vector<int> Framework::Group<Ts...>::filter(Compare compare, Attributes &&...attrs) const
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::filter makes no sense without specifying attributes!!");

  auto iA = (has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Group::filter some of the requested attributes are not within the group!!" );

  return filter_helper( compare, inquire(attrs)... );
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_less(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {return data < value;}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_less_equal(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {return data <= value;}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_greater(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {return data > value;}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_greater_equal(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {return data >= value;}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_equal(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {return data == value;}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_not(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {return data != value;}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_bit_and(const std::string &name, Number value) const
{
  return filter([&value] (auto &data) {
      if constexpr(std::is_integral_v<std::remove_cv_t<std::remove_reference_t<decltype(data)>>> and 
                   std::is_integral_v<std::remove_cv_t<std::remove_reference_t<Number>>>)
                    return (data & value);
      else
        return false;
    }, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_in(const std::string &name, Number min, Number max) const
{
  return filter([&min, &max] (auto &data) {return (data > min and data < max);}, name);
}



template <typename ...Ts>
template <typename Number>
std::vector<int> Framework::Group<Ts...>::filter_out(const std::string &name, Number min, Number max) const
{
  return filter([&min, &max] (auto &data) {return (data < min and data > max);}, name);
}



template <typename ...Ts>
template <typename Compare, typename ...Attributes>
int Framework::Group<Ts...>::count(Compare compare, Attributes &&...attrs) const
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::count makes no sense without specifying attributes!!");

  auto iA = (has_attribute(attrs) and ...);
  if (!iA)
    throw std::invalid_argument( "ERROR: Group::count some of the requested attributes are not within the group!!" );

  return filter_helper( compare, inquire(attrs)... ).size();
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_less(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {return data < value;}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_less_equal(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {return data <= value;}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_greater(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {return data > value;}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_greater_equal(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {return data >= value;}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_equal(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {return data == value;}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_not(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {return data != value;}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_bit_and(const std::string &name, Number value) const
{
  return count([&value] (auto &data) {
      if constexpr(std::is_integral_v<std::remove_cv_t<std::remove_reference_t<decltype(data)>>> and 
                   std::is_integral_v<std::remove_cv_t<std::remove_reference_t<Number>>>)
                    return (data & value);
      else
        return false;
    }, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_in(const std::string &name, Number min, Number max) const
{
  return count([&min, &max] (auto &data) {return (data > min and data < max);}, name);
}



template <typename ...Ts>
template <typename Number>
int Framework::Group<Ts...>::count_out(const std::string &name, Number min, Number max) const
{
  return count([&min, &max] (auto &data) {return (data < min and data > max);}, name);
}



template <typename ...Ts>
template <typename Compare>
std::vector<int> Framework::Group<Ts...>::sort(Compare compare, const std::string &name) const
{
  auto iA = inquire(name);
  if (iA == -1)
    throw std::invalid_argument( "ERROR: Group::sort: some of the requested attributes are not within the group!!" );

  return sort_helper( compare, iA );
}



template <typename ...Ts>
std::vector<int> Framework::Group<Ts...>::sort_ascending(const std::string &name) const
{
  return sort([] (const auto &p1, const auto &p2) { return (p1.second < p2.second); }, name);
}



template <typename ...Ts>
std::vector<int> Framework::Group<Ts...>::sort_descending(const std::string &name) const
{
  return sort([] (const auto &p1, const auto &p2) { return (p1.second > p2.second); }, name);
}



template <typename ...Ts>
std::vector<int> Framework::Group<Ts...>::sort_absolute_ascending(const std::string &name) const
{
  return sort([] (const auto &p1, const auto &p2) { return (std::abs(p1.second) < std::abs(p2.second)); }, name);
}



template <typename ...Ts>
std::vector<int> Framework::Group<Ts...>::sort_absolute_descending(const std::string &name) const
{
  return sort([] (const auto &p1, const auto &p2) { return (std::abs(p1.second) > std::abs(p2.second)); }, name);
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
template <typename Compare, typename ...Attributes>
std::vector<int> Framework::Group<Ts...>::filter_helper(Compare &compare, Attributes &&...attrs) const
{
  std::vector<int> v_idx;
  std::visit([this, &v_idx, &compare] (const auto &...vec) {
      for (auto &index : this->v_index) {
        if (compare(vec[index]...))
          v_idx.emplace_back(index);
      }
    }, v_data[attrs]...);

  return v_idx;
}



template <typename ...Ts>
template <typename Compare>
std::vector<int> Framework::Group<Ts...>::sort_helper(Compare &compare, int attr) const
{
  std::vector<int> v_idx;
  std::visit([this, &v_idx, &compare] (const auto &vec) {
      using VT = typename std::decay_t<decltype(vec)>::value_type;
      std::vector<std::pair<int, VT>> v_zip;
      for (auto &index : this->v_index)
        v_zip.emplace_back(index, vec[index]);

      std::sort(std::begin(v_zip), std::end(v_zip), compare);

      for (auto &zip : v_zip)
        v_idx.emplace_back(zip.first);
    }, v_data[attr]);

  return v_idx;
}
