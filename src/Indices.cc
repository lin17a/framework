// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename Group>
Framework::Indices<Group>::Indices(const Group *group, const std::vector<int> &idxs) :
ref(group),
v_index(idxs)
{ /*static int idx = 0; std::cout << "indices construction " << ++idx << "\n";*/ }



template <typename Group>
int& Framework::Indices<Group>::operator[](int idx)
{
  return v_index[idx];
}



template <typename Group>
const int& Framework::Indices<Group>::operator[](int idx) const
{
  return v_index[idx];
}



template <typename Group>
std::vector<int>::iterator Framework::Indices<Group>::begin() noexcept
{
  return v_index.begin();
}



template <typename Group>
std::vector<int>::const_iterator Framework::Indices<Group>::begin() const noexcept
{
  return v_index.begin();
}


template <typename Group>
std::vector<int>::const_iterator Framework::Indices<Group>::cbegin() const noexcept
{
  return v_index.cbegin();
}



template <typename Group>
std::vector<int>::iterator Framework::Indices<Group>::end() noexcept
{
  return v_index.end();
}



template <typename Group>
std::vector<int>::const_iterator Framework::Indices<Group>::end() const noexcept
{
  return v_index.end();
}



template <typename Group>
std::vector<int>::const_iterator Framework::Indices<Group>::cend() const noexcept
{
  return v_index.cend();
}



template <typename Group>
int Framework::Indices<Group>::size() const noexcept
{
  return v_index.size();
}



template <typename Group>
bool Framework::Indices<Group>::empty() const noexcept
{
  return this->size() == 0;
}



template <typename Group>
void Framework::Indices<Group>::reserve(int capacity)
{
  v_index.reserve(capacity);
}



template <typename Group>
int Framework::Indices<Group>::capacity() const noexcept
{
  return v_index.capacity();
}



template <typename Group>
void Framework::Indices<Group>::clear() noexcept
{
  v_index.clear();
}



template <typename Group>
void Framework::Indices<Group>::emplace_back(int idx)
{
  v_index.emplace_back(idx);
}



template <typename Group>
bool Framework::Indices<Group>::operator==(const Indices &other) const
{
  return (ref == other.ref) ? v_index == other.v_index : false;
}



template <typename Group>
template <typename Compare, typename ...Attributes>
Framework::Indices<Group> Framework::Indices<Group>::filter(Compare &compare, Attributes &&...attrs) const
{
  return ref->filter(*this, compare, attrs...);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_less(const std::string &name, Number value) const
{
  return ref->filter_less(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_less_equal(const std::string &name, Number value) const
{
  return ref->filter_less_equal(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_greater(const std::string &name, Number value) const
{
  return ref->filter_greater(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_greater_equal(const std::string &name, Number value) const
{
  return ref->filter_greater_equal(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_equal(const std::string &name, Number value) const
{
  return ref->filter_equal(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_not(const std::string &name, Number value) const
{
  return ref->filter_not(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_absolute_equal(const std::string &name, Number value) const
{
  return ref->filter_absolute_equal(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_absolute_not(const std::string &name, Number value) const
{
  return ref->filter_absolute_not(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_bit_and(const std::string &name, Number value) const
{
  return ref->filter_bit_and(name, value, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_in(const std::string &name, Number min, Number max) const
{
  return ref->filter_in(name, min, max, *this);
}



template <typename Group>
template <typename Number>
Framework::Indices<Group> Framework::Indices<Group>::filter_out(const std::string &name, Number min, Number max) const
{
  return ref->filter_out(name, min, max, *this);
}



template <typename Group>
template <typename Compare>
Framework::Indices<Group> Framework::Indices<Group>::sort(Compare compare, const std::string &name) const
{
  return ref->sort(*this, compare, name);
}



template <typename Group>
Framework::Indices<Group> Framework::Indices<Group>::sort_ascending(const std::string &name) const
{
  return ref->sort_ascending(name, *this);
}



template <typename Group>
Framework::Indices<Group> Framework::Indices<Group>::sort_descending(const std::string &name) const
{
  return ref->sort_descending(name, *this);
}



template <typename Group>
Framework::Indices<Group> Framework::Indices<Group>::sort_absolute_ascending(const std::string &name) const
{
  return ref->sort_absolute_ascending(name, *this);
}



template <typename Group>
Framework::Indices<Group> Framework::Indices<Group>::sort_absolute_descending(const std::string &name) const
{
  return ref->sort_absolute_descending(name, *this);
}



template <typename Idx, typename ...Idxs>
Idx Framework::merge(const Idx &index, const Idxs &...indices)
{
  static_assert(sizeof...(indices) > 0, 
                "ERROR: merge takes at least 2 index sets!!");

  static_assert(std::conjunction_v<std::is_same<Idx, Idxs>...>, 
                "ERROR: merge: requesting the merging of incompatible index set types!!");

  if (((index.ref != indices.ref) or ...))
    return Idx();

  Idx result(index);
  auto check_and_put = [&result, &index] (const auto &idxs) {
    for (const auto &idx : idxs) {
      if (!std::count(std::begin(index), std::end(index), idx))
        result.emplace_back(idx);
    }
  };
  (( check_and_put(indices) ), ...);

  return result;
}
