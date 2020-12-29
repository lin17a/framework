#ifndef FWK_GROUP_H
#define FWK_GROUP_H

// -*- C++ -*-
// author: afiq anuar
// short: group of objects i.e. sets of attributes e.g. a particle has pt eta phi...
// note: attributes that are self-referencing e.g. GenPart_motherIdx can't be handled by iterate(); for this one needs to use operator()

#include "Heap.h"
#include "Indices.h"

// https://stackoverflow.com/questions/670308/alternative-to-vectorbool
class boolean {
 public:
   boolean(): value() {}
   boolean(bool value_) : value(value_) {}

   operator bool() const {return value;}

   /// the following operators are to allow bool* b = &v[0]; (v is a vector here)
   bool* operator& () { return &value; }
   const bool* operator& () const { return &value; }

 private:
   bool value;
};

namespace Framework {
  template <typename ...Ts>
  class Group {
    static_assert(unique_types<Ts...>, "ERROR: a Group must be initialized with unique types!");
    static_assert(!contained_in<bool, Ts...>, "ERROR: Group relies on a few features that are incompatible with standard C++ "
                  "when bool is among the included types. For boolean attributes please use the 'boolean' type instead, which is a drop-in "
                  "replacement provided precisely to avoid this quirk of the standard.");

  public:
    using base = Group<Ts...>;
    using idxs = Indices<base>;

    /// no default constructor
    Group() = delete;

    /// constructor
    Group(const std::string &name_, int counter_);

    /// number of currently held elements
    int n_elements() const;

    /// ref instead of copy of the above
    const int& ref_to_n_elements() const;

    /// a mutable ref version
    /// can't be const if it's to be used to write TTree...
    /// might be worth considering to write TTree using copies rather than in-place references?
    int& mref_to_n_elements();

    /// number of currently held attributes
    int n_attributes() const;

    /// as it says on the tin
    bool has_attribute(const std::string &name) const;

    /// reserve the space for expected number of attributes
    void reserve(int attr);

    /// transform a group of internal attributes into another attribute
    /// the transformation is done element-wise on every element of held data
    /// as this adds a new attribute, its name has to be unique
    template <typename Function, typename ...Attributes>
    bool transform_attribute(const std::string &attr, Function function, Attributes &&...attrs);

    /// list of attributes
    std::vector<std::string> attributes() const;

    /// reference to container of elements
    const std::vector<std::variant<std::vector<Ts>...>>& data() const;

    /// reference to single attribute array - variant version
    const std::variant<std::vector<Ts>...>& operator()(const std::string &name) const;

    /// mutable version of the above
    /// only one version provided, intended for use by Tree only
    std::variant<std::vector<Ts>...>& mref_to_attribute(const std::string &name);

    /// reference to single attribute array - typed version
    template <typename T>
    const std::vector<T>& get(const std::string &name) const;

    /// reference to single element in an attribute - typed version
    /// equivalent to get<T>(attr)[v_index[index]]
    /// i.e. gives the nth element as per the current Group state accounting for previous update_indices calls 
    template <typename T>
    const T& get(const std::string &name, int index) const;

    /// the associated indices to be used with the above
    idxs indices() const;

    /// ref instead of copy of the above
    const idxs& ref_to_indices() const;

    /// and individual index access i.e. v_index[order] with bounds checking
    int index(int order) const;

    /// update the indices with another set e.g. output of filter/sort
    /// no checking if the indices are actually legit
    void update_indices(const idxs &v_idx);

    /// element iterator taking a function and runs the visitor over it
    /// can be type-dependent or otherwise
    /// other args are the attribute and element count to iterate over
    /// -1 or anything > v_index.size() to iterate over everything
    template <typename Function, typename ...Attributes>
    void iterate(Function function, const idxs &v_idx, Attributes &&...attrs) const;

    /// filter the elements in the collection by some criteria on a given attribute
    /// custom filter needs a function returning a bool and taking two args of type Number
    /// the first arg is the indices to restrict the filtering over
    /// if it equals idxs(), no restriction is imposed i.e. the current Group indices are used
    /// second arg is the function running on (single elements of) the attributes
    /// following args are the attributes themselves
    /// returns the indices after filtering
    /// this method DOES NOT update the index in-place
    template <typename Compare, typename ...Attributes>
    idxs filter(const idxs &v_idx, Compare compare, Attributes &&...attrs) const;

    /// common filters
    template <typename Number>
    idxs filter_less(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_less_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_greater(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_greater_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_not(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_absolute_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_absolute_not(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_bit_and(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    /// both are min and max exclusive
    template <typename Number>
    idxs filter_in(const std::string &name, Number min, Number max, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_out(const std::string &name, Number min, Number max, const idxs &v_idx = idxs()) const;

    /// sort the elements in the collection by a given attribute
    /// custom sorter needs a function returning a bool and taking two args, both of std::pair<int, decltype(data)>
    /// FIXME prepare a more convenient implementation
    /// returns the sorted indices
    template <typename Compare>
    idxs sort(const idxs &v_idx, Compare compare, const std::string &name) const;

    /// common sorts
    idxs sort_ascending(const std::string &name, const idxs &v_idx = idxs()) const;

    idxs sort_descending(const std::string &name, const idxs &v_idx = idxs()) const;

    idxs sort_absolute_ascending(const std::string &name, const idxs &v_idx = idxs()) const;

    idxs sort_absolute_descending(const std::string &name, const idxs &v_idx = idxs()) const;

    /// returns the index where an attribute occurs
    int inquire(const std::string &name) const;

    /// populate the Group data
    virtual void populate(long long entry) = 0;

    /// reorder the group data such that selected elements occur in front
    /// selected elements are those whose index is in v_index
    void reorder();

    /// name of the group
    std::string name;

  protected:
    /// this method ensures that all attributes have the proper capacity
    void initialize(int init);

    /// helper that actually does the filtering
    template <typename Compare, typename ...Attributes>
    idxs filter_helper(const idxs &v_idx, Compare &compare, Attributes &&...attrs) const;

    /// helper that actually does the sorting
    template <typename Compare>
    idxs sort_helper(const idxs &v_idx, Compare &compare, int attr) const;

    /// element counter before prefiltering
    int counter;

    /// element counter after prefiltering i.e. v_index.size()
    int selected;

    /// element indices in the group
    idxs v_index;

    /// register of the currently available attributes
    /// first string is attribute alias
    /// second function is for the element-wise transformation from other attributes
    std::vector<std::pair<std::string, std::function<void()>>> v_attr;

    /// attribute storage
    std::vector<std::variant<std::vector<Ts>...>> v_data;
  };
}

#include "Group.cc"

#endif
