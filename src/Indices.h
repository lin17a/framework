#ifndef FWK_INDICES_H
#define FWK_INDICES_H

// -*- C++ -*-
// author: afiq anuar
// short: handling of indices; a dedicated struct is needed if we're to make filter calls cuter

namespace Framework {
  template <typename Group>
  class Indices {
    static_assert(std::is_class_v<typename Group::base>, "ERROR: Indices is meant to index a Group!!");

  public:
    /// constructor
    /// the default indexer is the one used by Group's as dummy
    Indices(const Group *group = nullptr, const std::vector<int> &idxs = {});

    /// many of these are just reproducing the neccessary vector interface that we use...
    /// as such only Indices-specific items will be elaborated upon
    /// being that one isn't really meant to edit the indices etc
    /// that it's not fully reproducing vector is deliberate
    int& operator[](int idx);

    const int& operator[](int idx) const;

    std::vector<int>::iterator begin() noexcept;

    std::vector<int>::const_iterator begin() const noexcept;

    std::vector<int>::const_iterator cbegin() const noexcept;

    std::vector<int>::iterator end() noexcept;

    std::vector<int>::const_iterator end() const noexcept;

    std::vector<int>::const_iterator cend() const noexcept;

    int size() const noexcept;

    bool empty() const noexcept;

    void reserve(int capacity);

    int capacity() const noexcept;

    void clear() noexcept;

    void emplace_back(int idx);

    bool operator==(const Indices &other) const;

    /// such that if (Indices) is a well-defined expression
    explicit operator bool() const { return this->ref != nullptr and !(this->empty()); };

    /// proxy methods to allow chaining of filter calls
    /// unlike Group::filter which allows arbitrary Indices,
    /// Indices::filter can only pass itself
    template <typename Compare, typename ...Attributes>
    Indices filter(Compare &compare, Attributes &&...attrs) const;

    template <typename Number>
    Indices filter_less(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_less_equal(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_greater(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_greater_equal(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_equal(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_not(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_absolute_equal(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_absolute_not(const std::string &name, Number value) const;

    template <typename Number>
    Indices filter_bit_and(const std::string &name, Number value) const;

    /// both are min and max exclusive
    template <typename Number>
    Indices filter_in(const std::string &name, Number min, Number max) const;

    template <typename Number>
    Indices filter_out(const std::string &name, Number min, Number max) const;

    template <typename Compare>
    Indices sort(Compare compare, const std::string &name) const;

    Indices sort_ascending(const std::string &name) const;

    Indices sort_descending(const std::string &name) const;

    Indices sort_absolute_ascending(const std::string &name) const;

    Indices sort_absolute_descending(const std::string &name) const;

    /// reference to the thing this object is indexing
    const Group *ref;

  protected:
    /// the indices themselves
    std::vector<int> v_index;
  };

  // merge
  // i.e. getting the OR of two or more filter results
  template <typename Idx, typename ...Idxs>
  Idx merge(Idx &index, const Idxs &...indices);
}

#include "Indices.cc"

#endif
