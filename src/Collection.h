// -*- C++ -*-
// author: afiq anuar
// short: collections are groups whose input come from ROOT files - i.e. it is the group-dataset interface

#ifndef FWK_COLLECTION_H
#define FWK_COLLECTION_H

#include "Group.h"
#include "Dataset.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace Framework {
  template <typename ...Ts>
  class Collection : public Group<Ts...> {

  public:
    /// constructor
    /// no default constructor
    Collection() = delete;

    /// this is for single element collections
    /// reserve_ is an estimate of the number of attributes the collection will hold
    Collection(const std::string &name_, int reserve_);

    /// this is for array collections
    /// counter_branch is the presumed branch in the dataset that specifies element count
    /// reserve_ is an estimate of the number of attributes the collection will hold
    /// init is an estimate for the number of elements the collection will hold
    /// 4 is just a "reasonable default"; the code reallocates if capacity is hit
    /// that said, reallocating too often is not good for performance 
    Collection(const std::string &name_, const std::string &counter_name_, int reserve_, int init = 4);

    /// destructor
    ~Collection();

    /// reserve the space for expected number of attributes
    void reserve(int attr);

    /// dummy function so as to not break backward compatibility
    /// FIXME remove when people are done with it
    /// FIXME Lina, Charlotte and who else?
    template <typename Number>
    bool add_attribute(const std::string &attr, const std::string &branch, Number);

    /// add an attribute into the collection
    /// returns false upon failure to add the attribute
    /// this can happen if the data type is inconsistent with the collection
    /// returns true upon a successful registration
    bool add_attribute(const std::string &attr, const std::string &branch);

    /// transform a group of internal attributes into another attribute
    /// the transformation is done element-wise on every element of held data
    /// as this adds a new attribute, its name has to be unique
    template <typename Function, typename ...Attributes>
    bool transform_attribute(const std::string &attr, Function function, Attributes &&...attrs);

    /// associate the attributes to relevant branches in a Dataset
    template <typename Tree>
    void associate(Dataset<Tree> &dataset);

    /// reassociate the attributes to relevant branches in a Dataset
    /// ie to be ran when a reallocation is needed due to capacity hit
    void reassociate();

    /// populate the data with information read from the branches
    void populate(long long entry) override;

  protected:
    /// detach the branches
    void detach();

    /// TBranch isn't enough for realloc, so Collection needs to remember its tree
    TTree* tree;

    /// branch name which gets the counter
    std::string counter_name;

    /// branch associated to the counter
    TBranch *counter_branch;

    /// attribute branches
    std::vector<std::pair<std::string, TBranch *>> v_branch;
  };
}

#include "Collection.cc"

#endif
