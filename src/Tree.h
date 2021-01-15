// -*- C++ -*-
// author: afiq anuar
// short: an interface for creating output trees from groups

#ifndef FWK_TREE_H
#define FWK_TREE_H

// note: the class relies on Group::reorder so as to avoid copying, which means any branches whose data assumes a particular order of its elements
// note: e.g. GenPart_motherIdx will inevitably break i.e. their data becoming meaningless

#include "Group.h"

#include "TFile.h"
#include "TTree.h"

namespace Framework {
  class Tree {
  public:
    /// no default constructor
    Tree() = delete;

    /// constructor
    /// https://root.cern.ch/doc/master/Compression_8h_source.html
    /// main compression choices are 209 (size 39M, run time 4m 41s) or 505 (size 48M, run time 3m18s)
    /// note the whole run time including analysis (some 2M gen level) is counted in test
    /// default ROOT 6.18 compression is 3m15s or so
    /// test results in a tree with single and array branches with 300k events
    Tree(const std::string &filename, const std::string &treename, int compression = 505);

    /// make single type branches i.e. non-array
    /// given that group is always arrays underneath
    /// this just makes branches for the data[0] elements
    /// without an extra branch for the counter
    template <typename Group, typename ...Attributes>
    bool make_single_branches(Group &group, Attributes &&...attrs);

    /// make branches of the array type
    /// i.e. branches for all requested data
    /// and a branch keeping track of the current element counter
    template <typename Group, typename ...Attributes>
    bool make_array_branches(Group &group, Attributes &&...attrs);

    /// fill the tree - reallocate the branches if needed
    void fill();

    /// save the tree into a ROOT file
    void save(/*const std::string &name*/) const;

  protected:
    // FIXME https://root-forum.cern.ch/t/follow-up-on-ram-vs-disk-resident-ttree-compression-bug/40775
    // for the moment use disk-resident tree as workaround
    std::unique_ptr<TFile> file;

    /// ptr to the tree
    TTree *ptr;

    /// storage of source groups and the branches they contribute
    /// tl;dr in a given tree, each group must be unique (checked by name) and either one of make_x_branches is called for it exactly once
    /// more info: source group stored as functions (with lambda captures) to avoid shenanigans with differing group types
    /// the function is called whenever reallocation is needed
    /// the string is the group name, to ensure that for a given group, only the first call to make_x_branches is registered
    /// this is to make no-op calls like so:
    /// make_single_branches(group, attrs1...); make_array_branches(group, attrs2...);
    /// which, even if technically legit, is gonna be very confusing
    /// similarly, don't do:
    /// make_array_branches(group, attrs1...); make_array_branches(group, attrs2...);
    /// which, even though the counter branch clash is checkable, is more headache than it's worth
    std::vector<std::tuple<std::string, std::function<void()>, std::vector<TBranch *>>> v_branch;
  };
}

#include "Tree.cc"

#endif
