// -*- C++ -*-
// author: afiq anuar
// short: handling of datasets i.e. sets of files to be treated as single units

#ifndef FWK_DATASET_H
#define FWK_DATASET_H

#include "Allocator.h"
#include "TTree.h"
#include "TChain.h"

#include <iostream>

// accepted template types are:
// TChain for flat ROOT files analysis
// TTree for column txt csv etc files analysis
// any other type are currently not implemented

namespace Framework {
  template <typename Tree>
  class Dataset {
  public:
    /// constructor
    Dataset(const std::string &name_, const std::string &tree_name_, 
            const std::string &tree_struct_ = "", char tree_delim_ = ' ', 
            const std::vector<std::string> &v_file_ = {});

    /// destructor
    ~Dataset();

    /// setter methods
    void set_tree(const std::string &tree_name_);

    void set_structure(const std::string &tree_struct_);

    void set_delimiter(char tree_delim_);

    /// force_replace is when one wants to entirely replace the file list
    /// if false, and v_file is already filled, the method does nothing
    /// if true, tree is reset and then readded with the new set of files
    /// nfile is how many files one reads from the given vector, defaults to all
    void set_files(const std::vector<std::string> &v_file_, int nfile = -1, const bool force_replace = false);

    /// add a new file to v_file
    /// true if successful i.e. no duplicate files
    /// duplicate is defined by a simple string comparison of the filenames
    /// false when there is a duplicate file, nothing happens
    bool add_file(const std::string &file);

    /// return the number of currently held filenames
    int n_files();

    /// add a new weight to v_weight
    /// true if successful i.e. no weight already has that name
    /// false otherwise, so weight is ignored
    /// TODO include methods that sets weight based on tree content
    bool add_weight(const std::string &wgt_name, const double wgt);

    /// split a dataset into two parts
    /// also renames itself to indicate uniqueness
    /// splitting is file based, so returned dataset will be empty if <= 1 file is held
    /// both datasets will keep the same set of weights
    Dataset split();

    /// getter methods
    long long current_entry(long long entry) const;

    const std::unique_ptr<Tree>& tree() const;

    double get_weight(const std::string &wgt_name) const;

    /// add the files to the Tree and evaluate entries
    /// argument index can be 0 for evaluate everything, or -1 for evaluate only the last v_file element
    void evaluate(int index = 0);

    /// take all the Collections to associate to the tree and allocate resources
    template <typename ...Collections>
    void associate(Collections &...colls);

    /// provide the event analysis function
    /// signature: one argument corresponding to entry number and no return value
    /// the function should contain the relevant analysis steps for one event
    template <typename Analyzer>
    void set_analyzer(Analyzer analyzer_);

    /// perform the analysis
    /// can also cap the total events ran, or skip some
    void analyze(long long total = -1LL, long long skip = -1LL) const;

    /// reset Tree state, but keep the info strings
    void reset();

    /// name of the dataset
    std::string name;

  private:
    /// name of the tree
    std::string tree_name;

    /// tree structure and delimiter when reading text files
    /// ignored when the Tree is TChain
    std::string tree_struct;

    char tree_delim;

    /// filenames
    std::vector<std::string> v_file;

    /// ptr to the tree
    std::unique_ptr<Tree> tree_ptr;

    /// allocator function to be ran at each file change
    Allocator allocator;

    /// event analyzing function
    /// see Dataset::set_analyzer above for more info
    std::function<void(long long)> analyzer;

    /// weights associated to the dataset
    /// mainly in view of MC samples: xsec and such
    std::vector<std::pair<std::string, double>> v_weight;
  };
}

#include "Dataset.cc"

#endif
