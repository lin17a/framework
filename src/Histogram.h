// -*- C++ -*-
// author: afiq anuar
// short: an interface for creating, filling and saving of histograms from groups

#ifndef FWK_HISTOGRAM_H
#define FWK_HISTOGRAM_H

// note: clone() from another Histogram was attempted, 
// note: abandoned since there was no obvious way to clone the filler accordingly
// note: since that requires accessing copying lambdas and updating its captures

#include "Heap.h"
#include "misc/string_io.h"

#include "TFile.h"

#include "TH1.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"

#include "TH2.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TH2D.h"

#include "TH3.h"
#include "TH3I.h"
#include "TH3F.h"
#include "TH3D.h"

namespace Framework {
  class Histogram {
  public:
    using histfunc = std::pair<std::unique_ptr<TH1>, std::function<void()>>;

    /// constructor
    Histogram();

    /// provide the weight function
    /// signature: no argument and returns a double for the weight
    template <typename Weighter>
    void set_weighter(Weighter weighter_);

    /// make a histogram and its filling function
    /// the filling function takes two arguments
    /// a pointer to the histogram, and the associated weight to be filled
    /// both of which are handled by this class
    template <typename Hist, typename Filler, typename ...Args>
    bool make_histogram(Filler filler, const std::string &name, Args &&...args);

    /// compute the weight and fill all held histograms
    void fill();

    /// unroll 2D or 3D histograms into 1D, starting from the x axis
    /// the output will always be a TH1D, not really worth the hassle to preserve the X in THNX
    /// since it's likely the main histogram type will move away from ROOT
    /// bool argument is to recover the uoflows during the unrolling
    /// unrolled histogram will not have uoflows
    void unroll(bool add_under_overflow = true);

    /// save all held histograms into a ROOT file
    void save_as(const std::string &name) const;

    /// provide reference to held histograms
    const std::vector<histfunc>& histograms() const;

  protected:
    /// the weight to be used when filling the histograms
    double weight;

    /// how to compute the weights
    std::function<double()> weighter;

    /// all histograms and its filling function
    std::vector<histfunc> v_hist;
  };

  template <typename ...Hists>
  void save_all_as(const std::string &name, const Hists &...hists);
}

#include "Histogram.cc"

#endif
