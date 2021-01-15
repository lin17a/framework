// -*- C++ -*-
// author: afiq anuar
// short: handling of the Collection to Dataset linking through TObject's Notify mechanism

#ifndef FWK_ALLOCATOR_H
#define FWK_ALLOCATOR_H

// note: more details in https://root-forum.cern.ch/t/crash-when-reading-too-many-branches/39892/17

#include "Heap.h"
#include "TObject.h"

namespace Framework {
  class Allocator : public TObject {
  public:
    /// provide the allocating function
    template <typename Allocate>
    void set_allocator(Allocate allocate_);

    /// to check allocator validity
    explicit operator bool() const noexcept { return bool(allocate); };

  private:
    /// overriding TObject::Notify
    bool Notify() override;

    /// the fuction to be ran
    std::function<void()> allocate;
  };
}

#include "Allocator.cc"

#endif
