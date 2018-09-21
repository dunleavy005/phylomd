#include "EdgeSet.h"

#include <cstddef>


//
// This source file defines the EdgeSet-related functions.
//


Vector<EdgeSet> create_edge_sets(const VectorVector<int>& esets_inp) {
  Vector<EdgeSet> esets;
  esets.reserve(esets_inp.size());

  for (std::size_t i = 0; i < esets_inp.size(); ++i) {
    esets.emplace_back(i, esets_inp[i]);
  }

  return esets;
}