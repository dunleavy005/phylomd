#ifndef PHYLOMD_EDGESET_H_
#define PHYLOMD_EDGESET_H_

#include <algorithm>
#include <utility>

#include "phylomd_types.h"


//
// This header file defines the EdgeSet class and declares the EdgeSet-related
// functions.
//


class EdgeSet {
 private:
  int label_;
  Vector<int> elems_;

 public:
  EdgeSet(int label, Vector<int>&& elems)
      : label_(label), elems_(std::move(elems)) {
    for (auto& elem : elems_) elem -= 1;
    std::sort(elems_.begin(), elems_.end());
  }

  int label() const { return label_; }
  const Vector<int>& elems() const { return elems_; }
};


inline bool operator<(const EdgeSet& lhs, const EdgeSet& rhs) {
  return lhs.label() < rhs.label();
}
inline bool operator==(const EdgeSet& lhs, const EdgeSet& rhs) {
  return lhs.label() == rhs.label();
}
inline bool operator>(const EdgeSet& lhs, const EdgeSet& rhs) {
  return rhs < lhs;
}
inline bool operator<=(const EdgeSet& lhs, const EdgeSet& rhs) {
  return !(rhs < lhs);
}
inline bool operator>=(const EdgeSet& lhs, const EdgeSet& rhs) {
  return !(lhs < rhs);
}
inline bool operator!=(const EdgeSet& lhs, const EdgeSet& rhs) {
  return !(lhs == rhs);
}


Vector<EdgeSet> create_edge_sets(VectorVector<int>& esets_inp);


#endif  // PHYLOMD_EDGESET_H_