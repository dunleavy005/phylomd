#ifndef PHYLOMD_PARTITIONSET_H_
#define PHYLOMD_PARTITIONSET_H_

#include <cstddef>

#include "EdgeSet.h"
#include "phylomd_types.h"


//
// This header file defines the PartitionSet class and declares the
// PartitionSet-related functions.
//


class PartitionSet {
 private:
  Vector<int> label_;
  Vector<int> elems_;

 public:
  PartitionSet(const Vector<int>& label, const Vector<int>& elems)
      : label_(label), elems_(elems) {}

  const Vector<int>& label() const { return label_; }
  const Vector<int>& elems() const { return elems_; }
};


inline bool operator<(const PartitionSet& lhs, const PartitionSet& rhs) {
  return lhs.label() < rhs.label();
}
inline bool operator==(const PartitionSet& lhs, const PartitionSet& rhs) {
  return lhs.label() == rhs.label();
}
inline bool operator>(const PartitionSet& lhs, const PartitionSet& rhs) {
  return rhs < lhs;
}
inline bool operator<=(const PartitionSet& lhs, const PartitionSet& rhs) {
  return !(rhs < lhs);
}
inline bool operator>=(const PartitionSet& lhs, const PartitionSet& rhs) {
  return !(lhs < rhs);
}
inline bool operator!=(const PartitionSet& lhs, const PartitionSet& rhs) {
  return !(lhs == rhs);
}


void partition_edge_sets_aux(const Vector<EdgeSet>& esets, std::size_t curr_ind,
                             const Vector<int>& curr_label,
                             const Vector<int>& curr_elems,
                             Vector<PartitionSet>& psets);


Vector<PartitionSet> partition_edge_sets(const Vector<EdgeSet>& esets);


Map<int, const PartitionSet&> create_edge_pset_map(
    const Vector<PartitionSet>& psets);


#endif  // PHYLOMD_PARTITIONSET_H_