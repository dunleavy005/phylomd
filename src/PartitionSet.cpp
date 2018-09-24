#include "PartitionSet.h"

#include <algorithm>
#include <cmath>
#include <iterator>


//
// This source file defines the PartitionSet-related functions.
//


void partition_edge_sets_aux(const Vector<EdgeSet>& esets, std::size_t curr_ind,
                             const Vector<int>& curr_label,
                             const Vector<int>& curr_elems,
                             Vector<PartitionSet>& psets) {
  // If the current partitioned set is empty, we can exit the function.
  // If we have traversed over all the edge sets, cache the current partitioned
  // edge set and exit the function.
  if (curr_elems.empty()) return;
  if (curr_ind == esets.size()) {
    psets.emplace_back(curr_label, curr_elems);
    return;
  }

  // Modify the current partitioned edge set by intersecting/differencing with
  // the current edge set.
  Vector<int> intersect_label;
  intersect_label.reserve(curr_label.size() + 1);
  intersect_label.insert(intersect_label.end(), curr_label.begin(),
                         curr_label.end());
  intersect_label.push_back(esets[curr_ind].label());

  Vector<int> intersect_elems;
  intersect_elems.reserve(curr_elems.size());
  std::set_intersection(
      curr_elems.begin(), curr_elems.end(), esets[curr_ind].elems().begin(),
      esets[curr_ind].elems().end(), std::back_inserter(intersect_elems));

  Vector<int> diff_elems;
  diff_elems.reserve(curr_elems.size() - intersect_elems.size());
  std::set_difference(
      curr_elems.begin(), curr_elems.end(), esets[curr_ind].elems().begin(),
      esets[curr_ind].elems().end(), std::back_inserter(diff_elems));

  // Recurse over the next edge set.
  curr_ind += 1;
  partition_edge_sets_aux(esets, curr_ind, intersect_label, intersect_elems,
                          psets);
  partition_edge_sets_aux(esets, curr_ind, curr_label, diff_elems, psets);
}


Vector<PartitionSet> partition_edge_sets(const Vector<EdgeSet>& esets) {
  // This partitioning algorithm is summarized at
  // https://bosker.wordpress.com/2013/07/10/venn-diagram-partitioning/.

  // Union all the edge set elements.
  Vector<int> union_elems;
  for (const auto& eset : esets) {
    union_elems.insert(union_elems.end(), eset.elems().begin(),
                       eset.elems().end());
  }
  std::sort(union_elems.begin(), union_elems.end());
  auto end_it = std::unique(union_elems.begin(), union_elems.end());
  union_elems.resize(end_it - union_elems.begin());

  // Create the partitioned edge sets.
  Vector<PartitionSet> psets;
  psets.reserve(
      std::min((int)std::pow(2, esets.size()) - 1, (int)union_elems.size()));
  partition_edge_sets_aux(esets, 0, {}, union_elems, psets);
  std::sort(psets.begin(), psets.end());

  return psets;
}


Map<int, const PartitionSet&> create_edge_pset_map(
    const Vector<PartitionSet>& psets) {
  Map<int, const PartitionSet&> edge_psets;

  // Loop through the partition sets and incrementally add (`edge_ind`, `pset`)
  // pairs to the output map.
  for (const auto& pset : psets) {
    for (auto edge_ind : pset.elems()) {
      edge_psets.emplace(edge_ind, pset);
    }
  }

  return edge_psets;
}