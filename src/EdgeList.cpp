#include "EdgeList.h"

#include <algorithm>


//
// This source file defines the EdgeList-related member functions.
//


void EdgeList::init_recursion_info_aux(
    const Vector<ListID>& ids, const arma::mat& choose,
    Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info) const {
  // Loop through the unique ID elements in the edge list ID and cache the
  // corresponding recursion information 3-tuples.
  for (auto curr_it = id_.begin(); curr_it != id_.end();) {
    // Determine the node list ID elements associated with the child node.
    Vector<ListIDElement> nlist_id_elems;
    nlist_id_elems.reserve(id_.size() - 1);

    for (auto it = id_.begin(); it != id_.end(); ++it) {
      // The node list ID elements should not include the current ID element.
      // (Note: we intentionally compare iterators for equality.)
      if (it != curr_it) nlist_id_elems.push_back(*it);
    }

    // Compute the node list index associated with the child node.
    int nlist_ind = find_list_id_index(ids, nlist_id_elems);

    // Cache the current recursion information 3-tuple.
    const ListIDElement& curr_id_elem = *curr_it++;
    int choose_coef =
        choose(id_.sum_orders().at(curr_id_elem.set()), curr_id_elem.order());
    recursion_info.emplace_back(curr_id_elem, nlist_ind, choose_coef);

    // Find the next unique ID element in the edge list ID.
    curr_it = find_next_id_elem(curr_it, id_.end(), curr_id_elem);
  }
}


Vector<std::tuple<const ListIDElement&, int, int>>
EdgeList::init_recursion_info(const Vector<ListID>& ids,
                              const arma::mat& choose) const {
  Vector<std::tuple<const ListIDElement&, int, int>> recursion_info;
  recursion_info.reserve(id_.size());
  init_recursion_info_aux(ids, choose, recursion_info);

  return recursion_info;
}


Map<Ref<const PartitionSet>, std::pair<int, int>>
EdgeList::init_recursion_info_inds() const {
  Map<Ref<const PartitionSet>, std::pair<int, int>> recursion_info_inds;

  // Loop through the recursion information 3-tuples and cache the start/end
  // indices associated with each unique partition set.
  for (auto begin_it = recursion_info_.begin();
       begin_it != recursion_info_.end();) {
    // Extract the unique partition set.
    const PartitionSet& curr_pset = std::get<0>(*begin_it).set();

    // Find the next unique partition set.
    // (Note: the current unique partition set ends where the next unique
    // partition set starts.)
    auto end_it = std::find_if(
        begin_it, recursion_info_.end(),
        [&](const std::tuple<const ListIDElement&, int, int>& ri_tup) {
          return std::get<0>(ri_tup).set() != curr_pset;
        });

    // Cache the start/end indices associated with the current unique partition
    // set.
    int begin_ind = begin_it - recursion_info_.begin();
    int end_ind = end_it - recursion_info_.begin();
    recursion_info_inds.emplace(curr_pset,
                                std::pair<int, int>({begin_ind, end_ind}));

    // Prepare for the next loop iteration.
    begin_it = end_it;
  }

  return recursion_info_inds;
}