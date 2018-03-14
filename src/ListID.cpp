#include "ListID.h"

#include <algorithm>


//
// This source file defines the ListID-related functions.
//


void get_list_ids_aux(const Vector<PartitionSet>& psets,
                      std::size_t curr_ps_ind, int max_order,
                      const Vector<ListIDElement>& curr_elems,
                      const Map<Ref<const PartitionSet>, int>& curr_sum_orders,
                      Vector<ListID>& ids) {
  // Cache the current list ID.
  ids.emplace_back(curr_elems, curr_sum_orders);

  // Find the next list IDs.
  for (std::size_t ps_ind = curr_ps_ind; ps_ind < psets.size(); ++ps_ind) {
    // Will the next list ID elements use the current partition set?
    // If so, we iterate over a subset of orders.
    // Otherwise, we iterate over all possible orders.
    int begin_order = (ps_ind == curr_ps_ind) ? curr_elems.back().order() : 1;

    for (int order = begin_order; order <= max_order; ++order) {
      // Form the next list ID by concatenating the current list ID elements
      // with the next ID element and creating a new sum-order map by updating
      // the current one.
      Vector<ListIDElement> next_elems;
      next_elems.reserve(curr_elems.size() + 1);
      next_elems.insert(next_elems.end(), curr_elems.begin(), curr_elems.end());
      next_elems.emplace_back(psets[ps_ind], order);

      Map<Ref<const PartitionSet>, int> next_sum_orders = curr_sum_orders;
      auto insert_results = next_sum_orders.emplace(psets[ps_ind], order);
      if (!insert_results.second) insert_results.first->second += order;

      // Recurse over the next possible ID elements.
      get_list_ids_aux(psets, ps_ind, max_order - order, next_elems,
                       next_sum_orders, ids);
    }
  }
}


Vector<ListID> get_list_ids(const Vector<PartitionSet>& psets, int max_order) {
  Vector<ListID> ids;

  // Manually create the "empty" list ID (i.e. the partial likelihood list ID).
  ids.emplace_back(Vector<ListIDElement>(),
                   Map<Ref<const PartitionSet>, int>());

  // Start the recursion over the multiple partition sets and orders.
  for (std::size_t ps_ind = 0; ps_ind < psets.size(); ++ps_ind) {
    for (int order = 1; order <= max_order; ++order) {
      get_list_ids_aux(psets, ps_ind, max_order - order,
                       {ListIDElement(psets[ps_ind], order)},
                       {{psets[ps_ind], order}}, ids);
    }
  }

  return ids;
}


int find_list_id_index(const Vector<ListID>& ids,
                       const Vector<ListIDElement>& id_elems) {
  return std::lower_bound(
             ids.begin(), ids.end(), id_elems,
             [](const ListID& id, const Vector<ListIDElement>& id_elems) {
               return id.elems() < id_elems;
             }) -
         ids.begin();
}


ListID::const_iterator find_next_id_elem(ListID::const_iterator begin_it,
                                         ListID::const_iterator end_it,
                                         const ListIDElement& curr_id_elem) {
  return std::find_if(begin_it, end_it, [&](const ListIDElement& id_elem) {
    return id_elem != curr_id_elem;
  });
}