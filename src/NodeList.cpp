#include "NodeList.h"

#include <algorithm>
#include <cmath>
#include <iterator>


//
// This source file defines the NodeList-related member functions.
//


void NodeList::init_recursion_info_aux(
    const Vector<ListID>& ids, const arma::mat& choose,
    ListID::const_iterator curr_it,
    const Vector<ListIDElement>& left_elist_id_elems,
    Vector<std::tuple<int, int, int>>& recursion_info) const {
  // Determine the edge list ID elements associated with the right child edge.
  Vector<ListIDElement> right_elist_id_elems;
  right_elist_id_elems.reserve(id_.size() - left_elist_id_elems.size());
  std::set_difference(id_.begin(), id_.end(), left_elist_id_elems.begin(),
                      left_elist_id_elems.end(),
                      std::back_inserter(right_elist_id_elems));

  // Compute the edge list indices associated with both child edges.
  int left_elist_ind = find_list_id_index(ids, left_elist_id_elems);
  int right_elist_ind = find_list_id_index(ids, right_elist_id_elems);

  // Calculate the related counting coefficient.
  int choose_coef = 1;
  for (auto it = id_.sum_orders().begin(); it != id_.sum_orders().end(); ++it) {
    auto find_it = ids[left_elist_ind].sum_orders().find(it->first);
    if (find_it != ids[left_elist_ind].sum_orders().end() &&
        find_it->second != it->second) {
      choose_coef *= choose(it->second, find_it->second);
    }
  }

  // Cache the current recursion information 3-tuple.
  recursion_info.emplace_back(left_elist_ind, right_elist_ind, choose_coef);

  // Loop through the remaining unique ID elements in the node list ID.
  for (auto it = curr_it; it != id_.end();) {
    // Determine the next possible edge list ID elements associated with the
    // left child edge.
    Vector<ListIDElement> next_left_elist_id_elems;
    next_left_elist_id_elems.reserve(left_elist_id_elems.size() + 1);
    next_left_elist_id_elems.insert(next_left_elist_id_elems.end(),
                                    left_elist_id_elems.begin(),
                                    left_elist_id_elems.end());
    const ListIDElement& curr_id_elem = *it++;
    next_left_elist_id_elems.push_back(curr_id_elem);

    // Recurse over the previously computed edge list ID elements.
    init_recursion_info_aux(ids, choose, it, next_left_elist_id_elems,
                            recursion_info);

    // Find the next unique ID element in the node list ID.
    it = find_next_id_elem(it, id_.end(), curr_id_elem);
  }
}


Vector<std::tuple<int, int, int>> NodeList::init_recursion_info(
    const Vector<ListID>& ids, const arma::mat& choose) const {
  Vector<std::tuple<int, int, int>> recursion_info;
  recursion_info.reserve(std::pow(2, id_.size()));
  init_recursion_info_aux(ids, choose, id_.begin(), {}, recursion_info);

  return recursion_info;
}