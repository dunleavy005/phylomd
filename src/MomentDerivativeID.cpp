#include "MomentDerivativeID.h"

#include <algorithm>


//
// This source file defines the MomentDerivativeID-related functions.
//


void get_moment_derivative_ids_aux(const Vector<EdgeSet>& esets,
                                   std::size_t curr_es_ind, int max_order,
                                   const MomentDerivativeID& curr_id,
                                   Vector<MomentDerivativeID>& ids) {
  // Cache the current moment/derivative ID.
  ids.push_back(curr_id);

  // Find the next moment/derivative IDs.
  for (std::size_t es_ind = curr_es_ind; es_ind < esets.size(); ++es_ind) {
    for (int order = 1; order <= max_order; ++order) {
      // Form the next moment/derivative ID by concatenating the current
      // moment/derivative ID elements with the next ID element.
      MomentDerivativeID next_id;
      next_id.reserve(curr_id.size() + 1);
      next_id.insert(next_id.end(), curr_id.begin(), curr_id.end());
      next_id.emplace_back(esets[es_ind], order);

      // Recurse over the next possible ID elements.
      get_moment_derivative_ids_aux(esets, es_ind + 1, max_order - order,
                                    next_id, ids);
    }
  }
}


Vector<MomentDerivativeID> get_moment_derivative_ids(
    const Vector<EdgeSet>& esets, int max_order) {
  Vector<MomentDerivativeID> ids;

  // Manually create the "empty" moment/derivative ID (i.e. the partial
  // likelihood ID).
  ids.emplace_back();

  // Start the recursion over the multiple edge sets and orders.
  for (std::size_t es_ind = 0; es_ind < esets.size(); ++es_ind) {
    for (int order = 1; order <= max_order; ++order) {
      get_moment_derivative_ids_aux(
          esets, es_ind + 1, max_order - order,
          {MomentDerivativeIDElement(esets[es_ind], order)}, ids);
    }
  }

  return ids;
}


int find_moment_derivative_id_index(const Vector<MomentDerivativeID>& ids,
                                    const MomentDerivativeID& id) {
  return std::lower_bound(ids.begin(), ids.end(), id) - ids.begin();
}


std::string create_moment_derivative_id_label(const MomentDerivativeID& md_id) {
  std::string md_id_label;

  // Loop through the ID elements and append label information to the output
  // string.
  for (auto it = md_id.begin(); it != md_id.end(); ++it) {
    md_id_label += "(";
    md_id_label += std::to_string(it->set().label() + 1);
    md_id_label += "):";
    md_id_label += std::to_string(it->order());
    if (it != md_id.end() - 1) md_id_label += "--";
  }

  return md_id_label;
}