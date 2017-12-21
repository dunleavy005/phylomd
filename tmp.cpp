#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


template <typename T>
using Vector = std::vector<T>;

template <typename T>
using VectorVector = std::vector<std::vector<T>>;

template <typename T1, typename T2>
using Map = std::map<T1, T2>;


class PartitionSet {
 private:
  Vector<int> label_;
  Vector<int> set_;

 public:
  PartitionSet(const Vector<int>& label, const Vector<int>& set)
      : label_(label), set_(set) {}

  const Vector<int>& label() const { return label_; }
  const Vector<int>& set() const { return set_; }
};

// This operator is needed for sum-order map indexing within a given list ID.
bool operator<(const PartitionSet& lhs, const PartitionSet& rhs) {
  return lhs.label() < rhs.label();
}
// This operator is needed for finding unique ID elements in list IDs.
bool operator==(const PartitionSet& lhs, const PartitionSet& rhs) {
  return lhs.label() == rhs.label();
}
bool operator>(const PartitionSet& lhs, const PartitionSet& rhs) {
  return rhs < lhs;
}
bool operator<=(const PartitionSet& lhs, const PartitionSet& rhs) {
  return !(rhs < lhs);
}
bool operator>=(const PartitionSet& lhs, const PartitionSet& rhs) {
  return !(lhs < rhs);
}
bool operator!=(const PartitionSet& lhs, const PartitionSet& rhs) {
  return !(lhs == rhs);
}


class IDElement {
 private:
  const PartitionSet* pset_;
  int order_;

 public:
  IDElement(const PartitionSet& pset, int order)
      : pset_(&pset), order_(order) {}

  const PartitionSet& pset() const { return *pset_; }
  int order() const { return order_; }
};

// This operator is needed for sorting list IDs.
bool operator<(const IDElement& lhs, const IDElement& rhs) {
  return std::forward_as_tuple(lhs.pset(), lhs.order()) <
         std::forward_as_tuple(rhs.pset(), rhs.order());
}
// This operator is needed for finding unique ID elements in list IDs.
bool operator==(const IDElement& lhs, const IDElement& rhs) {
  return std::forward_as_tuple(lhs.pset(), lhs.order()) ==
         std::forward_as_tuple(rhs.pset(), rhs.order());
}
bool operator>(const IDElement& lhs, const IDElement& rhs) { return rhs < lhs; }
bool operator<=(const IDElement& lhs, const IDElement& rhs) {
  return !(rhs < lhs);
}
bool operator>=(const IDElement& lhs, const IDElement& rhs) {
  return !(lhs < rhs);
}
bool operator!=(const IDElement& lhs, const IDElement& rhs) {
  return !(lhs == rhs);
}


class ListID {
 private:
  Vector<IDElement> elems_;
  Map<PartitionSet, int> sum_orders_;

 public:
  ListID(const Vector<IDElement>& elems,
         const Map<PartitionSet, int>& sum_orders)
      : elems_(elems), sum_orders_(sum_orders) {
    std::sort(elems_.begin(), elems_.end());
  }

  typedef Vector<IDElement>::const_iterator const_iterator;

  const Vector<IDElement>& elems() const { return elems_; }
  const Map<PartitionSet, int>& sum_orders() const { return sum_orders_; }
};


class EdgeList {
 private:
  const ListID& id_;
  // This is a vector of 3-tuples that summarizes the recurrence relation
  // between this edge list and the node lists.
  //
  // For any given edge, the first tuple element represents a possible
  // assignment of an ID element, the second element contains the corresponding
  // node list index associated with the child node, and the third element
  // stores the related counting coefficient.
  Vector<std::tuple<const IDElement&, int, int>> recursion_info_;

  void init_recursion_info_aux(
      const Vector<ListID>& list_ids, const arma::mat& choose,
      Vector<std::tuple<const IDElement&, int, int>>& recursion_info) const;
  Vector<std::tuple<const IDElement&, int, int>> init_recursion_info(
      const Vector<ListID>& list_ids, const arma::mat& choose) const;

 public:
  EdgeList(const ListID& id, const Vector<ListID>& list_ids,
           const arma::mat& choose)
      : id_(id), recursion_info_(init_recursion_info(list_ids, choose)) {}

  const ListID& id() const { return id_; }
  const Vector<std::tuple<const IDElement&, int, int>>& recursion_info() const {
    return recursion_info_;
  }
};


class NodeList {
 private:
  const ListID& id_;
  // This is a vector of 3-tuples that summarizes the recurrence relation
  // between this node list and the edge lists.
  //
  // For any given node, the first tuple element contains the edge list index
  // associated with the left child edge, the second element represents the edge
  // list index associated with the right child edge, and the third element
  // stores the related counting coefficient.
  Vector<std::tuple<int, int, int>> recursion_info_;

  void init_recursion_info_aux(
      const Vector<ListID>& list_ids, const arma::mat& choose,
      ListID::const_iterator curr_it,
      const Vector<IDElement>& left_elist_id_elems,
      Vector<std::tuple<int, int, int>>& recursion_info) const;
  Vector<std::tuple<int, int, int>> init_recursion_info(
      const Vector<ListID>& list_ids, const arma::mat& choose) const;

 public:
  NodeList(const ListID& id, const Vector<ListID>& list_ids,
           const arma::mat& choose)
      : id_(id), recursion_info_(init_recursion_info(list_ids, choose)) {}

  const ListID& id() const { return id_; }
  const Vector<std::tuple<int, int, int>>& recursion_info() const {
    return recursion_info_;
  }
};


// helper functions

int find_list_id_index(const Vector<ListID>& list_ids,
                       const Vector<IDElement>& id_elems) {
  return std::lower_bound(
             list_ids.begin(), list_ids.end(), id_elems,
             [](const ListID& list_id, const Vector<IDElement>& id_elems) {
               return list_id.elems() < id_elems;
             }) -
         list_ids.begin();
}

ListID::const_iterator find_next_id_elem(ListID::const_iterator begin_it,
                                         ListID::const_iterator end_it,
                                         const IDElement& curr_id_elem) {
  return std::find_if(begin_it, end_it, [&](const IDElement& id_elem) {
    return id_elem != curr_id_elem;
  });
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////// PARTITION SETS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void partition_edge_sets_aux(const VectorVector<int>& edge_sets, int curr_ind,
                             const Vector<int>& curr_label,
                             const Vector<int>& curr_set,
                             Vector<PartitionSet>& partition_sets) {
  // 1) If the current partitioned set is empty, we can exit the function.
  // Otherwise, the current partitioned set is non-empty.
  // 2) Furthermore, if we have traversed over all the edge sets, cache the
  // current partitioned edge set and exit the function.
  if (curr_set.empty()) {
    return;
  } else if ((std::size_t)curr_ind == edge_sets.size()) {
    partition_sets.emplace_back(curr_label, curr_set);
    return;
  }

  // Modify the current partitioned edge set by intersecting/differencing with
  // the current edge set.
  Vector<int> intersect_label;
  intersect_label.reserve(curr_label.size() + 1);
  intersect_label.insert(intersect_label.end(), curr_label.begin(),
                         curr_label.end());
  intersect_label.push_back(curr_ind);

  Vector<int> intersect_set;
  intersect_set.reserve(curr_set.size());
  std::set_intersection(curr_set.begin(), curr_set.end(),
                        edge_sets[curr_ind].begin(), edge_sets[curr_ind].end(),
                        std::back_inserter(intersect_set));

  Vector<int> diff_set;
  diff_set.reserve(curr_set.size());
  std::set_difference(curr_set.begin(), curr_set.end(),
                      edge_sets[curr_ind].begin(), edge_sets[curr_ind].end(),
                      std::back_inserter(diff_set));

  // Recurse over the next edge set.
  ++curr_ind;
  partition_edge_sets_aux(edge_sets, curr_ind, intersect_label, intersect_set,
                          partition_sets);
  partition_edge_sets_aux(edge_sets, curr_ind, curr_label, diff_set,
                          partition_sets);
}

Vector<PartitionSet> partition_edge_sets(VectorVector<int>& edge_sets) {
  // This partitioning algorithm is summarized at
  // https://bosker.wordpress.com/2013/07/10/venn-diagram-partitioning/.

  // Sort the edge sets and incrementally union them together.
  Vector<int> union_set;
  for (auto& es : edge_sets) {
    std::sort(es.begin(), es.end());
    union_set.insert(union_set.end(), es.begin(), es.end());
  }
  std::sort(union_set.begin(), union_set.end());
  auto end_it = std::unique(union_set.begin(), union_set.end());
  union_set.resize(end_it - union_set.begin());

  // Create the partitioned edge sets.
  Vector<PartitionSet> partition_sets;
  partition_sets.reserve(
      std::min((int)std::pow(2, edge_sets.size()) - 1, (int)union_set.size()));
  partition_edge_sets_aux(edge_sets, 0, {}, union_set, partition_sets);

  return partition_sets;
}

void print_partition_set_label(const PartitionSet& pset) {
  Rcpp::Rcout << "(";
  for (std::size_t i = 0; i < pset.label().size(); ++i) {
    Rcpp::Rcout << pset.label()[i];
    if (i < pset.label().size() - 1) Rcpp::Rcout << ",";
  }
  Rcpp::Rcout << ")";
}

void print_partition_set_entries(const PartitionSet& pset) {
  for (std::size_t i = 0; i < pset.set().size(); ++i) {
    Rcpp::Rcout << pset.set()[i];
    if (i < pset.set().size() - 1) Rcpp::Rcout << ",";
  }
}

// [[Rcpp::export]]
int print_partition_sets(VectorVector<int> edge_sets) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);

  for (const auto& pset : psets) {
    print_partition_set_label(pset);
    Rcpp::Rcout << " - ";
    print_partition_set_entries(pset);
    Rcpp::Rcout << std::endl;
  }

  return psets.size();
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// LIST IDS //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void get_list_ids_aux(const Vector<PartitionSet>& psets,
                      std::size_t curr_ps_ind, int max_order,
                      const Vector<IDElement>& curr_elems,
                      const Map<PartitionSet, int>& curr_sum_orders,
                      Vector<ListID>& ids) {
  // Cache the current list ID.
  ids.emplace_back(curr_elems, curr_sum_orders);

  // Find the next list IDs.
  // (Note: the loop body will not be entered if there are no more list IDs to
  // find.)
  for (int order = curr_elems.back().order(); order <= max_order; ++order) {
    for (std::size_t ps_ind = curr_ps_ind; ps_ind < psets.size(); ++ps_ind) {
      // Form the next list ID by concatenating the current list ID elements
      // with the next ID element and creating a new sum-order map by updating
      // the current one.
      Vector<IDElement> next_elems;
      next_elems.reserve(curr_elems.size() + 1);
      next_elems.insert(next_elems.end(), curr_elems.begin(), curr_elems.end());
      next_elems.emplace_back(psets[ps_ind], order);

      Map<PartitionSet, int> next_sum_orders = curr_sum_orders;
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
  ids.emplace_back(Vector<IDElement>(), Map<PartitionSet, int>());

  // Start the recursion over the multiple partition sets and orders.
  for (int order = 1; order <= max_order; ++order) {
    for (std::size_t ps_ind = 0; ps_ind < psets.size(); ++ps_ind) {
      get_list_ids_aux(psets, ps_ind, max_order - order,
                       {IDElement(psets[ps_ind], order)},
                       {{psets[ps_ind], order}}, ids);
    }
  }

  // Sort the list IDs.
  std::sort(ids.begin(), ids.end(),
            [](const ListID& left_id, const ListID& right_id) {
              return left_id.elems() < right_id.elems();
            });

  return ids;
}

void print_id_element(const IDElement& id_elem) {
  Rcpp::Rcout << "[";
  print_partition_set_label(id_elem.pset());
  Rcpp::Rcout << "-" << id_elem.order() << "]";
}

void print_list_id_elems(const ListID& id) {
  Rcpp::Rcout << "< ";
  for (std::size_t i = 0; i < id.elems().size(); ++i) {
    print_id_element(id.elems()[i]);
    if (i < id.elems().size() - 1) Rcpp::Rcout << " , ";
  }
  Rcpp::Rcout << " >";
}

void print_list_id_sum_orders(const ListID& id) {
  Rcpp::Rcout << "{";
  for (auto it = id.sum_orders().begin(); it != id.sum_orders().end(); ++it) {
    print_partition_set_label(it->first);
    Rcpp::Rcout << ":" << it->second;
    if (std::distance(it, id.sum_orders().end()) > 1) Rcpp::Rcout << " , ";
  }
  Rcpp::Rcout << "}";
}

// [[Rcpp::export]]
int print_list_ids(VectorVector<int> edge_sets, int max_order) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);
  Vector<ListID> ids = get_list_ids(psets, max_order);

  for (const auto& id : ids) {
    print_list_id_elems(id);
    Rcpp::Rcout << " ----- ";
    print_list_id_sum_orders(id);
    Rcpp::Rcout << std::endl;
  }

  return ids.size();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// EDGE LISTS //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void EdgeList::init_recursion_info_aux(
    const Vector<ListID>& list_ids, const arma::mat& choose,
    Vector<std::tuple<const IDElement&, int, int>>& recursion_info) const {
  // Loop through the unique ID elements in the edge list ID and cache the
  // corresponding recursion information 3-tuples.
  for (auto curr_it = id_.elems().begin(); curr_it != id_.elems().end();) {
    // Determine the node list ID elements associated with the child node.
    Vector<IDElement> nlist_id_elems;
    nlist_id_elems.reserve(id_.elems().size() - 1);

    for (auto it = id_.elems().begin(); it != id_.elems().end(); ++it) {
      // The node list ID elements should not include the current ID element.
      // (Note: we intentionally compare iterators for equality.)
      if (it != curr_it) nlist_id_elems.push_back(*it);
    }

    // Compute the node list index associated with the child node.
    int nlist_ind = find_list_id_index(list_ids, nlist_id_elems);

    // Cache the current recursion information 3-tuple.
    const IDElement& curr_id_elem = *curr_it++;
    int choose_coef =
        choose(id_.sum_orders().at(curr_id_elem.pset()), curr_id_elem.order());
    recursion_info.emplace_back(curr_id_elem, nlist_ind, choose_coef);

    // Find the next unique ID element in the edge list ID.
    curr_it = find_next_id_elem(curr_it, id_.elems().end(), curr_id_elem);
  }
}

Vector<std::tuple<const IDElement&, int, int>> EdgeList::init_recursion_info(
    const Vector<ListID>& list_ids, const arma::mat& choose) const {
  Vector<std::tuple<const IDElement&, int, int>> recursion_info;
  recursion_info.reserve(id_.elems().size());
  init_recursion_info_aux(list_ids, choose, recursion_info);

  return recursion_info;
}

// [[Rcpp::export]]
int print_elist_recursion_info(VectorVector<int> edge_sets, int max_order,
                               int elist_ind) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);
  Vector<ListID> ids = get_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  EdgeList elist(ids[elist_ind], ids, choose);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(ids[elist_ind]);
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& tup : elist.recursion_info()) {
    const IDElement& id_elem = std::get<0>(tup);
    int nlist_ind = std::get<1>(tup);
    int choose_coef = std::get<2>(tup);

    print_id_element(id_elem);
    Rcpp::Rcout << " ----- ";
    print_list_id_elems(ids[nlist_ind]);
    Rcpp::Rcout << " ----- ";
    Rcpp::Rcout << choose_coef << std::endl;
  }

  return elist.recursion_info().size();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// NODE LISTS //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void NodeList::init_recursion_info_aux(
    const Vector<ListID>& list_ids, const arma::mat& choose,
    ListID::const_iterator curr_it,
    const Vector<IDElement>& left_elist_id_elems,
    Vector<std::tuple<int, int, int>>& recursion_info) const {
  // Determine the edge list ID elements associated with the right child edge.
  Vector<IDElement> right_elist_id_elems;
  right_elist_id_elems.reserve(id_.elems().size() - left_elist_id_elems.size());
  std::set_difference(id_.elems().begin(), id_.elems().end(),
                      left_elist_id_elems.begin(), left_elist_id_elems.end(),
                      std::back_inserter(right_elist_id_elems));

  // Compute the edge list indices associated with both child edges.
  int left_elist_ind = find_list_id_index(list_ids, left_elist_id_elems);
  int right_elist_ind = find_list_id_index(list_ids, right_elist_id_elems);

  // Calculate the related counting coefficient.
  int choose_coef = 1;
  for (auto it = id_.sum_orders().begin(); it != id_.sum_orders().end(); ++it) {
    auto find_it = list_ids[left_elist_ind].sum_orders().find(it->first);
    if (find_it != list_ids[left_elist_ind].sum_orders().end() &&
        find_it->second != it->second) {
      choose_coef *= choose(it->second, find_it->second);
    }
  }

  // Cache the current recursion information 3-tuple.
  recursion_info.emplace_back(left_elist_ind, right_elist_ind, choose_coef);

  // Loop through the remaining unique ID elements in the node list ID.
  for (auto it = curr_it; it != id_.elems().end();) {
    // Determine the next possible edge list ID elements associated with the
    // left child edge.
    Vector<IDElement> next_left_elist_id_elems;
    next_left_elist_id_elems.reserve(left_elist_id_elems.size() + 1);
    next_left_elist_id_elems.insert(next_left_elist_id_elems.end(),
                                    left_elist_id_elems.begin(),
                                    left_elist_id_elems.end());
    const IDElement& next_id_elem = *it++;
    next_left_elist_id_elems.push_back(next_id_elem);

    // Recurse over the previously computed edge list ID elements.
    init_recursion_info_aux(list_ids, choose, it, next_left_elist_id_elems,
                            recursion_info);

    // Find the next unique ID element in the node list ID.
    it = find_next_id_elem(it, id_.elems().end(), next_id_elem);
  }
}

Vector<std::tuple<int, int, int>> NodeList::init_recursion_info(
    const Vector<ListID>& list_ids, const arma::mat& choose) const {
  Vector<std::tuple<int, int, int>> recursion_info;
  recursion_info.reserve(std::pow(2, id_.elems().size()));
  init_recursion_info_aux(list_ids, choose, id_.elems().begin(), {},
                          recursion_info);

  return recursion_info;
}

// [[Rcpp::export]]
int print_nlist_recursion_info(VectorVector<int> edge_sets, int max_order,
                               int nlist_ind) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);
  Vector<ListID> ids = get_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  NodeList nlist(ids[nlist_ind], ids, choose);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(ids[nlist_ind]);
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& tup : nlist.recursion_info()) {
    int left_elist_ind = std::get<0>(tup);
    int right_elist_ind = std::get<1>(tup);
    int choose_coef = std::get<2>(tup);

    print_list_id_elems(ids[left_elist_ind]);
    Rcpp::Rcout << " ----- ";
    print_list_id_elems(ids[right_elist_ind]);
    Rcpp::Rcout << " ----- ";
    Rcpp::Rcout << choose_coef << std::endl;
  }

  return nlist.recursion_info().size();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// DO NOT DELETE! //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//
//void find_key_subsets_aux(std::vector<int>::const_iterator curr_it, std::vector<int>::const_iterator end_it,
//                          const std::vector<int>& curr_subset, std::vector<std::vector<int>>& subsets) {
//  
//  subsets.push_back(curr_subset);
//  if(curr_it == end_it) return;
//  
//  std::vector<int> uniq_nums(end_it - curr_it);
//  std::unique_copy(curr_it, end_it, uniq_nums.begin());
//  
//  std::vector<int> next_subset = curr_subset;
//  next_subset.resize(next_subset.size() + 1);
//  
//  for (int i = 0; uniq_nums[i] != 0 && i < uniq_nums.size(); ++i) {
//    std::vector<int>::const_iterator next_it = std::find(curr_it, end_it, uniq_nums[i]);
//    next_subset.back() = uniq_nums[i];
//    find_key_subsets_aux(++next_it, end_it, next_subset, subsets);
//  }
//}
//
//// [[Rcpp::export]]
//std::vector<std::vector<int>> find_key_subsets(const std::vector<int>& key) {
//  std::vector<std::vector<int>> subsets;
//  find_key_subsets_aux(key.begin(), key.end(), {}, subsets);
//  return subsets;
//}




//void initialize_cache_lists_aux(int order) {
//  std::vector<CacheList> cache_lists;
//  std::vector<std::vector<int>> list_ids;
//  
//  // zero case
//  list_ids.push_back({});
//  cache_lists.emplace_back({}, list_ids);
//  
//  for (int i = 1; i <= order; ++i) {
//    find_cache_list_ids_aux(i, order, 1, list_ids);
//    for (int j = 0; j)
//    cache_lists.emplace_back()
//  }
//}
//
//std::vector<CacheList> initialize_cache_lists(int order) {
//  std::vector<CacheList> cache_lists;
//  initialize_cache_lists_aux(order, cache_lists);
//  
//  
//  for (int i = 0; i < list_ids.size(); ++i) {
//    outp.emplace_back(list_ids[0], list_ids);
//  }
//}
//
//CacheList::CacheList(const std::vector<int>& list_id, const std::vector<std::vector<int>>& list_ids) {
//  
//}
