#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


template <typename T>
using Vector = std::vector<T>;

template <typename T>
using VectorVector = std::vector<std::vector<T>>;


class PartitionSet {
 private:
  Vector<int> label_;
  Vector<int> set_;

 public:
  PartitionSet(const Vector<int>& label, const Vector<int>& set)
      : label_(label), set_(set) {}

  // This operator is needed for sum-order map indexing within a given list ID.
  bool operator<(const PartitionSet& rhs) const { return label_ < rhs.label_; }
  // This operator is needed for finding unique ID elements in edge list IDs.
  bool operator==(const PartitionSet& rhs) const { return label_ == rhs.label_; }

  const Vector<int>& label() const { return label_; }
  const Vector<int>& set() const { return set_; }
};


class IDElement {
 private:
  // This class does not "own" the partition set resources.
  const PartitionSet* pset_;
  int order_;

 public:
  IDElement(const PartitionSet* pset, int order) : pset_(pset), order_(order) {}

  // This operator is needed for sorting list IDs.
  bool operator<(const IDElement& rhs) const {
    return std::tie(*pset_, order_) < std::tie(*rhs.pset_, rhs.order_);
  }
  // This operator is needed for finding unique ID elements in edge list IDs.
  bool operator==(const IDElement& rhs) const {
    return std::tie(*pset_, order_) == std::tie(*rhs.pset_, rhs.order_);
  }

  const PartitionSet& pset() const { return *pset_; }
  int order() const { return order_; }
};


class ListID {
 private:
  Vector<IDElement> elems_;
  std::map<PartitionSet, int> sum_orders_;

 public:
  ListID(const Vector<IDElement>& elems,
         const std::map<PartitionSet, int>& sum_orders)
      : elems_(elems), sum_orders_(sum_orders) {
    std::sort(elems_.begin(), elems_.end());
  }

  // This operator is needed for sorting list IDs.
  bool operator<(const ListID& rhs) const { return elems_ < rhs.elems_; }

  const Vector<IDElement>& elems() const { return elems_; }
  const std::map<PartitionSet, int>& sum_orders() const { return sum_orders_; }
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

  Vector<std::tuple<int, int, int>> init_recursion_info(
      const ListID& id, const Vector<ListID>& list_ids,
      const arma::mat& choose) const;

 public:
  NodeList(const ListID& id, const Vector<ListID>& list_ids,
           const arma::mat& choose)
      : id_(id), recursion_info_(init_recursion_info(id, list_ids, choose)) {}

  const ListID& id() const { return id_; }
  const Vector<std::tuple<int, int, int>>& recursion_info() const {
    return recursion_info_;
  }
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

  Vector<std::tuple<const IDElement&, int, int>> init_recursion_info(
      const ListID& id, const Vector<ListID>& list_ids,
      const arma::mat& choose) const;

 public:
  EdgeList(const ListID& id, const Vector<ListID>& list_ids,
           const arma::mat& choose)
      : id_(id), recursion_info_(init_recursion_info(id, list_ids, choose)) {}

  const ListID& id() const { return id_; }
  const Vector<std::tuple<const IDElement&, int, int>>& recursion_info() const {
    return recursion_info_;
  }
};



///////////////////////////////////////////////////////////////////////////////
//////////////////////////// PARTITION SETS ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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
  intersect_label.insert(intersect_label.end(), curr_label.begin(), curr_label.end());
  intersect_label.emplace_back(curr_ind);

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
  partition_edge_sets_aux(edge_sets, curr_ind, intersect_label, intersect_set, partition_sets);
  partition_edge_sets_aux(edge_sets, curr_ind, curr_label, diff_set, partition_sets);
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
  partition_sets.reserve(std::pow(2, edge_sets.size()) - 1);
  partition_edge_sets_aux(edge_sets, 0, {}, union_set, partition_sets);

  return partition_sets;
}

void print_partition_set_label(const Vector<int>& label) {
  Rcpp::Rcout << "(";
  for (std::size_t i = 0; i < label.size(); ++i) {
    Rcpp::Rcout << label[i];
    if (i < label.size() - 1) Rcpp::Rcout << ",";
  }
  Rcpp::Rcout << ")";
}

void print_partition_set_entries(const Vector<int>& set) {
  for (std::size_t i = 0; i < set.size(); ++i) {
    Rcpp::Rcout << set[i];
    if (i < set.size() - 1) Rcpp::Rcout << ",";
  }
}

// [[Rcpp::export]]
int print_partition_sets(VectorVector<int> edge_sets) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);

  for (const auto& pset : psets) {
    print_partition_set_label(pset.label());
    Rcpp::Rcout << " - ";
    print_partition_set_entries(pset.set());
    Rcpp::Rcout << std::endl;
  }

  return psets.size();
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// LIST IDS /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void find_list_ids_aux(const Vector<PartitionSet>& psets,
                       std::size_t curr_ps_ind, int max_order,
                       const Vector<IDElement>& curr_elems,
                       const std::map<PartitionSet, int>& curr_sum_orders,
                       Vector<ListID>& ids) {
  // Cache the current list ID.
  ids.emplace_back(curr_elems, curr_sum_orders);

  // Find the next list IDs.
  // (Note: the loop body will not be entered if there are no more list IDs to find.)
  for (int order = curr_elems.back().order(); order <= max_order; ++order) {
    for (std::size_t ps_ind = curr_ps_ind; ps_ind < psets.size(); ++ps_ind) {
      // Form the next list ID by concatenating the current list ID elements
      // with the next ID element and creating a new sum-order map by updating
      // the current one.
      Vector<IDElement> next_elems;
      next_elems.reserve(curr_elems.size() + 1);
      next_elems.insert(next_elems.end(), curr_elems.begin(), curr_elems.end());
      next_elems.emplace_back(&psets[ps_ind], order);

      std::map<PartitionSet, int> next_sum_orders = curr_sum_orders;
      auto insert_results = next_sum_orders.emplace(psets[ps_ind], order);
      if (!insert_results.second) insert_results.first->second += order;

      // Recurse over the next possible ID elements.
      find_list_ids_aux(psets, ps_ind, max_order - order, next_elems, next_sum_orders, ids);
    }
  }
}

Vector<ListID> find_list_ids(const Vector<PartitionSet>& psets, int max_order) {
  Vector<ListID> ids;

  // Manually create the "empty" list ID (i.e. the partial likelihood list ID).
  ids.emplace_back(Vector<IDElement>(), std::map<PartitionSet, int>());

  // Start the recursion over the multiple partition sets and orders.
  for (int order = 1; order <= max_order; ++order) {
    for (std::size_t ps_ind = 0; ps_ind < psets.size(); ++ps_ind) {
      find_list_ids_aux(psets, ps_ind, max_order - order,
                        {IDElement(&psets[ps_ind], order)},
                        {{psets[ps_ind], order}}, ids);
    }
  }

  // Sort the list IDs.
  std::sort(ids.begin(), ids.end());

  return ids;
}

void print_id_element(const IDElement& id_elem) {
  Rcpp::Rcout << "[";
  print_partition_set_label(id_elem.pset().label());
  Rcpp::Rcout << "-" << id_elem.order() << "]";
}

void print_list_id_elems(const Vector<IDElement>& id_elems) {
  Rcpp::Rcout << "< ";
  for (std::size_t i = 0; i < id_elems.size(); ++i) {
    print_id_element(id_elems[i]);
    if (i < id_elems.size() - 1) Rcpp::Rcout << " , ";
  }
  Rcpp::Rcout << " >";
}

void print_list_id_sum_orders(const std::map<PartitionSet, int>& sum_orders) {
  Rcpp::Rcout << "{";
  for (auto it = sum_orders.begin(); it != sum_orders.end(); ++it) {
    print_partition_set_label(it->first.label());
    Rcpp::Rcout << ":" << it->second;
    if (std::distance(it, sum_orders.end()) > 1) Rcpp::Rcout << " , ";
  }
  Rcpp::Rcout << "}";
}

// [[Rcpp::export]]
int print_list_ids(VectorVector<int> edge_sets, int max_order) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);
  Vector<ListID> ids = find_list_ids(psets, max_order);

  for (const auto& id : ids) {
    print_list_id_elems(id.elems());
    Rcpp::Rcout << " ----- ";
    print_list_id_sum_orders(id.sum_orders());
    Rcpp::Rcout << std::endl;
  }

  return ids.size();
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// EDGE LISTS /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Vector<std::tuple<const IDElement&, int, int>> EdgeList::init_recursion_info(
    const ListID& id, const Vector<ListID>& list_ids,
    const arma::mat& choose) const {
  Vector<std::tuple<const IDElement&, int, int>> recursion_info;
  recursion_info.reserve(id.elems().size());

  // Loop through the unique ID elements in the edge list ID and cache the
  // corresponding recursion information 3-tuples.
  for (auto curr_it = id.elems().begin(); curr_it != id.elems().end();) {
    // Determine the node list ID elements associated with the child node.
    Vector<IDElement> node_list_id_elems;
    node_list_id_elems.reserve(id.elems().size() - 1);

    for (auto it = id.elems().begin(); it != id.elems().end(); ++it) {
      // The node list ID elements should not include the current ID element.
      // (Note: we intentionally compare iterators for equality.)
      if (it != curr_it) node_list_id_elems.emplace_back(&it->pset(), it->order());
    }

    // Compute the node list index associated with the child node.
    auto node_list_id_it = std::lower_bound(
        list_ids.begin(), list_ids.end(), node_list_id_elems,
        [](const ListID& left_id, const Vector<IDElement>& right_id_elems) {
          return left_id.elems() < right_id_elems;
        });
    int node_list_ind = node_list_id_it - list_ids.begin();

    // Cache the current recursion information 3-tuple.
    const IDElement& curr_id_elem = *curr_it;
    int choose_coef = choose(id.sum_orders().at(curr_it->pset()), curr_it->order());
    recursion_info.emplace_back(curr_id_elem, node_list_ind, choose_coef);

    // Find the next unique ID element in the edge list ID.
    curr_it = std::find_if_not(++curr_it, id.elems().end(),
                               [&curr_id_elem](const IDElement& id_elem) {
                                 return id_elem == curr_id_elem;
                               });
  }

  return recursion_info;
}

// [[Rcpp::export]]
int print_elist_recursion_info(VectorVector<int> edge_sets, int max_order,
                               int elist_ind) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);
  Vector<ListID> ids = find_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  EdgeList elist(ids[elist_ind], ids, choose);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(ids[elist_ind].elems());
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& tup : elist.recursion_info()) {
    const IDElement& id_elem = std::get<0>(tup);
    int node_list_ind = std::get<1>(tup);
    int choose_coef = std::get<2>(tup);

    print_id_element(id_elem);
    Rcpp::Rcout << " ----- ";
    print_list_id_elems(ids[node_list_ind].elems());
    Rcpp::Rcout << " ----- ";
    Rcpp::Rcout << choose_coef << std::endl;
  }

  return elist.recursion_info().size();
}



///////////////////////////////////////////////////////////////////////////////
////////////////////////////// DO NOT DELETE! /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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
