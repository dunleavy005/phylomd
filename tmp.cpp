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




///////////////////////////////////////////////////////////////////////////////
//////////////////////////// SET PARTITIONS ///////////////////////////////////
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

// [[Rcpp::export]]
int print_part_sets(VectorVector<int> edge_sets) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);

  for (std::size_t i = 0; i < psets.size(); ++i) {
    Rcpp::Rcout << "(";
    for (std::size_t j = 0; j < psets[i].label().size(); ++j) {
      Rcpp::Rcout << psets[i].label()[j];
      if (j < psets[i].label().size() - 1) Rcpp::Rcout << ",";
    }
    Rcpp::Rcout << ") - ";
    for (std::size_t j = 0; j < psets[i].set().size(); ++j) {
      Rcpp::Rcout << psets[i].set()[j];
      if (j < psets[i].set().size() - 1) Rcpp::Rcout << ",";
    }
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

// [[Rcpp::export]]
int print_list_ids(VectorVector<int> edge_sets, int max_order) {
  Vector<PartitionSet> psets = partition_edge_sets(edge_sets);
  Vector<ListID> ids = find_list_ids(psets, max_order);

  for (std::size_t i = 0; i < ids.size(); ++i) {
    Rcpp::Rcout << "<";
    for (std::size_t j = 0; j < ids[i].elems().size(); ++j) {
      Rcpp::Rcout << "[(";
      for (std::size_t k = 0; k < ids[i].elems()[j].pset().label().size(); ++k) {
        Rcpp::Rcout << ids[i].elems()[j].pset().label()[k];
        if (k < ids[i].elems()[j].pset().label().size() - 1) Rcpp::Rcout << ",";
      }
      Rcpp::Rcout << ")-" << ids[i].elems()[j].order() << "]";
      if (j < ids[i].elems().size() - 1) Rcpp::Rcout << " , ";
    }
    Rcpp::Rcout << "> ----- {";
    for (auto it = ids[i].sum_orders().begin(); it != ids[i].sum_orders().end(); ++it) {
      Rcpp::Rcout << "(";
      for (std::size_t j = 0; j < it->first.label().size(); ++j) {
        Rcpp::Rcout << it->first.label()[j];
        if (j < it->first.label().size() - 1) Rcpp::Rcout << ",";
      }
      Rcpp::Rcout << "):" << it->second;
      if (std::distance(it, ids[i].sum_orders().end()) > 1) Rcpp::Rcout << " , ";
    }
    Rcpp::Rcout << "}" << std::endl;
  }

  return ids.size();
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
