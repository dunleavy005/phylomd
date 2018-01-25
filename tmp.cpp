#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


template <typename T>
using Vector = std::vector<T>;

template <typename T>
using VectorVector = std::vector<std::vector<T>>;

template <typename T1, typename T2>
using Map = std::map<T1, T2>;

template <typename T>
using Ref = std::reference_wrapper<T>;


class EdgeSet {
 private:
  int label_;
  Vector<int> elems_;

 public:
  EdgeSet(int label, Vector<int>&& elems)
      : label_(label), elems_(std::move(elems)) {
    std::sort(elems_.begin(), elems_.end());
  }

  int label() const { return label_; }
  const Vector<int>& elems() const { return elems_; }
};


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


template <typename Set>
class IDElement {
 private:
  Ref<const Set> set_;
  int order_;

 public:
  IDElement(const Set& set, int order) : set_(set), order_(order) {}

  const Set& set() const { return set_; }
  int order() const { return order_; }
};

typedef IDElement<PartitionSet> ListIDElement;

// This operator is needed for sorting list IDs.
template <typename Set>
bool operator<(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return std::forward_as_tuple(lhs.set(), lhs.order()) <
         std::forward_as_tuple(rhs.set(), rhs.order());
}
// This operator is needed for finding unique ID elements in list IDs.
template <typename Set>
bool operator==(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return std::forward_as_tuple(lhs.set(), lhs.order()) ==
         std::forward_as_tuple(rhs.set(), rhs.order());
}
template <typename Set>
bool operator>(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return rhs < lhs;
}
template <typename Set>
bool operator<=(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return !(rhs < lhs);
}
template <typename Set>
bool operator>=(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return !(lhs < rhs);
}
template <typename Set>
bool operator!=(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return !(lhs == rhs);
}


class ListID {
 private:
  Vector<ListIDElement> elems_;
  Map<Ref<const PartitionSet>, int> sum_orders_;

 public:
  ListID(const Vector<ListIDElement>& elems,
         const Map<Ref<const PartitionSet>, int>& sum_orders)
      : elems_(elems), sum_orders_(sum_orders) {}

  typedef Vector<ListIDElement>::const_iterator const_iterator;

  const Vector<ListIDElement>& elems() const { return elems_; }
  const Map<Ref<const PartitionSet>, int>& sum_orders() const {
    return sum_orders_;
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
  Vector<std::tuple<const ListIDElement&, int, int>> recursion_info_;

  void init_recursion_info_aux(
      const Vector<ListID>& ids, const arma::mat& choose,
      Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info) const;
  Vector<std::tuple<const ListIDElement&, int, int>> init_recursion_info(
      const Vector<ListID>& ids, const arma::mat& choose) const;

 public:
  EdgeList(const ListID& id, const Vector<ListID>& ids, const arma::mat& choose)
      : id_(id), recursion_info_(init_recursion_info(ids, choose)) {}

  const ListID& id() const { return id_; }
  const Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info()
      const {
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
      const Vector<ListID>& ids, const arma::mat& choose,
      ListID::const_iterator curr_it,
      const Vector<ListIDElement>& left_elist_id_elems,
      Vector<std::tuple<int, int, int>>& recursion_info) const;
  Vector<std::tuple<int, int, int>> init_recursion_info(
      const Vector<ListID>& ids, const arma::mat& choose) const;

 public:
  NodeList(const ListID& id, const Vector<ListID>& ids, const arma::mat& choose)
      : id_(id), recursion_info_(init_recursion_info(ids, choose)) {}

  const ListID& id() const { return id_; }
  const Vector<std::tuple<int, int, int>>& recursion_info() const {
    return recursion_info_;
  }
};


// helper functions

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


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// EDGE SETS //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Vector<EdgeSet> create_edge_sets(VectorVector<int>& esets_inp) {
  Vector<EdgeSet> esets;
  esets.reserve(esets_inp.size());

  for (std::size_t i = 0; i < esets_inp.size(); ++i) {
    esets.emplace_back(i, std::move(esets_inp[i]));
  }

  return esets;
}

void print_set_label(const EdgeSet& eset) {
  Rcpp::Rcout << "(" << eset.label() << ")";
}

template <typename Set>
void print_set_elems(const Set& set) {
  for (std::size_t i = 0; i < set.elems().size(); ++i) {
    Rcpp::Rcout << set.elems()[i];
    if (i < set.elems().size() - 1) Rcpp::Rcout << ",";
  }
}

// [[Rcpp::export]]
int print_esets(VectorVector<int> esets_inp) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);

  for (const auto& eset : esets) {
    print_set_label(eset);
    Rcpp::Rcout << " - ";
    print_set_elems(eset);
    Rcpp::Rcout << std::endl;
  }

  return esets.size();
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// PARTITION SETS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void partition_edge_sets_aux(const Vector<EdgeSet>& esets, std::size_t curr_ind,
                             const Vector<int>& curr_label,
                             const Vector<int>& curr_elems,
                             Vector<PartitionSet>& psets) {
  // 1) If the current partitioned set is empty, we can exit the function.
  // Otherwise, the current partitioned set is non-empty.
  // 2) Furthermore, if we have traversed over all the edge sets, cache the
  // current partitioned edge set and exit the function.
  if (curr_elems.empty()) {
    return;
  } else if (curr_ind == esets.size()) {
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
  ++curr_ind;
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

void print_set_label(const PartitionSet& pset) {
  Rcpp::Rcout << "(";
  for (std::size_t i = 0; i < pset.label().size(); ++i) {
    Rcpp::Rcout << pset.label()[i];
    if (i < pset.label().size() - 1) Rcpp::Rcout << ",";
  }
  Rcpp::Rcout << ")";
}

// [[Rcpp::export]]
int print_psets(VectorVector<int> esets_inp) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);

  for (const auto& pset : psets) {
    print_set_label(pset);
    Rcpp::Rcout << " - ";
    print_set_elems(pset);
    Rcpp::Rcout << std::endl;
  }

  return psets.size();
}

// [[Rcpp::export]]
int print_edge_pset_map(VectorVector<int> esets_inp) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
  Map<int, const PartitionSet&> edge_psets = create_edge_pset_map(psets);

  for (const auto& edge_pset : edge_psets) {
    Rcpp::Rcout << edge_pset.first;
    Rcpp::Rcout << " - ";
    print_set_label(edge_pset.second);
    Rcpp::Rcout << std::endl;
  }

  return edge_psets.size();
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// LIST IDS //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void get_list_ids_aux(const Vector<PartitionSet>& psets,
                      std::size_t curr_ps_ind, int max_order,
                      const Vector<ListIDElement>& curr_elems,
                      const Map<Ref<const PartitionSet>, int>& curr_sum_orders,
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
  for (int order = 1; order <= max_order; ++order) {
    for (std::size_t ps_ind = 0; ps_ind < psets.size(); ++ps_ind) {
      get_list_ids_aux(psets, ps_ind, max_order - order,
                       {ListIDElement(psets[ps_ind], order)},
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

template <typename Set>
void print_id_element(const IDElement<Set>& id_elem) {
  Rcpp::Rcout << "[";
  print_set_label(id_elem.set());
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
    print_set_label(it->first);
    Rcpp::Rcout << ":" << it->second;
    if (std::distance(it, id.sum_orders().end()) > 1) Rcpp::Rcout << " , ";
  }
  Rcpp::Rcout << "}";
}

// [[Rcpp::export]]
int print_list_ids(VectorVector<int> esets_inp, int max_order) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
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
    const Vector<ListID>& ids, const arma::mat& choose,
    Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info) const {
  // Loop through the unique ID elements in the edge list ID and cache the
  // corresponding recursion information 3-tuples.
  for (auto curr_it = id_.elems().begin(); curr_it != id_.elems().end();) {
    // Determine the node list ID elements associated with the child node.
    Vector<ListIDElement> nlist_id_elems;
    nlist_id_elems.reserve(id_.elems().size() - 1);

    for (auto it = id_.elems().begin(); it != id_.elems().end(); ++it) {
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
    curr_it = find_next_id_elem(curr_it, id_.elems().end(), curr_id_elem);
  }
}

Vector<std::tuple<const ListIDElement&, int, int>>
EdgeList::init_recursion_info(const Vector<ListID>& ids,
                              const arma::mat& choose) const {
  Vector<std::tuple<const ListIDElement&, int, int>> recursion_info;
  recursion_info.reserve(id_.elems().size());
  init_recursion_info_aux(ids, choose, recursion_info);

  return recursion_info;
}

// [[Rcpp::export]]
int print_elist_recursion_info(VectorVector<int> esets_inp, int max_order,
                               int elist_ind) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
  Vector<ListID> ids = get_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  EdgeList elist(ids[elist_ind], ids, choose);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(elist.id());
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& tup : elist.recursion_info()) {
    const ListIDElement& curr_id_elem = std::get<0>(tup);
    int nlist_ind = std::get<1>(tup);
    int choose_coef = std::get<2>(tup);

    print_id_element(curr_id_elem);
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
    const Vector<ListID>& ids, const arma::mat& choose,
    ListID::const_iterator curr_it,
    const Vector<ListIDElement>& left_elist_id_elems,
    Vector<std::tuple<int, int, int>>& recursion_info) const {
  // Determine the edge list ID elements associated with the right child edge.
  Vector<ListIDElement> right_elist_id_elems;
  right_elist_id_elems.reserve(id_.elems().size() - left_elist_id_elems.size());
  std::set_difference(id_.elems().begin(), id_.elems().end(),
                      left_elist_id_elems.begin(), left_elist_id_elems.end(),
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
  for (auto it = curr_it; it != id_.elems().end();) {
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
    it = find_next_id_elem(it, id_.elems().end(), curr_id_elem);
  }
}

Vector<std::tuple<int, int, int>> NodeList::init_recursion_info(
    const Vector<ListID>& ids, const arma::mat& choose) const {
  Vector<std::tuple<int, int, int>> recursion_info;
  recursion_info.reserve(std::pow(2, id_.elems().size()));
  init_recursion_info_aux(ids, choose, id_.elems().begin(), {}, recursion_info);

  return recursion_info;
}

// [[Rcpp::export]]
int print_nlist_recursion_info(VectorVector<int> esets_inp, int max_order,
                               int nlist_ind) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
  Vector<ListID> ids = get_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  NodeList nlist(ids[nlist_ind], ids, choose);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(nlist.id());
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

