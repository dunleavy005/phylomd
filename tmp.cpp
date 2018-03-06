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
    for (auto& elem : elems_) elem -= 1;
    std::sort(elems_.begin(), elems_.end());
  }

  int label() const { return label_; }
  const Vector<int>& elems() const { return elems_; }
};

bool operator<(const EdgeSet& lhs, const EdgeSet& rhs) {
  return lhs.label() < rhs.label();
}
bool operator==(const EdgeSet& lhs, const EdgeSet& rhs) {
  return lhs.label() == rhs.label();
}
bool operator>(const EdgeSet& lhs, const EdgeSet& rhs) { return rhs < lhs; }
bool operator<=(const EdgeSet& lhs, const EdgeSet& rhs) { return !(rhs < lhs); }
bool operator>=(const EdgeSet& lhs, const EdgeSet& rhs) { return !(lhs < rhs); }
bool operator!=(const EdgeSet& lhs, const EdgeSet& rhs) {
  return !(lhs == rhs);
}


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

template <typename Set>
class FlatIDElement {
 private:
  Ref<const Set> set_;

 public:
  FlatIDElement(const Set& set) : set_(set) {}

  const Set& set() const { return set_; }
};

typedef IDElement<EdgeSet> MomentDerivativeIDElement;
typedef FlatIDElement<EdgeSet> FlatMomentDerivativeIDElement;
typedef IDElement<PartitionSet> ListIDElement;
typedef FlatIDElement<PartitionSet> FlatListIDElement;

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

template <typename Set>
bool operator<(const FlatIDElement<Set>& lhs, const FlatIDElement<Set>& rhs) {
  return lhs.set() < rhs.set();
}
template <typename Set>
bool operator==(const FlatIDElement<Set>& lhs, const FlatIDElement<Set>& rhs) {
  return lhs.set() == rhs.set();
}
template <typename Set>
bool operator>(const FlatIDElement<Set>& lhs, const FlatIDElement<Set>& rhs) {
  return rhs < lhs;
}
template <typename Set>
bool operator<=(const FlatIDElement<Set>& lhs, const FlatIDElement<Set>& rhs) {
  return !(rhs < lhs);
}
template <typename Set>
bool operator>=(const FlatIDElement<Set>& lhs, const FlatIDElement<Set>& rhs) {
  return !(lhs < rhs);
}
template <typename Set>
bool operator!=(const FlatIDElement<Set>& lhs, const FlatIDElement<Set>& rhs) {
  return !(lhs == rhs);
}


typedef Vector<MomentDerivativeIDElement> MomentDerivativeID;

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

typedef Vector<FlatMomentDerivativeIDElement> FlatMomentDerivativeID;
typedef Vector<FlatListIDElement> FlatListID;


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
  Map<Ref<const PartitionSet>, std::pair<int, int>> recursion_info_inds_;
  arma::mat elems_;

  void init_recursion_info_aux(
      const Vector<ListID>& ids, const arma::mat& choose,
      Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info) const;
  Vector<std::tuple<const ListIDElement&, int, int>> init_recursion_info(
      const Vector<ListID>& ids, const arma::mat& choose) const;
  Map<Ref<const PartitionSet>, std::pair<int, int>> init_recursion_info_inds()
      const;

 public:
  EdgeList(const ListID& id, const Vector<ListID>& ids, const arma::mat& choose,
           int nrow, int ncol)
      : id_(id),
        recursion_info_(init_recursion_info(ids, choose)),
        recursion_info_inds_(init_recursion_info_inds()),
        elems_(arma::mat(nrow, ncol, arma::fill::zeros)) {}

  arma::mat& elems() { return elems_; }
  const arma::mat& elems() const { return elems_; }

  const ListID& id() const { return id_; }
  const Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info()
      const {
    return recursion_info_;
  }
  const Map<Ref<const PartitionSet>, std::pair<int, int>>& recursion_info_inds()
      const {
    return recursion_info_inds_;
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
  arma::mat elems_;

  void init_recursion_info_aux(
      const Vector<ListID>& ids, const arma::mat& choose,
      ListID::const_iterator curr_it,
      const Vector<ListIDElement>& left_elist_id_elems,
      Vector<std::tuple<int, int, int>>& recursion_info) const;
  Vector<std::tuple<int, int, int>> init_recursion_info(
      const Vector<ListID>& ids, const arma::mat& choose) const;

 public:
  NodeList(const ListID& id, const Vector<ListID>& ids, const arma::mat& choose,
           int nrow, int ncol)
      : id_(id),
        recursion_info_(init_recursion_info(ids, choose)),
        elems_(arma::mat(nrow, ncol, arma::fill::zeros)) {}

  arma::mat& elems() { return elems_; }
  const arma::mat& elems() const { return elems_; }

  const ListID& id() const { return id_; }
  const Vector<std::tuple<int, int, int>>& recursion_info() const {
    return recursion_info_;
  }
};


// helper functions

int find_moment_derivative_id_index(const Vector<MomentDerivativeID>& ids,
                                    const MomentDerivativeID& id) {
  return std::lower_bound(ids.begin(), ids.end(), id) - ids.begin();
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
/////////////////////////// MOMENT/DERIVATIVE IDS //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void get_moment_derivative_ids_aux(const Vector<EdgeSet>& esets,
                                   std::size_t curr_es_ind, int max_order,
                                   const MomentDerivativeID& curr_id,
                                   Vector<MomentDerivativeID>& ids) {
  // Cache the current moment/derivative ID.
  ids.push_back(curr_id);

  // Find the next moment/derivative IDs.
  // (Note: the loop body will not be entered if there are no more
  // moment/derivative IDs to find.)
  for (int order = 1; order <= max_order; ++order) {
    for (std::size_t es_ind = curr_es_ind; es_ind < esets.size(); ++es_ind) {
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

  // Sort the moment/derivative IDs.
  std::sort(ids.begin(), ids.end());

  return ids;
}

template <typename Set>
void print_id_element(const IDElement<Set>& id_elem) {
  Rcpp::Rcout << "[";
  print_set_label(id_elem.set());
  Rcpp::Rcout << "-" << id_elem.order() << "]";
}

void print_moment_derivative_id(const MomentDerivativeID& id) {
  Rcpp::Rcout << "< ";
  for (std::size_t i = 0; i < id.size(); ++i) {
    print_id_element(id[i]);
    if (i < id.size() - 1) Rcpp::Rcout << " , ";
  }
  Rcpp::Rcout << " >";
}

// [[Rcpp::export]]
int print_moment_derivative_ids(VectorVector<int> esets_inp, int max_order) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<MomentDerivativeID> ids = get_moment_derivative_ids(esets, max_order);

  for (const auto& id : ids) {
    print_moment_derivative_id(id);
    Rcpp::Rcout << std::endl;
  }

  return ids.size();
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

  // Sort the list IDs.
  std::sort(ids.begin(), ids.end(),
            [](const ListID& left_id, const ListID& right_id) {
              return left_id.elems() < right_id.elems();
            });

  return ids;
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

// [[Rcpp::export]]
int print_elist_recursion_info(VectorVector<int> esets_inp, int max_order,
                               int elist_ind) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
  Vector<ListID> ids = get_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  EdgeList elist(ids[elist_ind], ids, choose, 0, 0);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(elist.id());
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& ri_tup : elist.recursion_info()) {
    const ListIDElement& curr_id_elem = std::get<0>(ri_tup);
    int nlist_ind = std::get<1>(ri_tup);
    int choose_coef = std::get<2>(ri_tup);

    print_id_element(curr_id_elem);
    Rcpp::Rcout << " ----- ";
    print_list_id_elems(ids[nlist_ind]);
    Rcpp::Rcout << " ----- ";
    Rcpp::Rcout << choose_coef << std::endl;
  }

  return elist.recursion_info().size();
}

// [[Rcpp::export]]
int print_elist_recursion_info_inds(VectorVector<int> esets_inp, int max_order,
                                    int elist_ind) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
  Vector<ListID> ids = get_list_ids(psets, max_order);
  arma::mat choose = choose_table(max_order);
  EdgeList elist(ids[elist_ind], ids, choose, 0, 0);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(elist.id());
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& pset_inds : elist.recursion_info_inds()) {
    const PartitionSet& pset = pset_inds.first;
    int begin_ind, end_ind;
    std::tie(begin_ind, end_ind) = pset_inds.second;

    print_set_label(pset);
    Rcpp::Rcout << " ----- ";
    Rcpp::Rcout << begin_ind << "," << end_ind << std::endl;
  }

  return elist.recursion_info_inds().size();
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
  NodeList nlist(ids[nlist_ind], ids, choose, 0, 0);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(nlist.id());
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& ri_tup : nlist.recursion_info()) {
    int left_elist_ind = std::get<0>(ri_tup);
    int right_elist_ind = std::get<1>(ri_tup);
    int choose_coef = std::get<2>(ri_tup);

    print_list_id_elems(ids[left_elist_ind]);
    Rcpp::Rcout << " ----- ";
    print_list_id_elems(ids[right_elist_ind]);
    Rcpp::Rcout << " ----- ";
    Rcpp::Rcout << choose_coef << std::endl;
  }

  return nlist.recursion_info().size();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// MISCELLANEOUS ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void find_connected_moment_derivative_ids_aux(
    const ListID& list_id, const Vector<EdgeSet>& esets,
    ListID::const_iterator curr_it, int curr_order,
    const Map<int, int>& eset_labels_orders,
    Vector<MomentDerivativeID>& md_ids) {
  // If we have traversed over all the list ID elements, then cache the current
  // moment/derivative ID and exit the function.
  if (curr_it == list_id.elems().end()) {
    MomentDerivativeID curr_md_id;
    curr_md_id.reserve(eset_labels_orders.size());

    for (const auto& eset_label_order : eset_labels_orders) {
      curr_md_id.emplace_back(esets[eset_label_order.first],
                              eset_label_order.second);
    }

    md_ids.push_back(std::move(curr_md_id));
    return;
  }

  // Otherwise, loop over the possible edge set labels in the current list ID
  // element.
  for (auto eset_label : curr_it->set().label()) {
    // Construct the next (`eset_label`, `order`) map by inserting the current
    // (`eset_label`, 1) pair into the current map.
    Map<int, int> next_eset_labels_orders = eset_labels_orders;
    auto insert_results = next_eset_labels_orders.emplace(eset_label, 1);
    if (!insert_results.second) insert_results.first->second += 1;

    // Is the current order equal to the order of the current list ID element?
    if (curr_order == curr_it->order()) {
      // If so, recurse over the next list ID element.
      find_connected_moment_derivative_ids_aux(list_id, esets, curr_it + 1, 1,
                                               next_eset_labels_orders, md_ids);
    } else {
      // Otherwise, recurse over the next order.
      find_connected_moment_derivative_ids_aux(list_id, esets, curr_it,
                                               curr_order + 1,
                                               next_eset_labels_orders, md_ids);
    }
  }
}

Vector<MomentDerivativeID> find_connected_moment_derivative_ids(
    const ListID& list_id, const Vector<EdgeSet>& esets) {
  Vector<MomentDerivativeID> md_ids;
  find_connected_moment_derivative_ids_aux(
      list_id, esets, list_id.elems().begin(), 1, {}, md_ids);

  // Keep only the unique moment/derivative IDs.
  std::sort(md_ids.begin(), md_ids.end());
  auto end_it = std::unique(md_ids.begin(), md_ids.end());
  md_ids.resize(end_it - md_ids.begin());

  return md_ids;
}

// [[Rcpp::export]]
int print_connected_moment_derivative_ids(VectorVector<int> esets_inp,
                                          int max_order, int list_ind) {
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);
  Vector<ListID> list_ids = get_list_ids(psets, max_order);
  Vector<MomentDerivativeID> md_ids =
      find_connected_moment_derivative_ids(list_ids[list_ind], esets);

  Rcpp::Rcout << "List ID Elements: ";
  print_list_id_elems(list_ids[list_ind]);
  Rcpp::Rcout << "\n" << std::endl;

  for (const auto& md_id : md_ids) {
    print_moment_derivative_id(md_id);
    Rcpp::Rcout << std::endl;
  }

  return md_ids.size();
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////// PHYLO MOMENTS/DERIVATIVES ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//if (!tree.inherits("phylo")) Rcpp::stop("The tree object must be of class 'phylo'.");
//if (tree.attr("order") != "cladewise") Rcpp::stop("The edge matrix must be in 'cladewise' order.");
//if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

void update_node_lists(int child_node_ind, const arma::uvec& child_edge_inds,
                       const Vector<EdgeList>& elists,
                       Vector<NodeList>& nlists) {
  for (std::size_t nlist_ind = 0; nlist_ind < nlists.size(); ++nlist_ind) {
    for (const auto& ri_tup : nlists[nlist_ind].recursion_info()) {
      // Extract the recursion information 3-tuple elements.
      int left_elist_ind = std::get<0>(ri_tup);
      int right_elist_ind = std::get<1>(ri_tup);
      int choose_coef = std::get<2>(ri_tup);

      // Update the current node list.
      nlists[nlist_ind].elems().col(child_node_ind) +=
          elists[left_elist_ind].elems().col(child_edge_inds(0)) %
          elists[right_elist_ind].elems().col(child_edge_inds(1)) * choose_coef;
    }
  }
}

void update_edge_lists(
    int edge_ind, int child_node_ind, const arma::cube& edge_mds,
    Map<int, const PartitionSet&>::const_iterator edge_pset_it,
    Map<int, const PartitionSet&>::const_iterator edge_psets_end_it,
    const Vector<NodeList>& nlists, Vector<EdgeList>& elists) {
  for (std::size_t elist_ind = 0; elist_ind < elists.size(); ++elist_ind) {
    // Is the current edge in any of the partition sets?
    if (edge_pset_it != edge_psets_end_it) {
      // If so, is this partition set in the current edge list ID?
      auto pset_inds_it =
          elists[elist_ind].recursion_info_inds().find(edge_pset_it->second);
      if (pset_inds_it != elists[elist_ind].recursion_info_inds().end()) {
        // If so, we must account for the possible assignments of ID elements to
        // the current edge.
        int begin_ind, end_ind;
        std::tie(begin_ind, end_ind) = pset_inds_it->second;

        for (int ri_ind = begin_ind; ri_ind < end_ind; ++ri_ind) {
          // Extract the recursion information 3-tuple elements.
          const auto& ri_tup = elists[elist_ind].recursion_info()[ri_ind];
          const ListIDElement& curr_id_elem = std::get<0>(ri_tup);
          int nlist_ind = std::get<1>(ri_tup);
          int choose_coef = std::get<2>(ri_tup);

          // Update the current edge list.
          elists[elist_ind].elems().col(edge_ind) +=
              edge_mds.slice(curr_id_elem.order()) *
              nlists[nlist_ind].elems().col(child_node_ind) * choose_coef;
        }
      }
    }

    // Update the current edge list.
    // (Note: this update accounts for the case where no ID element is assigned
    // to the current edge.)
    elists[elist_ind].elems().col(edge_ind) +=
        edge_mds.slice(0) * nlists[elist_ind].elems().col(child_node_ind);
  }
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

Map<std::string, double> phylo_moments_derivatives(
    const Rcpp::List& tree, const arma::mat& Q, const arma::mat& B,
    const arma::vec& pi, VectorVector<int>& esets_inp, int max_order, Mode mode,
    const arma::ivec& tip_data) {
  // Extract the necessary elements from the tree object.
  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  const arma::vec& edge_lengths = tree["edge.length"];
  int num_int_nodes = tree["Nnode"];

  edge -= 1;
  int num_edges = edge.n_rows;
  int num_term_nodes = tip_labels.size();
  int root_node_ind = num_term_nodes;

  // Initialize the counting table.
  arma::mat choose = choose_table(max_order);

  // Construct the edge sets and partition sets.
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);

  // Form the (edge index, partition set) inverse map.
  Map<int, const PartitionSet&> edge_psets = create_edge_pset_map(psets);

  // Generate the moment/derivative IDs and list IDs.
  Vector<MomentDerivativeID> md_ids =
      get_moment_derivative_ids(esets, max_order);
  Vector<ListID> list_ids = get_list_ids(psets, max_order);

  // Create the node lists and edge lists.
  Vector<NodeList> nlists;
  nlists.reserve(list_ids.size());
  Vector<EdgeList> elists;
  elists.reserve(list_ids.size());

  for (const auto& list_id : list_ids) {
    nlists.emplace_back(list_id, list_ids, choose, Q.n_rows,
                        num_int_nodes + num_term_nodes);
    elists.emplace_back(list_id, list_ids, choose, Q.n_rows, num_edges);
  }

  // Perform the postorder recursion.
  // (Note: we iterate up the edge matrix because, by default, it is arranged in
  // a preorder format.)
  for (int edge_ind = num_edges - 1; edge_ind >= 0; --edge_ind) {
    int child_node_ind = edge(edge_ind, 1);

    ////
    //// First, we consider the node lists at the current child node.
    ////

    // Is the current child node an internal node?
    if (child_node_ind >= root_node_ind) {
      // If so, update each node list using the appropriate edge lists at the
      // two branches below the current child node.
      arma::uvec child_edge_inds = arma::find(edge.col(0) == child_node_ind);
      update_node_lists(child_node_ind, child_edge_inds, elists, nlists);
    } else {
      // Otherwise, initialize the "empty" (i.e. partial likelihood) node list
      // at the given terminal node.
      // (Note: the "empty" node list is always the first element in the node
      // list vector.)

      // Do we have observed data at the terminal nodes?
      if (tip_data.n_elem > 0) {
        nlists[0].elems()(tip_data(child_node_ind), child_node_ind) = 1;
      } else {
        nlists[0].elems().col(child_node_ind).fill(1);
      }
    }

    ////
    //// Second, we consider the edge lists at the current edge.
    ////

    // Is the current edge in any of the partition sets?
    auto edge_pset_it = edge_psets.find(edge_ind);

    // If so, compute the necessary CTMC moments/derivatives at the current
    // edge.
    // Otherwise, compute only the zeroth CTMC moment/derivative (i.e. the
    // transition probability matrix) at the current edge.
    int edge_max_order = (edge_pset_it != edge_psets.end()) ? max_order : 0;
    arma::cube edge_mds = ctmc_moments_derivatives(edge_lengths(edge_ind), Q, B,
                                                   edge_max_order, mode);

    // Update the edge lists at the current edge.
    update_edge_lists(edge_ind, child_node_ind, edge_mds, edge_pset_it,
                      edge_psets.end(), nlists, elists);
  }

  // Update the node lists at the root node.
  arma::uvec root_edge_inds = arma::find(edge.col(0) == root_node_ind);
  update_node_lists(root_node_ind, root_edge_inds, elists, nlists);

  // Connect the list IDs to the moment/derivative IDs and calculate the
  // moments/derivatives by accumulating dot products between the root
  // distribution and the node lists at the root node.
  arma::vec phylo_mds(md_ids.size(), arma::fill::zeros);
  for (const auto& nlist : nlists) {
    double md_part = arma::dot(pi, nlist.elems().col(root_node_ind));
    Vector<MomentDerivativeID> conn_md_ids =
        find_connected_moment_derivative_ids(nlist.id(), esets);

    for (const auto& conn_md_id : conn_md_ids) {
      int conn_md_id_ind = find_moment_derivative_id_index(md_ids, conn_md_id);
      phylo_mds(conn_md_id_ind) += md_part;
    }
  }

  // Prepare the output map.
  Map<std::string, double> phylo_mds_map;
  for (std::size_t i = 0; i < md_ids.size(); ++i) {
    // Create the moment/derivative ID label.
    std::string md_id_label = create_moment_derivative_id_label(md_ids[i]);

    // Cache the moment/derivative.
    phylo_mds_map.emplace(std::move(md_id_label), phylo_mds(i));
  }

  return phylo_mds_map;
}

// [[Rcpp::export]]
Map<std::string, double> phylomd_wrap(const Rcpp::List& tree,
                                      const arma::mat& Q, const arma::mat& B,
                                      const arma::vec& pi,
                                      VectorVector<int> esets_inp,
                                      int max_order,
                                      const arma::ivec& tip_data) {
  return phylo_moments_derivatives(tree, Q, B, pi, esets_inp, max_order,
                                   Mode::MOMENTS, tip_data);
}
