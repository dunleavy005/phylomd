#ifndef PHYLOMD_LISTID_H_
#define PHYLOMD_LISTID_H_

#include <cstddef>

#include "IDElement.h"
#include "PartitionSet.h"
#include "phylomd_types.h"


//
// This header file defines the ListID and FlatListID classes and declares the
// ListID-related functions.
//


typedef IDElement<PartitionSet> ListIDElement;


class ListID {
 private:
  Vector<ListIDElement> elems_;
  Map<Ref<const PartitionSet>, int> sum_orders_;

 public:
  ListID(const Vector<ListIDElement>& elems,
         const Map<Ref<const PartitionSet>, int>& sum_orders)
      : elems_(elems), sum_orders_(sum_orders) {}

  typedef Vector<ListIDElement>::const_iterator const_iterator;
  const_iterator begin() const { return elems_.begin(); }
  const_iterator end() const { return elems_.end(); }
  std::size_t size() const { return elems_.size(); }

  const Vector<ListIDElement>& elems() const { return elems_; }
  const Map<Ref<const PartitionSet>, int>& sum_orders() const {
    return sum_orders_;
  }
};


typedef FlatIDElement<PartitionSet> FlatListIDElement;


typedef Vector<FlatListIDElement> FlatListID;


void get_list_ids_aux(const Vector<PartitionSet>& psets,
                      std::size_t curr_ps_ind, int max_order,
                      const Vector<ListIDElement>& curr_elems,
                      const Map<Ref<const PartitionSet>, int>& curr_sum_orders,
                      Vector<ListID>& ids);


Vector<ListID> get_list_ids(const Vector<PartitionSet>& psets, int max_order);


int find_list_id_index(const Vector<ListID>& ids,
                       const Vector<ListIDElement>& id_elems);


ListID::const_iterator find_next_id_elem(ListID::const_iterator begin_it,
                                         ListID::const_iterator end_it,
                                         const ListIDElement& curr_id_elem);


#endif  // PHYLOMD_LISTID_H_