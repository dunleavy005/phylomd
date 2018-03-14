#ifndef PHYLOMD_IDELEMENT_H_
#define PHYLOMD_IDELEMENT_H_

#include <tuple>

#include "phylomd_types.h"


//
// This header file defines the IDElement and FlatIDElement class templates.
//


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


template <typename Set>
inline bool operator<(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return std::forward_as_tuple(lhs.set(), lhs.order()) <
         std::forward_as_tuple(rhs.set(), rhs.order());
}
template <typename Set>
inline bool operator==(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return std::forward_as_tuple(lhs.set(), lhs.order()) ==
         std::forward_as_tuple(rhs.set(), rhs.order());
}
template <typename Set>
inline bool operator>(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return rhs < lhs;
}
template <typename Set>
inline bool operator<=(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return !(rhs < lhs);
}
template <typename Set>
inline bool operator>=(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return !(lhs < rhs);
}
template <typename Set>
inline bool operator!=(const IDElement<Set>& lhs, const IDElement<Set>& rhs) {
  return !(lhs == rhs);
}


template <typename Set>
inline bool operator<(const FlatIDElement<Set>& lhs,
                      const FlatIDElement<Set>& rhs) {
  return lhs.set() < rhs.set();
}
template <typename Set>
inline bool operator==(const FlatIDElement<Set>& lhs,
                       const FlatIDElement<Set>& rhs) {
  return lhs.set() == rhs.set();
}
template <typename Set>
inline bool operator>(const FlatIDElement<Set>& lhs,
                      const FlatIDElement<Set>& rhs) {
  return rhs < lhs;
}
template <typename Set>
inline bool operator<=(const FlatIDElement<Set>& lhs,
                       const FlatIDElement<Set>& rhs) {
  return !(rhs < lhs);
}
template <typename Set>
inline bool operator>=(const FlatIDElement<Set>& lhs,
                       const FlatIDElement<Set>& rhs) {
  return !(lhs < rhs);
}
template <typename Set>
inline bool operator!=(const FlatIDElement<Set>& lhs,
                       const FlatIDElement<Set>& rhs) {
  return !(lhs == rhs);
}


#endif  // PHYLOMD_IDELEMENT_H_