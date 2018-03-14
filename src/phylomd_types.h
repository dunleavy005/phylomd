#ifndef PHYLOMD_PHYLOMD_TYPES_H_
#define PHYLOMD_PHYLOMD_TYPES_H_

#include <functional>
#include <map>
#include <vector>


//
// This header file defines some useful typedefs of STL class templates.
//


template <typename T1, typename T2>
using Map = std::map<T1, T2>;


template <typename T>
using Ref = std::reference_wrapper<T>;


template <typename T>
using Vector = std::vector<T>;


template <typename T>
using VectorVector = std::vector<std::vector<T>>;


#endif  // PHYLOMD_PHYLOMD_TYPES_H_