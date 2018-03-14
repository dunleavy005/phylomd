#ifndef PHYLOMD_MOMENTDERIVATIVEID_H_
#define PHYLOMD_MOMENTDERIVATIVEID_H_

#include <cstddef>
#include <string>

#include "EdgeSet.h"
#include "IDElement.h"
#include "phylomd_types.h"


//
// This header file defines the MomentDerivativeID and FlatMomentDerivativeID
// classes and declares the MomentDerivativeID-related functions.
//


typedef IDElement<EdgeSet> MomentDerivativeIDElement;


typedef Vector<MomentDerivativeIDElement> MomentDerivativeID;


typedef FlatIDElement<EdgeSet> FlatMomentDerivativeIDElement;


typedef Vector<FlatMomentDerivativeIDElement> FlatMomentDerivativeID;


void get_moment_derivative_ids_aux(const Vector<EdgeSet>& esets,
                                   std::size_t curr_es_ind, int max_order,
                                   const MomentDerivativeID& curr_id,
                                   Vector<MomentDerivativeID>& ids);


Vector<MomentDerivativeID> get_moment_derivative_ids(
    const Vector<EdgeSet>& esets, int max_order);


int find_moment_derivative_id_index(const Vector<MomentDerivativeID>& ids,
                                    const MomentDerivativeID& id);


std::string create_moment_derivative_id_label(const MomentDerivativeID& md_id);


#endif  // PHYLOMD_MOMENTDERIVATIVEID_H_