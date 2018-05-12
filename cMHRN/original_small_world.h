#ifndef __original_small_world_H__
#define __original_small_world_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

pair < size_t, vector < pair < size_t, size_t > > > original_small_world_edge_list(
        size_t N,
        size_t k,
        double p,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        size_t seed
        );

vector < set < size_t > * > original_small_world_neighbor_set(
        size_t N,
        size_t k,
        double p,
        bool use_giant_component,
        size_t seed
        );

tuple < size_t, vector <size_t>, vector<size_t> > original_small_world_coord_lists(
        size_t N,
        double k,
        double p,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        size_t seed
        );
#endif
