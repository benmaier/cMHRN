#ifndef __KLEINBERG_H__
#define __KLEINBERG_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

tuple < size_t, vector <size_t>, vector<size_t> > kleinberg_coord_lists(
        size_t N,
        double k,
        double mu,
        bool use_giant_component = false,
        bool delete_non_giant_component_nodes = true,
        size_t seed = 0
        );

pair < size_t, vector < pair < size_t, size_t > > > kleinberg_edge_list(
        size_t N,
        double k,
        double mu,
        bool use_giant_component = false,
        bool delete_non_giant_component_nodes = true,
        size_t seed = 0
        );

vector < set < size_t > * > kleinberg_neighbor_set(
        size_t N,
        double k,
        double mu,
        bool use_giant_component = false,
        size_t seed = 0
        );
#endif
