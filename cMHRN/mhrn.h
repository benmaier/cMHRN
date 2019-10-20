#ifndef __MHRN_H__
#define __MHRN_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

double f(const size_t &B, const size_t &L, const double &k, double xi);

double fP(const size_t &B, const size_t &L, const double &k, double xi);

double find_xi_min(const size_t &B, const size_t &L, const double &k);

tuple < size_t, vector <size_t>, vector<size_t> > fast_mhrn_coord_lists(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component = false,
        bool delete_non_giant_component_nodes = true,
        bool allow_probability_redistribution = false,
        size_t seed = 0
        );

pair < size_t, vector < pair < size_t, size_t > > > fast_mhrn_edge_list(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component = false,
        bool delete_non_giant_component_nodes = true,
        bool allow_probability_redistribution = false,
        size_t seed = 0
        );

vector < set < size_t > * > fast_mhrn_neighbor_set(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component = false,
        bool allow_probability_redistribution = false,
        size_t seed = 0
        );

vector < pair < size_t, size_t > > fast_gnp(
        size_t N_,
        double p,
        size_t start_node = 0,
        size_t seed = 0
        );
#endif
