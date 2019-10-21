#ifndef __random_geometric_kleinberg_H__
#define __random_geometric_kleinberg_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

typedef tuple < pair<size_t,size_t>, double, double > edge_distance;

bool compare_distance(
                    const edge_distance &p1,
                    const edge_distance &p2
                 );

pair < size_t, vector < pair < size_t, size_t > > > random_geometric_kleinberg_edge_list(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        bool use_theory_algorithm = false,
        size_t seed = 0,
        double epsilon = 0
        );

vector < set < size_t > * > theoretical_random_geometric_kleinberg_neighbor_set(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        size_t seed = 0,
        double epsilon = 0
        );

vector < set < size_t > * > random_geometric_kleinberg_neighbor_set(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        size_t seed = 0
        );

tuple < size_t, vector <size_t>, vector<size_t> > random_geometric_kleinberg_coord_lists(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        bool use_theory_algorithm = false,
        size_t seed = 0,
        double epsilon = 0
        );

double f1(const size_t &N, const double &k, const double &kappa, const double &epsilon, double R);

double f1P(const size_t &N, const double &k, const double &kappa, const double &epsilon, double R);

double f2(const size_t &N, const double &k, const double &epsilon, double R);

double f2P(const size_t &N, const double &k, const double &epsilon, double R);

double find_R_min(const size_t &N, const double &k, const double &kappa, double &epsilon);

#endif
