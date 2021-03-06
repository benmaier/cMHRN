/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */

#include "Utilities.h"
#include "mhrn.h"

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <random>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <tuple>
#include <assert.h>
#include <sstream>

using namespace std;

double f(const size_t &B, const size_t &L, const double &k, double xi)
{
    return (B-1.0) * (1.0-pow(xi,L)) - k * (1.0-xi);
}

double fP(const size_t &B, const size_t &L, const double &k, double xi)
{
    return - (B-1.0) * (double(L) * pow(xi,L-1.0)) + k;
}

double find_xi_min(const size_t &B, const size_t &L, const double &k)
{
    double x1 = 0.1;
    double x, fx, fx1;

    double eps = 1e-14;
    do
    {
        x = x1;
        fx = f(B,L,k,x);
        fx1 = fP(B,L,k,x);
        x1 = x - (fx/fx1);

    } while(fabs(x1-x) >= eps);

    return x1;
}

tuple < size_t, vector <size_t>, vector<size_t> > fast_mhrn_coord_lists(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        bool allow_probability_redistribution,
        size_t seed
        )
{
    vector < set < size_t > * > G = fast_mhrn_neighbor_set(B,L,k,xi,
                                                           use_giant_component,
                                                           allow_probability_redistribution,
                                                           seed);
    size_t N = pow(B,L);
    size_t new_N = N;
    vector < size_t > rows;
    vector < size_t > cols;

    if ( use_giant_component && delete_non_giant_component_nodes )
    {
        vector < size_t > map_to_new_ids(N);
        size_t current_id = 0;
        for(size_t u = 0; u < N; u++)
            if (G[u]->size()>0)
            {
                map_to_new_ids[u] = current_id;
                current_id++;
            }

        new_N = current_id;

        for(size_t u = 0; u < N; u++)
        {
            for( auto const& v: *G[u] )
            {
                rows.push_back(map_to_new_ids[u]);
                cols.push_back(map_to_new_ids[v]);
            }
            delete G[u];
        }
    }
    else
    {
        for(size_t u = 0; u < N; u++)
        {
            for( auto const& v: *G[u] )
            {
                rows.push_back(u);
                cols.push_back(v);
            }
            delete G[u];
        }
    }
    
    return make_tuple(new_N,rows,cols);
}

pair < size_t, vector < pair < size_t, size_t > > > fast_mhrn_edge_list(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        bool allow_probability_redistribution,
        size_t seed
        )
{
    size_t N = pow(B,L);
    size_t new_N = N;
    vector < set < size_t > * > G = fast_mhrn_neighbor_set(B,L,k,xi,
                                                           use_giant_component,
                                                           allow_probability_redistribution,
                                                           seed);
    vector < pair < size_t, size_t > > edge_list;

    if ( use_giant_component && delete_non_giant_component_nodes )
    {
        vector < size_t > map_to_new_ids(N);
        size_t current_id = 0;
        for(size_t u = 0; u < N; u++)
            if (G[u]->size()>0)
            {
                map_to_new_ids[u] = current_id;
                current_id++;
            }

        new_N = current_id;

        for(size_t u = 0; u < N; u++)
        {
            size_t u_ = map_to_new_ids[u];
            for( auto const& v: *G[u] )
            {
                size_t v_ = map_to_new_ids[v];
                if (u_<v_)
                {
                    edge_list.push_back( make_pair( u_, v_ ) );
                }
            }
            delete G[u];
        }
    }
    else
    {
        for(size_t u = 0; u < N; u++)
        {
            for( auto const& v: *G[u] )
            {
                if (u<v)
                {
                    edge_list.push_back( make_pair(u,v) );
                }
            }
            delete G[u];
        }
    }
    
    return make_pair(new_N,edge_list);
}

vector < set < size_t > * > fast_mhrn_neighbor_set(
        size_t B,
        size_t L,
        double k,
        double xi,
        bool use_giant_component,
        bool allow_probability_redistribution,
        size_t seed
        )
{

    assert(L>1);
    assert(k>0);
    assert(k<=pow(B,L)-1);
    assert(B>1);
    assert((xi>=0.0) && (xi<=B));

    //initialize random generators
    default_random_engine generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);
    //binomial_distribution<int> distribution(9,0.5)
    
    size_t N = pow(B,L);

    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
    {
        G.push_back( new set <size_t> );
    }

    double p1;
    if ( xi == 1.0 ) 
        p1 = k / double(B-1) / double(L);
    else
        p1 = k / double(B-1) * (1.0-xi) / (1.0-pow(xi,L));

    vector < double > p(L);
    p[0] = p1;
    for (size_t l=2; l<=L; l++)
        p[l-1] = p1 * pow(xi/double(B),l-1);

    if ((p1>1.0) and (not allow_probability_redistribution))
    {
        stringstream ss; 
        ss.precision(4);
        ss << "The lowest layer connection probability is >1.0. If you want to keep xi = ";
        ss << xi;

        double new_k, new_xi;

        if (xi == 1.0)
            new_k = (B-1)*L;
        else
            new_k = (B-1) * (1-pow(xi,L)) / (1-xi);

        new_xi = find_xi_min(B,L,k);

        ss << " change the mean degree to k <= ";
        ss << new_k;

        ss << "; If you want to keep k = " << k << " change the structural control parameter to xi > " << new_xi;

        ss << " You may also want to allow for redistribution of connection probability to higher layers. To this end, call the function with parameter `allow_probability_redistribution = True`."; 
            
        throw domain_error(ss.str());
    }
    else if ((p1>1.0) and (allow_probability_redistribution))
    {

        // redistribute non-creatable edges from lower to higher layers
        
        // compute number of edges per layer
        vector < double > m(L);
        for (size_t l=1; l<=L; l++)
            m[l-1] = (0.5 * N) * pow(B,l-1.0) * (B-1.0);

        size_t l = 0;
        double excess_edges = (p[0]-1.0)*m[0];

        while (excess_edges > 0.0)
        {
            l++;

            // this one is actuall already covered in assert(k<N-1) but let's keep it here for safety.
            if (l==L)
                throw domain_error("There's not enough node pairs to redistribute all excess probability.");

            excess_edges += (p[l]-1.0)*m[l];
        }

        excess_edges -= (p[l]-1.0)*m[l];

        vector < double > new_p(L);

        for (size_t j=0; j<l; j++)
            new_p[j] = 1.0;

        // Eq. (B.3) on page 235 of my dissertation
        new_p[l] = p[l] + (1.0/m[l]) * excess_edges;

        // a safety mechanism
        if (new_p[l] > 1.0)
            throw domain_error("There's something wrong with the redistribution algorithm. Please open an issue at https://github.com/benmaier/cMHRN/issues");

        for (size_t j=l+1; j<L; j++)
            new_p[j] = p[j];

        p = new_p;
    }

    //add ER graphs in lowest layer l=1
    for (size_t start_node = 0; start_node<N; start_node += B)
    {
        add_random_subgraph(B,p[0],G,generator,uni_distribution,start_node);
    }

    for(size_t l=2;l<=L; l++)
    {
        binomial_distribution<size_t> binom( size_t(0.5*pow(B,(L+l-1))*(B-1)),p[l-1]);
        size_t current_m_l = binom(generator);

        for (size_t m = 0; m< current_m_l; m++)
        {
            size_t B_l = pow(B,l);
            size_t B_lm1 = pow(B,l-1);
            bool already_contains_edge = true;
            
            do
            {
                size_t w = size_t(uni_distribution(generator)*N);
                size_t b = w / B_l;
                size_t b_lower = w / B_lm1;
                size_t v;
                uniform_int_distribution<size_t> randint(b*B_l,(b+1)*B_l-1);

                do
                    v = randint(generator);
                while (b_lower == v / B_lm1);

                already_contains_edge = G[w]->find(v) != G[w]->end();

                if ( not already_contains_edge)
                {
                    G[w]->insert(v);
                    G[v]->insert(w);
                }

            } while ( already_contains_edge );
        }
    }

    if (use_giant_component)
    {
        get_giant_component(G);
        return G;
    }
    else
    {
        return G;
    }

}

vector < pair < size_t, size_t > > fast_gnp(
        size_t N_,
        double p,
        size_t start_node,
        size_t seed
        )
{
    size_t N = N_ + start_node;

    //initialize random generators
    default_random_engine generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
        G.push_back( new set <size_t> );

    add_random_subgraph(N,p,G,generator,uni_distribution,start_node);
    vector < pair < size_t, size_t > > edge_list;

    for(size_t u = 0; u < N; u++)
        for( auto const& v: *G[u] )
        {
            if (u<v)
            {
                edge_list.push_back( make_pair(u,v) );
            }
        }
    
    return edge_list;
}
