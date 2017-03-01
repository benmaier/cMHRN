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
#include "kleinberg.h"

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

using namespace std;

tuple < size_t, vector <size_t>, vector<size_t> > kleinberg_coord_lists(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        size_t seed
        )
{
    vector < set < size_t > * > G = kleinberg_neighbor_set(N,k,mu,use_giant_component,seed);
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

pair < size_t, vector < pair < size_t, size_t > > > kleinberg_edge_list(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        size_t seed
        )
{
    size_t new_N = N;
    vector < set < size_t > * > G = kleinberg_neighbor_set(N,k,mu,use_giant_component,seed);
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

vector < set < size_t > * > kleinberg_neighbor_set(
        size_t N,
        double k,
        double mu,
        bool use_giant_component,
        size_t seed
        )
{

    assert(k>0);
    assert(N>1);
    assert(mu<=1.);

    //initialize random generators
    default_random_engine generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);
    //binomial_distribution<int> distribution(9,0.5)
    
    vector < set < size_t > * > G;
    for(size_t node=0; node<N; node++)
    {
        G.push_back( new set <size_t> );
    }

    //get probability mass function for lattice distance
    vector < double > pmf = get_kleinberg_pmf(N,k,mu);


    double k_meas = 0.;
    for(int i=0; i<N-1; i++)
        k_meas += pmf[i];

    //loop over pairs
    for (size_t u=0; u<N-1; u++)
    {
        for (size_t v=u+1; v<N; v++)
        {
            //assign edge according to probability given lattice distance of pair
            if (uni_distribution(generator) < pmf[v-u-1])
            {
                G[u]->insert(v);
                G[v]->insert(u);
            }
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

