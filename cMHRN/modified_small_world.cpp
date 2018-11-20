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
#include "modified_small_world.h"

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

tuple < size_t, vector <size_t>, vector<size_t> > modified_small_world_coord_lists(
        size_t N,
        double k,
        double p,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        bool use_fast_algorithm,
        size_t seed
        )
{
    vector < set < size_t > * > G;

    if (use_fast_algorithm)
        G = fast_modified_small_world_neighbor_set(N,k,p,use_giant_component,seed);
    else
        G = modified_small_world_neighbor_set(N,k,p,use_giant_component,seed);

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
    
pair < size_t, vector < pair < size_t, size_t > > > modified_small_world_edge_list(
        size_t N,
        size_t k,
        double p,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        bool use_fast_algorithm,
        size_t seed
        )
{
    size_t new_N = N;

    vector < set < size_t > * > G;

    if (use_fast_algorithm)
        G = fast_modified_small_world_neighbor_set(N,k,p,use_giant_component,seed);
    else
        G = modified_small_world_neighbor_set(N,k,p,use_giant_component,seed);

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

vector < set < size_t > * > modified_small_world_neighbor_set(
        size_t N,
        size_t k,
        double p,
        bool use_giant_component,
        size_t seed
        )
{

    assert(k>0);
    assert(N>1);
    assert(p <= 1.);
    assert(p >= 0.);
    assert(k % 2 == 0);
    assert(k < N);

    double _k = (double) k;

    double p0 = _k / (_k - p*_k + p*(N-1.0));
    double p1 = p0 * p;

    size_t max_neighbor = k / 2;

    //initialize random generators
    default_random_engine generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
   
    uniform_real_distribution<double> random_number(0., 1.);
    uniform_int_distribution<size_t> random_node(0, N-1);
    
    vector < set < size_t > * > G;
    for (size_t node = 0; node < N; ++node)
        G.push_back( new set < size_t >() );

    // loop over all pairs and draw according to the right probability
    // (this is a lazy slow algorithm running in O(N^2) time
    for (size_t i = 0; i < N-1; ++i)
    {
        for (size_t j = i+1; j < N; ++j)
        {
            size_t distance = j - i;
            double probability = p1;

            if (distance <= max_neighbor or (N - distance) <= max_neighbor)
                probability = p0;

            if (random_number(generator) < probability)
            {
                G[i]->insert(j);
                G[j]->insert(i);
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

vector < set < size_t > * > fast_modified_small_world_neighbor_set(
        size_t N,
        size_t k,
        double p,
        bool use_giant_component,
        size_t seed
        )
{

    assert(k>0);
    assert(N>1);
    assert(p <= 1.);
    assert(p >= 0.);
    assert(k % 2 == 0);
    assert(k < N);

    double _k = (double) k;

    double p0 = _k / (_k - p*_k + p*(N-1.0));
    double p1 = p0 * p;

    size_t max_neighbor = k / 2;

    //initialize random generators
    default_random_engine generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);
   
    uniform_real_distribution<double> random_number(0., 1.);
    uniform_int_distribution<size_t> random_node(0, N-1);
    uniform_int_distribution<size_t> random_long_range_neighbor(1, N-1-k);
    binomial_distribution<size_t> random_long_range_edges(N*(N-1-k) / 2, p1);
    
    vector < set < size_t > * > G;

    for (size_t node = 0; node < N; ++node)
        G.push_back( new set < size_t >() );

    // loop over all short-range pairs and draw according to the right probability
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = i+1; j < i + max_neighbor +1; ++j)
        {
            size_t v = j % N;
            if (random_number(generator) < p0)
            {
                G[i]->insert(v);
                G[v]->insert(i);
            }
        }
    }

    size_t mL = random_long_range_edges(generator);

    for(size_t m = 0; m < mL; ++m)
    {
        size_t u, v;

        do
        {
            u = random_node(generator);
            v = u + max_neighbor + random_long_range_neighbor(generator);
            v %= N;
        } while (G[u]->find(v) != G[u]->end());

        G[u]->insert(v);
        G[v]->insert(u);

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

