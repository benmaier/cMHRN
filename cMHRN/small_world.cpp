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
#include "small_world.h"

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

tuple < size_t, vector <size_t>, vector<size_t> > small_world_coord_lists(
        size_t N,
        double k,
        double p,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        size_t seed
        )
{
    vector < set < size_t > * > G = small_world_neighbor_set(N,k,p,use_giant_component,seed);
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
    
pair < size_t, vector < pair < size_t, size_t > > > small_world_edge_list(
        size_t N,
        size_t k,
        double p,
        bool use_giant_component,
        bool delete_non_giant_component_nodes,
        size_t seed
        )
{
    size_t new_N = N;
    vector < set < size_t > * > G = small_world_neighbor_set(N,k,p,use_giant_component,seed);
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

vector < set < size_t > * > small_world_neighbor_set(
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

    // create ring
    for (size_t node = 0; node < N; ++node)
    {
        for (size_t neighbor = node+1; neighbor <= node+max_neighbor; ++neighbor)
        {
            size_t target = neighbor % N;
            G[node]->insert(target);
            G[target]->insert(node);
        }
    }

    // loop through each edge
    // for each edge pick randomly one of the neighbors as the base node
    // and reconnect to someone else
    for (size_t node = 0; node < N; ++node)
    {
        for (size_t neighbor = node+1; neighbor <= node+max_neighbor; ++neighbor)
        {
            size_t target = neighbor % N;
            size_t base;

            if (random_number(generator) < 0.5)
                base = node;
            else
                base = target;

            if (random_number(generator) < p)
            {
                size_t new_neigh;
                do
                {
                    new_neigh = random_node(generator);
                } while ( (new_neigh == base) or
                          (G[base]->find(new_neigh) != G[base]->end())
                        );

                G[node]->erase(target);
                G[target]->erase(node);

                // choose original node as base node
                G[base]->insert(new_neigh);
                G[new_neigh]->insert(base);
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

