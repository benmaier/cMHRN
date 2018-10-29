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
#ifndef __MHRN_UTILITIES_H__
#define __MHRN_UTILITIES_H__

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <cstdlib>
#include <tuple>
#include <thread>
#include <iso646.h>

using namespace std;

void add_random_subgraph(
        size_t n,
        double p,
        vector < set < size_t > * > & G,
        default_random_engine & generator, 
        uniform_real_distribution<double> & distribution,
        size_t start_node = 0
        );

void add_nodes_belonging_to_this_component(
        size_t start_node,
        const vector < set < size_t > * > &G,
        set < size_t > * comp,
        vector < bool > &already_visited
       );

vector < set < size_t > * > get_components(
        const vector < set < size_t > * > &G
        );

void get_giant_component(
        vector < set < size_t > * > &G
        );

vector < set < size_t > * > get_components_from_edgelist(
        size_t N,
        vector < pair < size_t,size_t > > &edge_list
        );

vector < double > get_kleinberg_pmf(
        size_t N,
        double k,
        double mu
        );

void randomly_seed_engine(
        default_random_engine &generator
        );
#endif
