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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Utilities.h"
#include "mhrn.h"
#include "kleinberg.h"
#include "small_world.h"
//#include "ResultClasses.h"
//#include "test.h"

using namespace std;
namespace py = pybind11;

PYBIND11_MODULE(cMHRN, m) {
    m.doc() = "Module to generate MHRN, Kleinberg networks and small world networks in a fast manner.";
    
    m.def("fast_mhrn", &fast_mhrn_edge_list, "Returns an MHRN as edge list.",
            py::arg("B"),
            py::arg("L"),
            py::arg("k"),
            py::arg("xi"),
            py::arg("use_giant_component") = false,
            py::arg("delete_non_giant_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("fast_mhrn_coord_lists", &fast_mhrn_coord_lists, "Returns an MHRN as lists of adjacency matrix coordinates.",
            py::arg("B"),
            py::arg("L"),
            py::arg("k"),
            py::arg("xi"),
            py::arg("use_giant_component") = false,
            py::arg("delete_non_giant_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("kleinberg_network", &kleinberg_edge_list, "Returns a 1d Kleinberg network (with periodic boundary conditions) as edge list. Connection probability of two nodes u and v is ~ d(u,v)^(mu-1) where d(u,v) is the pair's lattice distance. If you want to map from an MHRN, bear in mind that N=B^L and mu=log(xi)/log(B).",
            py::arg("N"),
            py::arg("k"),
            py::arg("mu"),
            py::arg("use_giant_component") = false,
            py::arg("delete_non_giant_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("kleinberg_network_coord_lists", &kleinberg_coord_lists, "Returns a 1d Kleinberg network (with periodic boundary conditions) as lists of adjacency matrix coordinates. Connection probability of two nodes u and v is ~ d(u,v)^(mu-1) where d(u,v) is the pair's lattice distance. If you want to map from an MHRN, bear in mind that N=B^L and mu=log(xi)/log(B).",
            py::arg("N"),
            py::arg("k"),
            py::arg("mu"),
            py::arg("use_giant_component") = false,
            py::arg("delete_non_giant_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("small_world_network", &small_world_edge_list, "Returns a Watts-Strogatz small world network as a pair of (number_of_nodes, edge_list). The degree k has to be an even integer.",
            py::arg("N"),
            py::arg("k"),
            py::arg("p"),
            py::arg("use_giant_component") = false,
            py::arg("delete_non_giant_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("small_world_network_coord_lists", &small_world_coord_lists, "Returns a Watts-Strogatz small world network as lists of adjacency matrix coordinates. The degree k has to be an even integer",
            py::arg("N"),
            py::arg("k"),
            py::arg("p"),
            py::arg("use_giant_component") = false,
            py::arg("delete_non_giant_component_nodes") = true,
            py::arg("seed") = 0
            );

    m.def("fast_gnp", &fast_gnp, "Returns a G(N,p) random graph in O(N+m) time as described by Batagelj and Brandes (in edge list format).",
            py::arg("N"),
            py::arg("p"),
            py::arg("start_node") = 0,
            py::arg("seed") = 0
           );

    m.def("get_components", &get_components_from_edgelist, "Get a list of sets. Each list entry is a set of nodes which make up one component of the graph.",
            py::arg("N"),
            py::arg("edge_list")
          );

    m.def("get_kleinberg_connection_probability", &get_kleinberg_pmf, "get a list of connection probabilities by neighbor distance",
            py::arg("N"),
            py::arg("k"),
            py::arg("mu")
        );

}
