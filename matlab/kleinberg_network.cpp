#include "CastResult.h"
#include "kleinberg.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t N, seed;
    double k, mu;
    bool use_giant_component, delete_non_giant_component_nodes;
    
    if (nrhs!=6)
    {
        mexPrintf("Got %d input arguments. ", nrhs);
        throw length_error("Invalid number of input arguments. Has to be 6.");
    }

    if (nlhs!=2)
    {
        mexPrintf("Got %d output arguments. ", nrhs);
        throw length_error("Invalid number of output arguments. Has to be 1.");
    }


    read_single_value(prhs[0],N);
    read_single_value(prhs[1],k);
    read_single_value(prhs[2],mu);
    read_single_value(prhs[3],use_giant_component);
    read_single_value(prhs[4],delete_non_giant_component_nodes);
    read_single_value(prhs[5],seed);

    pair < size_t, vector < pair < size_t, size_t > > > 
         edge_list = kleinberg_edge_list(N,k,mu,
                                         use_giant_component,
                                         delete_non_giant_component_nodes,
                                         seed
                                        );

    plhs[0] = cast_single_value(edge_list.first);
    plhs[1] = cast_edgelist_to_matlab(edge_list.second.begin(), edge_list.second.end());
}
