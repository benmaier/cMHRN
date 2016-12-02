#include "CastResult.h"
#include "mhrn.h"
#include "math.h"
#include "matrix.h"
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t B, L, seed;
    double k, xi;
    bool use_giant_component, delete_non_giant_component_nodes;
    
    if (nrhs!=7)
    {
        mexPrintf("Got %d input arguments. ", nrhs);
        throw length_error("Invalid number of input arguments. Has to be 12.");
    }

    if (nlhs!=1)
    {
        mexPrintf("Got %d output arguments. ", nrhs);
        throw length_error("Invalid number of output arguments. Has to be 1.");
    }


    read_single_value(prhs[0],B);
    read_single_value(prhs[1],L);
    read_single_value(prhs[2],k);
    read_single_value(prhs[3],xi);
    read_single_value(prhs[4],use_giant_component);
    read_single_value(prhs[5],delete_non_giant_component_nodes);
    read_single_value(prhs[6],seed);

    vector < pair < size_t, size_t > > edge_list = fast_mhrn_edge_list(B,L,k,xi,
                                                                       use_giant_component,
                                                                       delete_non_giant_component_nodes,
                                                                       seed
                                                                      );

    plhs[0] = cast_edgelist_to_matlab(edge_list.begin(), edge_list.end());
}
