clear all;
close all;

B = 4;
L = 3;
k = 3;
xi = 0.3;
use_giant_component = false;
delete_nodes_non_giant = false;
random_seed = 3984576;

N, edge_list = fast_mhrn(B,L,k,xi,use_giant_component,delete_nodes_non_giant,random_seed);

rows = edge_list(:,1);
cols = edge_list(:,2);

G = graph(rows,cols);
plot(G);
