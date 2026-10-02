#include "netgraph/core/shortest_paths.hpp"
#include <iostream>
using namespace netgraph::core;
int main(){
 auto g=StrictMultiDiGraph::from_arrays(6,std::vector<NodeId>{0,0,2,0,0,4},std::vector<NodeId>{1,2,1,3,4,5},std::vector<Cap>{1,10,10,10,10,10},std::vector<Cost>{1,0,0,1,2,1});
 auto r=shortest_paths(g,0,3,true,EdgeSelection{});
 std::cout<<"target_cost="<<r.first[3]<<" beyond_target_cost="<<r.first[5]<<" (node 5 requires expanding distance-2 node 4)\n";
}
