#include "netgraph/core/strict_multidigraph.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
using namespace netgraph::core;
int main() {
  for (bool ordered: {false,true}) {
    const int n=16384,m=131072;
    std::vector<NodeId>s(m),d(m);std::vector<Cost>w(m);std::vector<Cap>a(m,1);
    std::mt19937 rng(11);
    for(int i=0;i<m;++i){s[i]=i%n;d[i]=rng()%n;w[i]=ordered?i+1:1+rng()%1048575;}
    std::vector<double> samples;
    std::uint64_t checksum=0;
    for(int rep=-2;rep<9;++rep){
      auto t=std::chrono::steady_clock::now();
      auto g=StrictMultiDiGraph::from_arrays(n,s,d,a,w);
      auto ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();
      std::uint64_t h=0;for(auto c:g.cost_view())h=h*31+c;for(auto c:g.adj_edge_index_view())h=h*31+c;
      if(rep==-2)checksum=h;else if(h!=checksum)return 2;
      if(rep>=0)samples.push_back(ms);
    }
    auto sorted=samples;std::sort(sorted.begin(),sorted.end());
    std::cout<<(ordered?"ordered":"random")<<','<<sorted[4]<<','<<checksum<<',';
    for(auto t:samples)std::cout<<t<<';';std::cout<<std::endl;
  }
}
