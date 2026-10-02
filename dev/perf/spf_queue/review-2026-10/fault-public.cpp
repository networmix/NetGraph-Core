#include <cstdlib>
#include <new>
#include <iostream>
static bool fail_four = false;
void* operator new(std::size_t n) {
  if (fail_four && n == 4) { fail_four=false; throw std::bad_alloc(); }
  if (void* p=std::malloc(n?n:1)) return p;
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
#include "netgraph/core/shortest_paths.hpp"
using namespace netgraph::core;
int main(int argc,char**argv) {
  std::string mode=argc>1?argv[1]:"fault";
  EdgeSelection sel;
  if(mode=="fault") {
    auto g=StrictMultiDiGraph::from_arrays(8,std::vector<NodeId>{1},std::vector<NodeId>{0},std::vector<Cap>{1},std::vector<Cost>{1});
    fail_four=true;
    try { (void)shortest_paths(g,0,std::nullopt,true,sel); }
    catch(const std::bad_alloc&) { std::cout<<"caught forward bad_alloc\n"; }
    auto r=shortest_paths(g,1,std::nullopt,true,sel);
    std::cout<<"forward distance[0]="<<r.first[0]<<" expected=1\n";
    auto rev=StrictMultiDiGraph::from_arrays(8,std::vector<NodeId>{0},std::vector<NodeId>{1},std::vector<Cap>{1},std::vector<Cost>{1});
    fail_four=true;
    try { (void)shortest_paths_to(rev,0,true,sel); }
    catch(const std::bad_alloc&) { std::cout<<"caught reverse bad_alloc\n"; }
    auto t=shortest_paths_to(rev,1,true,sel);
    std::cout<<"reverse distance[0]="<<t.first[0]<<" expected=1\n";
  }
}
