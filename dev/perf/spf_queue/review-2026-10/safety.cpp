#include "../spf_variants.hpp"
#include <iostream>
#include <memory>
#include <thread>
using namespace netgraph::core;
static bool eq(const auto&a,const auto&b) {return a.first==b.first&&a.second.parent_offsets==b.second.parent_offsets&&a.second.parents==b.second.parents&&a.second.via_edges==b.second.via_edges;}
static void checks() {
  for(Cost k: {Cost(0),Cost(1),Cost(65534),Cost(65535),Cost(65536),(Cost(1)<<60)-1}) {
    auto g=StrictMultiDiGraph::from_arrays(8,std::vector<NodeId>{0,0,1,2,3,4,6},std::vector<NodeId>{1,2,3,3,4,5,7},std::vector<Cap>{1,2,4,8,16,32,64},std::vector<Cost>{k,k,0,0,k,0,1});
    auto nm=std::make_unique<bool[]>(8);auto em=std::make_unique<bool[]>(7);
    std::fill_n(nm.get(),8,true);std::fill_n(em.get(),7,true);
    std::vector<Cap> res(g.capacity_view().begin(),g.capacity_view().end());
    for(int mp=0;mp<2;++mp)for(int me=0;me<2;++me)for(int tb=0;tb<2;++tb)for(int rc=0;rc<2;++rc)for(int mask=0;mask<2;++mask) {
      EdgeSelection sel;sel.multi_edge=me;sel.require_capacity=rc;sel.tie_break=tb?EdgeTieBreak::PreferHigherResidual:EdgeTieBreak::Deterministic;
      em[0]=!mask;nm[2]=!mask;
      for(int dest=-1;dest<8;++dest) {
        std::optional<NodeId> dst=dest<0?std::nullopt:std::optional<NodeId>(dest);
        auto a=shortest_paths(g,0,dst,mp,sel,res,{nm.get(),8},{em.get(),7},SpfQueue::Heap);
        auto b=shortest_paths(g,0,dst,mp,sel,res,{nm.get(),8},{em.get(),7},SpfQueue::Bucket);
        spfx::HeapQueue h;spfx::Workspace ws;auto c=spfx::spf_variant(g,0,dst,mp,sel,res,{nm.get(),8},{em.get(),7},h,ws,k,1);
        if(!eq(a,b)||!eq(a,c))throw std::runtime_error("forward mismatch");
        if(dest>=0){auto x=shortest_paths_to(g,dest,mp,sel,res,{nm.get(),8},{em.get(),7},{},SpfQueue::Heap);auto y=shortest_paths_to(g,dest,mp,sel,res,{nm.get(),8},{em.get(),7},{},SpfQueue::Bucket);if(!eq(x,y))throw std::runtime_error("reverse mismatch");}
      }
    }
  }
  for(int n:{0,1,8}){
    auto g=StrictMultiDiGraph::from_arrays(n,{},{},{},{});
    auto a=shortest_paths(g,0,std::nullopt,true,EdgeSelection{});
    auto b=shortest_paths_to(g,0,true,EdgeSelection{});
    if(a.first.size()!=std::size_t(n)||b.first.size()!=std::size_t(n)||a.second.parent_offsets.size()!=std::size_t(n+1))throw std::runtime_error("empty mismatch");
  }
}
int main(){checks();std::thread t(checks);checks();t.join();std::cout<<"PASS: boundary widths, near-int64-bound, zero costs, edgeless, forward/reverse, reused and concurrent threads\n";}
