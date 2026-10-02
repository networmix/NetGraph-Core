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
#ifndef REVIEW_SOURCE
#define REVIEW_SOURCE "baseline.cpp"
#endif
#include REVIEW_SOURCE
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
  } else {
    auto bytes=[] { std::size_t n=0;for(auto& b:tls_forward.bucket.buckets)n+=b.capacity()*sizeof(BucketQueue::Item);return n; };
    for(int c=1;c<=63;++c) {
      std::vector<NodeId> s,d;std::vector<Cap> a;std::vector<Cost> w;
      for(int v=1;v<=1024;++v){s.push_back(0);d.push_back(v);a.push_back(1);w.push_back(c);}
      s.push_back(1025);d.push_back(1026);a.push_back(1);w.push_back(63);
      s.push_back(1025);d.push_back(1026);a.push_back(1);w.push_back(1);
      auto g=StrictMultiDiGraph::from_arrays(1027,s,d,a,w);
      auto r=shortest_paths(g,0,std::nullopt,true,sel);
      if(r.first[1]!=c) return 2;
      if(c==1 || c==63)std::cout<<"after "<<c<<" graphs: bucket item capacity bytes="<<bytes()<<" N=1027 E=1026 peak_queue=1024\n";
    }
    auto g=StrictMultiDiGraph::from_arrays(2,std::vector<NodeId>{0},std::vector<NodeId>{1},std::vector<Cap>{1},std::vector<Cost>{1});
    (void)shortest_paths(g,0,std::nullopt,true,sel);
    std::cout<<"after N=2 graph: bucket item capacity bytes="<<bytes()<<" dist_capacity="<<tls_forward.dist.capacity()<<"\n";
  }
}
