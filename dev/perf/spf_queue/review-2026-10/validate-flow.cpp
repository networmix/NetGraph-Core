#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/max_flow.hpp"
#include <cmath>
#include <iostream>
using namespace netgraph::core;
int main(){
 auto g=StrictMultiDiGraph::from_arrays(4,std::vector<NodeId>{0,0,1,1,2},std::vector<NodeId>{1,2,2,3,3},std::vector<Cap>{3,4,1,2,4},std::vector<Cost>{1,4,2,2,1});
 bool em[5];auto s=g.edge_src_view();auto d=g.edge_dst_view();auto cap=g.capacity_view();
 for(int excluded=-1;excluded<5;++excluded){
  for(int e=0;e<5;++e)em[e]=e!=excluded;
  auto [total,sm]=calc_max_flow(g,0,3,FlowPlacement::Proportional,false,true,true,true,true,{},em);
  double balance[4]={};for(int e=0;e<5;++e){double f=sm.edge_flows[e];if(f<-1e-7||f>cap[e]+1e-7||(!em[e]&&std::abs(f)>1e-7))throw std::runtime_error("capacity/mask");balance[s[e]]+=f;balance[d[e]]-=f;}
  for(int v=0;v<4;++v)if(std::abs(balance[v]-(v==0?total:v==3?-total:0))>1e-7)throw std::runtime_error("conservation");
  double cut=1e30;for(int subset=1;subset<8;subset+=2){double k=0;for(int e=0;e<5;++e)if(em[e]&&(subset&(1<<s[e]))&&!(subset&(1<<d[e])))k+=cap[e];cut=std::min(cut,k);}
  if(std::abs(cut-total)>1e-7)throw std::runtime_error("optimality");
  double reported=0;for(auto e:sm.min_cut.edges)reported+=cap[e];if(std::abs(reported-total)>1e-7)throw std::runtime_error("reported cut");
  std::cout<<"excluded="<<excluded<<" flow=cut="<<total<<'\n';
 }
 auto algs=std::make_shared<Algorithms>(make_cpu_backend());auto gh=algs->build_graph(g);ExecutionContext ctx(algs,gh);
 for(double factor:{1.,2.}){
  FlowPolicyConfig cfg;cfg.max_path_cost_factor=factor;cfg.require_capacity=true;cfg.multipath=true;
  FlowPolicy p(ctx,cfg);FlowGraph fg(g);auto res=p.place_demand(fg,0,3,1,100);
  double expected=factor==1.?2.:6.;if(std::abs(res.first-expected)>1e-7)throw std::runtime_error("factor cap");
  std::cout<<"factor="<<factor<<" placed="<<res.first<<'\n';
 }
}
