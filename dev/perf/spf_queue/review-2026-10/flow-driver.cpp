#define main original_harness_main
#include "../harness.cpp"
#undef main
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/flow_state.hpp"
static void feasible(const StrictMultiDiGraph&g,NodeId s,NodeId t,double f,std::span<const Cap> ef) {
 std::vector<double> balance(g.num_nodes());auto src=g.edge_src_view();auto dst=g.edge_dst_view();auto cap=g.capacity_view();
 for(std::size_t e=0;e<ef.size();++e){if(ef[e]<-1e-7||ef[e]>cap[e]+1e-7)throw std::runtime_error("capacity");balance[src[e]]+=ef[e];balance[dst[e]]-=ef[e];}
 for(int v=0;v<g.num_nodes();++v)if(std::abs(balance[v]-(v==s?f:v==t?-f:0))>1e-6)throw std::runtime_error("conservation");
}
int main() {
 std::cout<<"case,median_ms,hash,samples_ms\n";
 for(bool masked:{false,true})for(bool state:{false,true}) {
  auto c=make_case("random",4096,8,masked);std::vector<double> samples;std::uint64_t hash=0;
  for(int q=0;q<8;++q){c.node_mask[q*17]=true;c.node_mask[2048+q]=true;}
  for(int rep=-2;rep<7;++rep){
   std::vector<std::pair<Flow,FlowSummary>> out;out.reserve(8);
   auto start=Clock::now();
   for(int q=0;q<8;++q){
    NodeId s=q*17,t=2048+q;
    if(state){FlowState fs(c.graph);double f=fs.place_max_flow(s,t,FlowPlacement::Proportional,false,false,c.nm(),c.em());FlowSummary sm;sm.edge_flows.assign(fs.edge_flow_view().begin(),fs.edge_flow_view().end());out.emplace_back(f,std::move(sm));}
    else out.push_back(calc_max_flow(c.graph,s,t,FlowPlacement::Proportional,false,false,true,false,false,c.nm(),c.em()));
   }
   auto elapsed=ms(start);std::uint64_t h=0;
   for(int q=0;q<8;++q){auto&[f,sm]=out[q];feasible(c.graph,q*17,2048+q,f,sm.edge_flows);h=h*31+std::bit_cast<std::uint64_t>(f);for(auto e:sm.edge_flows)h=h*31+std::bit_cast<std::uint64_t>(e);for(auto e:sm.costs)h=h*31+e;for(auto e:sm.flows)h=h*31+std::bit_cast<std::uint64_t>(e);}
   if(rep==-2)hash=h;else if(hash!=h)throw std::runtime_error("changed output");if(rep>=0)samples.push_back(elapsed);
  }
  auto sorted=samples;std::sort(sorted.begin(),sorted.end());std::cout<<(state?"place_max_flow":"calc_max_flow")<<(masked?"_masked":"")<<','<<sorted[3]<<','<<hash<<',';for(auto t:samples)std::cout<<t<<';';std::cout<<std::endl;
 }
}
