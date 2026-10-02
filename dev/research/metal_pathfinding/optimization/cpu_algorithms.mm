// Stronger CPU algorithm comparators for the same positive-cost result contract.
#define main phase_one_main
#include "../experiment.mm"
#undef main

struct PreparedCpu {
  const Case& c;
  std::vector<uint8_t> enabled;
  Cost max_weight=0;
  explicit PreparedCpu(const Case& graph):c(graph),enabled(graph.graph.num_edges()) {
    const auto& g=c.graph;
    for(int e=0;e<g.num_edges();++e) {
      enabled[e]=c.node_mask[g.edge_src_view()[e]] && c.node_mask[g.edge_dst_view()[e]] &&
                 c.edge_mask[e] && c.residual[e]>=kMinCap;
      max_weight=std::max(max_weight,g.cost_view()[e]);
    }
  }
  Result solve(NodeId source,bool dial) const {
    const auto& g=c.graph; int n=g.num_nodes();
    auto row=g.row_offsets_view(),col=g.col_indices_view(),ids=g.adj_edge_index_view();
    auto costs=g.cost_view();
    Result result; auto& dist=result.first; auto& dag=result.second;
    dist.assign(n,INF);
    using Item=std::pair<Cost,NodeId>;
    auto visit=[&](Cost d,NodeId u,auto&& push) {
      if(d!=dist[u]) return;
      for(int j=row[u];j<row[u+1];++j) {
        int v=col[j],e=ids[j]; if(!enabled[e]) continue;
        Cost next=d+costs[e];
        if(next<dist[v]) {dist[v]=next;push(Item{next,v});}
      }
    };
    if(c.node_mask[source]) {
      dist[source]=0;
      if(dial) {
        if(max_weight>1024) throw std::runtime_error("Dial comparison requires max weight <= 1024");
        size_t width=size_t(max_weight)+1,pending=1;
        std::vector<std::vector<Item>> buckets(width);
        buckets[0].emplace_back(0,source); Cost current=0;
        auto push=[&](Item x){buckets[size_t(x.first)%width].push_back(x);++pending;};
        while(pending) {
          auto& bucket=buckets[size_t(current)%width];
          if(bucket.empty()) {++current;continue;}
          auto [d,u]=bucket.back();bucket.pop_back();--pending;
          if(d!=current) throw std::runtime_error("Dial ordering invariant failed");
          visit(d,u,push);
        }
      } else {
        std::priority_queue<Item,std::vector<Item>,std::greater<Item>> heap;
        heap.emplace(0,source);
        auto push=[&](Item x){heap.push(x);};
        while(!heap.empty()) {auto [d,u]=heap.top();heap.pop();visit(d,u,push);}
      }
    }
    auto inrow=g.in_row_offsets_view(),incol=g.in_col_indices_view(),inids=g.in_adj_edge_index_view();
    dag.parent_offsets.resize(n+1);dag.parents.reserve(n);dag.via_edges.reserve(n);
    for(int v=0;v<n;++v) {
      dag.parent_offsets[v]=int(dag.parents.size());
      for(int j=inrow[v];j<inrow[v+1];++j) {
        int u=incol[j],e=inids[j];
        if(enabled[e] && dist[u]!=INF && dist[v]==dist[u]+costs[e]) {
          dag.parents.push_back(u);dag.via_edges.push_back(e);
        }
      }
    }
    dag.parent_offsets[n]=int(dag.parents.size());return result;
  }
};

int main(int argc,char** argv) {
  try {
    int samples=argc>1 ? std::stoi(argv[1]) : 7; bool correctness=argc>2;
    Pool pool(10);
    std::cout<<"case,n,e,queries,method,wall_ms,samples,correct,wall_samples_ms\n";
    auto run=[&](Case c,int q) {
      std::vector<Result> expected(q); std::vector<NodeId> sources(q);
      for(int i=0;i<q;++i) {sources[i]=(i*997)%c.graph.num_nodes();expected[i]=cpu(c,sources[i]);}
      PreparedCpu prepared(c);
      for(bool dial:{false,true}) {
        std::vector<double> times;
        for(int rep=-2;rep<samples;++rep) {
          auto start=Clock::now(); std::vector<Result> actual(q);
          auto query=[&](int i){actual[i]=prepared.solve(sources[i],dial);};
          if(q==1) query(0);else pool.run(q,query);
          double elapsed=ms(start);verify(expected,actual);if(rep>=0)times.push_back(elapsed);
        }
        std::cout<<c.name<<','<<c.graph.num_nodes()<<','<<c.graph.num_edges()<<','<<q<<','
          <<(dial ? "cpu_dial_full" : "cpu_heap_cached_full")<<','<<std::fixed<<std::setprecision(6)
          <<median(times)<<','<<samples<<",pass,";
        for(size_t j=0;j<times.size();++j)std::cout<<(j ? ";" : "")<<times[j];
        std::cout<<'\n'<<std::flush;
      }
    };
    if(correctness) {
      for(unsigned seed=1;seed<=12;++seed)run(make_case("random",37+seed*3,6,true,false,false,seed),8);
    } else {
      run(make_case("random",1024),1);run(make_case("random",4096),256);
      run(make_case("random",16384),64);run(make_case("random",65536),1);
      run(make_case("random",4096,32),64);run(make_case("clos",800),64);run(make_case("grid",10000),1);
    }
  } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
