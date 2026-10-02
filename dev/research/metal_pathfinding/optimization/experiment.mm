// Reuse the audited phase-one generators, CPU implementations and validation.
// The original source/kernel files remain unchanged, preserving their hashes.
#define main phase_one_main
#include "../experiment.mm"
#undef main
#include <map>

struct OptimizedMetal : Metal {
  std::map<std::string, id<MTLComputePipelineState>> pipelines;
  OptimizedMetal() : Metal("dev/research/metal_pathfinding/sssp.metal") {
    NSError* error = nil;
    NSString* code = [NSString stringWithContentsOfFile:@"dev/research/metal_pathfinding/optimization/kernels.metal"
                                             encoding:NSUTF8StringEncoding error:&error];
    auto lib = [device newLibraryWithSource:code options:nil error:&error];
    if (!lib) throw std::runtime_error(error.localizedDescription.UTF8String);
    for (const char* name : {"pull64","pull32","prepare","expand","bucket_reset",
                            "bucket_min","bucket_select","dag_count64","dag_fill64",
                            "dag_count32","dag_fill32"}) {
      auto fn = [lib newFunctionWithName:[NSString stringWithUTF8String:name]];
      auto pipeline = [device newComputePipelineStateWithFunction:fn error:&error];
      if (!pipeline) throw std::runtime_error(error.localizedDescription.UTF8String);
      if (pipeline.threadExecutionWidth != 32 || pipeline.maxTotalThreadsPerThreadgroup < 256)
        throw std::runtime_error("this research configuration requires SIMD32 and 256-thread groups");
      pipelines.emplace(name, pipeline);
    }
  }
};

struct Config {
  std::string name, solver;
  bool layout = false, narrow = false, gpu_dag = false;
  int lanes = 1, delta = 0;
};

struct OptimizedGraph {
  const Case& c;
  OptimizedMetal& m;
  Pool& pool;
  GpuGraph old;
  std::vector<NodeId> sources;
  std::vector<Cost> in_costs;
  std::vector<uint8_t> in_enabled;
  id<MTLBuffer> a64,b64,a32,b32,weights32;
  id<MTLBuffer> outrow,outcol,outcost,outenabled,queueA,queueB,counts,marks,args,work,minimum;
  id<MTLBuffer> inrow,incol,incost,inids,inenabled,offsets,parents,via;
  size_t n, q, e;
  bool fits32;
  double upload_ms=0, distance_ms=0, gpu_ms=0;
  uint64_t edge_visits=0;
  int rounds=0, waits=0;
  bool fallback=false;
  struct P { uint32_t n,q,round,layout,current,epoch,lanes,delta; };

  OptimizedGraph(const Case& testcase, OptimizedMetal& metal, Pool& workers,
                 std::vector<NodeId> query_sources)
      : c(testcase),m(metal),pool(workers),old(c,m,query_sources),sources(std::move(query_sources)),
        n(c.graph.num_nodes()),q(sources.size()),e(c.graph.num_edges()) {
    auto begin=Clock::now();
    if (!n || n*q>UINT32_MAX/32 || e*q>UINT32_MAX/8) throw std::runtime_error("experiment size limit");
    auto& g=c.graph;
    Cost largest=0;
    for (Cost cost:g.cost_view()) largest=std::max(largest,cost);
    // Strictly less than the sentinel; the multiplication is bounded in uint64.
    fits32=largest < Cost(UINT32_MAX) && uint64_t(largest)*n < UINT32_MAX;
    a64=m.buffer(nullptr,n*q*8); b64=m.buffer(nullptr,n*q*8);
    a32=m.buffer(nullptr,n*q*4); b32=m.buffer(nullptr,n*q*4);
    std::vector<uint32_t> ws(e);
    auto pull_ids=c.reverse ? g.adj_edge_index_view() : g.in_adj_edge_index_view();
    for(size_t j=0;j<e;++j) ws[j]=fits32 ? uint32_t(g.cost_view()[pull_ids[j]]) : 0;
    weights32=m.buffer(ws.data(),e*4);
    auto orow=c.reverse ? g.in_row_offsets_view() : g.row_offsets_view();
    auto ocol=c.reverse ? g.in_col_indices_view() : g.col_indices_view();
    auto oids=c.reverse ? g.in_adj_edge_index_view() : g.adj_edge_index_view();
    std::vector<uint8_t> oe(e);
    for(size_t j=0;j<e;++j) {
      size_t id=oids[j];
      ws[j]=fits32 ? uint32_t(g.cost_view()[id]) : 0;
      oe[j]=c.node_mask[g.edge_src_view()[id]] && c.node_mask[g.edge_dst_view()[id]] &&
            c.edge_mask[id] && c.residual[id]>=kMinCap;
    }
    outrow=m.buffer(orow.data(),orow.size_bytes()); outcol=m.buffer(ocol.data(),ocol.size_bytes());
    outcost=m.buffer(ws.data(),e*4); outenabled=m.buffer(oe.data(),e);
    queueA=m.buffer(nullptr,n*q*4); queueB=m.buffer(nullptr,n*q*4);
    counts=m.buffer(nullptr,8); marks=m.buffer(nullptr,n*q*4);
    args=m.buffer(nullptr,12); work=m.buffer(nullptr,32); minimum=m.buffer(nullptr,4);
    in_costs.resize(e); in_enabled.resize(e);
    auto ids=g.in_adj_edge_index_view();
    for(size_t j=0;j<e;++j) {
      size_t id=ids[j];
      in_costs[j]=g.cost_view()[id];
      in_enabled[j]=c.node_mask[g.edge_src_view()[id]] && c.node_mask[g.edge_dst_view()[id]] &&
                    c.edge_mask[id] && c.residual[id]>=kMinCap;
    }
    inrow=m.buffer(g.in_row_offsets_view().data(),g.in_row_offsets_view().size_bytes());
    incol=m.buffer(g.in_col_indices_view().data(),g.in_col_indices_view().size_bytes());
    inids=m.buffer(ids.data(),ids.size_bytes()); incost=m.buffer(in_costs.data(),e*8);
    inenabled=m.buffer(in_enabled.data(),e); offsets=m.buffer(nullptr,(n+1)*q*4);
    // Reserve output capacity once; only actual DAG entries are copied to host.
    parents=m.buffer(nullptr,e*q*4); via=m.buffer(nullptr,e*q*4);
    upload_ms=old.upload_ms+ms(begin);
  }

  void wait(id<MTLCommandBuffer> cb) {
    [cb commit]; [cb waitUntilCompleted]; ++waits;
    if(cb.status==MTLCommandBufferStatusError) throw std::runtime_error(cb.error.localizedDescription.UTF8String);
    gpu_ms+=(cb.GPUEndTime-cb.GPUStartTime)*1000;
  }
  void pipeline(id<MTLComputeCommandEncoder> enc,const char* name) {
    [enc setComputePipelineState:m.pipelines.at(name)];
  }
  void dispatch(id<MTLComputeCommandEncoder> enc,size_t size) {
    [enc dispatchThreads:MTLSizeMake(size,1,1) threadsPerThreadgroup:MTLSizeMake(256,1,1)];
    [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
  }
  void one(id<MTLComputeCommandEncoder> enc) {
    [enc dispatchThreadgroups:MTLSizeMake(1,1,1) threadsPerThreadgroup:MTLSizeMake(1,1,1)];
    [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
  }

  void cpu_dags(std::vector<Result>& results) {
    auto row=c.graph.in_row_offsets_view(); auto col=c.graph.in_col_indices_view();
    auto ids=c.graph.in_adj_edge_index_view();
    auto build=[&](int query) {
      auto& d=results[query].first; auto& dag=results[query].second;
      dag.parent_offsets.resize(n+1);
      // Append in incoming-CSR order in one scan. Eligibility was cached at upload.
      dag.parents.reserve(n); dag.via_edges.reserve(n);
      for(size_t v=0;v<n;++v) {
        dag.parent_offsets[v]=int(dag.parents.size());
        for(int j=row[v];j<row[v+1];++j) {
          int u=col[j]; Cost w=in_costs[j];
          bool tight=in_enabled[j] && (c.reverse ? (d[v]!=INF && d[u]==d[v]+w) :
                                                       (d[u]!=INF && d[v]==d[u]+w));
          if(tight) { dag.parents.push_back(u); dag.via_edges.push_back(ids[j]); }
        }
      }
      dag.parent_offsets[n]=int(dag.parents.size());
    };
    if(q==1) build(0); else pool.run(int(q),build);
  }

  id<MTLBuffer> pull(Config cfg) {
    id<MTLBuffer> before=cfg.narrow ? a32 : a64, after=cfg.narrow ? b32 : b64;
    if(cfg.narrow) std::fill_n(static_cast<uint32_t*>(before.contents),n*q,UINT32_MAX);
    else std::fill_n(static_cast<Cost*>(before.contents),n*q,INF);
    for(size_t i=0;i<q;++i) if(c.node_mask[sources[i]]) {
      size_t index=cfg.layout ? size_t(sources[i])*q+i : i*n+sources[i];
      if(cfg.narrow) static_cast<uint32_t*>(before.contents)[index]=0;
      else static_cast<Cost*>(before.contents)[index]=0;
    }
    P p{uint32_t(n),uint32_t(q),0,uint32_t(cfg.layout),0,0,1,0};
    while(rounds<int(n)) {
      int steps=std::min(8,int(n)-rounds);
      std::memset(old.changed.contents,0,32);
      auto cb=[m.queue commandBuffer]; auto enc=[cb computeCommandEncoder];
      pipeline(enc,cfg.narrow ? "pull32" : "pull64");
      [enc setBuffer:old.row offset:0 atIndex:0]; [enc setBuffer:old.col offset:0 atIndex:1];
      [enc setBuffer:cfg.narrow ? weights32 : old.cost offset:0 atIndex:2];
      [enc setBuffer:old.allowed offset:0 atIndex:3]; [enc setBuffer:old.changed offset:0 atIndex:6];
      for(int j=0;j<steps;++j) {
        p.round=j; [enc setBytes:&p length:sizeof(p) atIndex:7];
        [enc setBuffer:before offset:0 atIndex:4]; [enc setBuffer:after offset:0 atIndex:5];
        dispatch(enc,n*q); std::swap(before,after);
      }
      [enc endEncoding]; wait(cb); rounds+=steps; edge_visits+=uint64_t(steps)*e*q;
      if(static_cast<uint32_t*>(old.changed.contents)[steps-1]==0) break;
    }
    return before;
  }

  id<MTLBuffer> frontier(Config cfg) {
    std::fill_n(static_cast<uint32_t*>(a32.contents),n*q,UINT32_MAX);
    std::memset(marks.contents,0,n*q*4); std::memset(counts.contents,0,8);
    uint32_t active=0;
    for(size_t query=0;query<q;++query) if(c.node_mask[sources[query]]) {
      uint32_t index=uint32_t(query*n+sources[query]);
      static_cast<uint32_t*>(a32.contents)[index]=0;
      static_cast<uint32_t*>(queueA.contents)[active++]=index;
      if(cfg.delta) static_cast<uint32_t*>(marks.contents)[index]=1;
    }
    static_cast<uint32_t*>(counts.contents)[0]=active;
    uint32_t current=0;
    id<MTLBuffer> input=queueA, output=queueB;
    P p{uint32_t(n),uint32_t(q),0,0,0,1,uint32_t(cfg.lanes),uint32_t(cfg.delta)};
    // Fail rather than returning unconverged data if an implementation bug or
    // pathological scheduler exceeds this research harness's operational limit.
    while(rounds<1000000) {
      std::memset(work.contents,0,32);
      auto cb=[m.queue commandBuffer]; auto enc=[cb computeCommandEncoder];
      for(int j=0;j<8;++j) {
        p.round=j; p.current=cfg.delta ? 0 : current; p.epoch=uint32_t(rounds+j+1);
        [enc setBytes:&p length:sizeof(p) atIndex:7];
        if(cfg.delta) {
          input=queueA; output=queueB;
          pipeline(enc,"bucket_reset");
          [enc setBuffer:counts offset:0 atIndex:0]; [enc setBuffer:minimum offset:0 atIndex:1]; one(enc);
          pipeline(enc,"bucket_min");
          [enc setBuffer:a32 offset:0 atIndex:0]; [enc setBuffer:marks offset:0 atIndex:1];
          [enc setBuffer:minimum offset:0 atIndex:2]; dispatch(enc,n*q);
          pipeline(enc,"bucket_select");
          [enc setBuffer:a32 offset:0 atIndex:0]; [enc setBuffer:marks offset:0 atIndex:1];
          [enc setBuffer:minimum offset:0 atIndex:2]; [enc setBuffer:input offset:0 atIndex:3];
          [enc setBuffer:counts offset:0 atIndex:4]; dispatch(enc,n*q);
        }
        pipeline(enc,"prepare");
        [enc setBuffer:counts offset:0 atIndex:0]; [enc setBuffer:args offset:0 atIndex:1]; one(enc);
        pipeline(enc,"expand");
        [enc setBuffer:outrow offset:0 atIndex:0]; [enc setBuffer:outcol offset:0 atIndex:1];
        [enc setBuffer:outcost offset:0 atIndex:2]; [enc setBuffer:outenabled offset:0 atIndex:3];
        [enc setBuffer:a32 offset:0 atIndex:4]; [enc setBuffer:input offset:0 atIndex:5];
        [enc setBuffer:output offset:0 atIndex:6]; [enc setBuffer:counts offset:0 atIndex:8];
        [enc setBuffer:marks offset:0 atIndex:9]; [enc setBuffer:work offset:0 atIndex:10];
        [enc dispatchThreadgroupsWithIndirectBuffer:args indirectBufferOffset:0 threadsPerThreadgroup:MTLSizeMake(256,1,1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
        if(!cfg.delta) { std::swap(input,output); current=1-current; }
      }
      [enc endEncoding]; wait(cb); rounds+=8;
      for(int j=0;j<8;++j) edge_visits+=static_cast<uint32_t*>(work.contents)[j];
      if(static_cast<uint32_t*>(counts.contents)[cfg.delta ? 0 : current]==0) return a32;
    }
    throw std::runtime_error("frontier failed to converge");
  }

  void gpu_dags(std::vector<Result>& results,id<MTLBuffer> distances,bool narrow,bool layout) {
    struct DP { uint32_t n,q,e,layout,reverse,a,b,c; } p{uint32_t(n),uint32_t(q),uint32_t(e),uint32_t(layout),uint32_t(c.reverse),0,0,0};
    auto encode=[&](bool fill) {
      auto cb=[m.queue commandBuffer]; auto enc=[cb computeCommandEncoder];
      pipeline(enc,fill ? (narrow ? "dag_fill32" : "dag_fill64") : (narrow ? "dag_count32" : "dag_count64"));
      [enc setBuffer:inrow offset:0 atIndex:0]; [enc setBuffer:incol offset:0 atIndex:1];
      [enc setBuffer:incost offset:0 atIndex:2]; [enc setBuffer:inids offset:0 atIndex:3];
      [enc setBuffer:inenabled offset:0 atIndex:4]; [enc setBuffer:distances offset:0 atIndex:5];
      [enc setBuffer:offsets offset:0 atIndex:6]; [enc setBytes:&p length:sizeof(p) atIndex:7];
      [enc setBuffer:parents offset:0 atIndex:8]; [enc setBuffer:via offset:0 atIndex:9];
      dispatch(enc,n*q); [enc endEncoding]; wait(cb);
    };
    encode(false);
    pool.run(int(q),[&](int query) {
      auto ptr=static_cast<int*>(offsets.contents)+query*(n+1);
      ptr[0]=0; std::partial_sum(ptr,ptr+n+1,ptr);
      results[query].second.parent_offsets.assign(ptr,ptr+n+1);
    });
    encode(true);
    pool.run(int(q),[&](int query) {
      auto& dag=results[query].second; int size=dag.parent_offsets.back();
      auto ps=static_cast<int*>(parents.contents)+query*e;
      auto es=static_cast<int*>(via.contents)+query*e;
      dag.parents.assign(ps,ps+size); dag.via_edges.assign(es,es+size);
    });
  }

  std::vector<Result> run(Config cfg) {
    auto begin=Clock::now(); gpu_ms=0; waits=0; rounds=0; edge_visits=0; fallback=false;
    if(cfg.solver=="old_pull" || cfg.solver=="old_group") {
      auto results=old.run(cfg.solver=="old_group",cfg.name.ends_with("serial"));
      gpu_ms=old.gpu_ms; distance_ms=old.distance_ms; rounds=old.rounds;
      waits=cfg.solver=="old_pull" ? (rounds+7)/8 : 1;
      // For persistent groups, per-query round counts can differ.
      if(cfg.solver=="old_pull") edge_visits=uint64_t(rounds)*e*q;
      else for(size_t i=0;i<q;++i) edge_visits+=uint64_t(static_cast<uint32_t*>(old.changed.contents)[i])*e;
      if(!cfg.name.ends_with("serial")) cpu_dags(results);
      return results;
    }
    if(cfg.narrow && !fits32) { cfg.narrow=false; cfg.solver="pull"; fallback=true; }
    if(cfg.solver!="pull") cfg.layout=false;
    auto data=cfg.solver=="pull" ? pull(cfg) : frontier(cfg);
    std::vector<Result> results(q);
    // Hoist Objective-C accessors out of the per-distance conversion loop.
    // Calling .contents for every element dominated large-batch host time.
    auto data32=static_cast<const uint32_t*>(data.contents);
    auto data64=static_cast<const Cost*>(data.contents);
    auto copy=[&](int query) {
      auto& dist=results[query].first; dist.resize(n);
      if(!cfg.layout && !cfg.narrow) {
        auto ptr=data64+query*n;
        std::copy_n(ptr,n,dist.begin());
      } else for(size_t v=0;v<n;++v) {
        size_t index=cfg.layout ? v*q+query : query*n+v;
        if(cfg.narrow) {
          auto d=data32[index];
          dist[v]=d==UINT32_MAX ? INF : Cost(d);
        } else dist[v]=data64[index];
      }
    };
    if(q==1) copy(0); else pool.run(int(q),copy);
    distance_ms=ms(begin);
    if(cfg.gpu_dag) gpu_dags(results,data,cfg.narrow,cfg.layout); else cpu_dags(results);
    return results;
  }
};

std::vector<Config> configurations() {
  return {
    {"old_pull_serial","old_pull"}, {"old_group_serial","old_group"},
    {"old_pull_parallel_dag","old_pull"}, {"old_group_parallel_dag","old_group"},
    {"pull64_barrier","pull"}, {"pull64_interleaved","pull",true},
    {"pull32_interleaved","pull",true,true},
    {"pull32_gpu_dag","pull",true,true,true},
    {"frontier32_vertex","frontier",false,true,false,1},
    {"frontier32_simd","frontier",false,true,false,32},
    {"frontier32_gpu_dag","frontier",false,true,true,32},
    {"bucket32_d4","frontier",false,true,false,32,4},
    {"bucket32_d16","frontier",false,true,false,32,16},
    {"bucket32_d64","frontier",false,true,false,32,64}
  };
}

Case hub_case(int n) {
  std::vector<NodeId> a,b; std::vector<Cost> w; std::vector<Cap> cap;
  auto add=[&](int u,int v,int cost) {a.push_back(u);b.push_back(v);w.push_back(cost);cap.push_back(1);};
  for(int i=0;i<n;++i) add(i,(i+1)%n,1);
  for(int i=1;i<n;++i) {add(0,i,1+i%31);add(i,0,1+(i*13)%31);}
  auto c=make_case("random",n);
  c.name="hub"; c.graph=StrictMultiDiGraph::from_arrays(n,a,b,cap,w);
  c.edge_mask=std::make_unique<bool[]>(w.size()); std::fill_n(c.edge_mask.get(),w.size(),true);
  c.residual.assign(w.size(),1); return c;
}

// A stronger comparator for the uniform positive-weight grid/fabric workloads.
// This uses the same cached eligibility and one-pass parallel DAG builder as
// the optimized GPU variants. All output is checked against native Dijkstra.
std::vector<Cost> cpu_bfs(const Case& c, const OptimizedGraph& gpu, NodeId source) {
  int n=c.graph.num_nodes();
  auto row=c.reverse ? c.graph.in_row_offsets_view() : c.graph.row_offsets_view();
  auto col=c.reverse ? c.graph.in_col_indices_view() : c.graph.col_indices_view();
  auto enabled=static_cast<const uint8_t*>(gpu.outenabled.contents);
  Cost weight=c.graph.cost_view()[0];
  std::vector<Cost> dist(n,INF);
  if(!c.node_mask[source]) return dist;
  std::vector<NodeId> queue(n); size_t head=0,tail=0;
  dist[source]=0; queue[tail++]=source;
  while(head<tail) {
    NodeId u=queue[head++];
    for(int j=row[u];j<row[u+1];++j) {
      NodeId v=col[j];
      if(enabled[j] && dist[v]==INF) {dist[v]=dist[u]+weight;queue[tail++]=v;}
    }
  }
  return dist;
}

Case boundary_case(Cost weight, bool reverse, bool blocked) {
  constexpr int n=6;
  std::vector<NodeId> a{0,0,1,1,2,3,4,0,1,2},b{1,1,2,2,3,4,5,5,1,5};
  std::vector<Cost> w(a.size(),weight);
  std::vector<Cap> cap(a.size(),1);
  auto c=make_case("random",n);
  c.name="boundary_"+std::to_string(weight)+(reverse ? "_reverse" : "")+(blocked ? "_masked" : "");
  c.graph=StrictMultiDiGraph::from_arrays(n,a,b,cap,w);
  c.edge_mask=std::make_unique<bool[]>(w.size()); std::fill_n(c.edge_mask.get(),w.size(),true);
  c.residual.assign(w.size(),1); c.reverse=reverse;
  // Competing direct/parallel edges, a self-loop, threshold-adjacent residuals,
  // and blocked sources/unreachable nodes. Eligibility is evaluated in double.
  c.residual[1]=std::nextafter(kMinCap,0.0); c.residual[3]=kMinCap;
  c.edge_mask[7]=false; c.edge_mask[9]=false;
  if(blocked) c.node_mask[2]=false;
  return c;
}

void run_case(Case c,int queries,std::string mode,int samples,OptimizedMetal& metal,Pool& pool) {
  std::vector<NodeId> sources(queries);
  std::vector<Result> expected(queries);
  for(int i=0;i<queries;++i) { sources[i]=(i*997)%c.graph.num_nodes(); expected[i]=cpu(c,sources[i]); }
  OptimizedGraph gpu(c,metal,pool,sources);
  auto configs=configurations();
  if(mode=="base") configs.resize(2);
  else if(mode=="cpu") {
    configs={{"cpu_full","cpu"},{"cpu_distance","cpu"}};
    auto costs=c.graph.cost_view();
    if(!costs.empty() && costs[0]>0 && std::all_of(costs.begin(),costs.end(),[&](Cost w){return w==costs[0];})) {
      configs.push_back({"cpu_bfs_full","cpu"}); configs.push_back({"cpu_bfs_distance","cpu"});
    }
  }
  else if(mode=="optimized") configs.erase(configs.begin(),configs.begin()+2);
  else if(mode!="all") throw std::runtime_error("bad mode");
  for(auto cfg:configs) {
    std::vector<double> times,dt,gt,work,rounds,waits;
    int fallback=0;
    for(int rep=-2;rep<samples;++rep) @autoreleasepool {
      auto begin=Clock::now(); std::vector<Result> result;
      if(mode=="cpu") {
        result.resize(queries);
        auto query=[&](int i) {
          if(cfg.name.starts_with("cpu_bfs")) result[i].first=cpu_bfs(c,gpu,sources[i]);
          else if(cfg.name=="cpu_full") result[i]=cpu(c,sources[i]);
          else result[i].first=cpu_distance(c,sources[i]);
        };
        if(queries==1) query(0); else pool.run(queries,query);
        if(cfg.name=="cpu_bfs_full") gpu.cpu_dags(result);
      } else result=gpu.run(cfg);
      double elapsed=ms(begin);
      if(cfg.name.ends_with("distance")) {
        for(int i=0;i<queries;++i) if(result[i].first!=expected[i].first) throw std::runtime_error("CPU distance mismatch");
      } else verify(expected,result);
      if(rep>=0) {
        times.push_back(elapsed); dt.push_back(mode=="cpu" ? elapsed : gpu.distance_ms);
        gt.push_back(mode=="cpu" ? 0 : gpu.gpu_ms); work.push_back(mode=="cpu" ? 0 : double(gpu.edge_visits));
        rounds.push_back(mode=="cpu" ? 0 : gpu.rounds); waits.push_back(mode=="cpu" ? 0 : gpu.waits);
        fallback=mode=="cpu" ? 0 : gpu.fallback;
      }
    }
    std::cout<<c.name<<','<<c.graph.num_nodes()<<','<<c.graph.num_edges()<<','<<queries<<','<<cfg.name
      <<','<<std::fixed<<std::setprecision(6)<<median(times)<<','<<median(dt)<<','<<median(gt)
      <<','<<uint64_t(median(work))<<','<<int(median(rounds))<<','<<int(median(waits))<<','<<fallback
      <<','<<gpu.upload_ms<<','<<samples<<",pass,";
    for(size_t j=0;j<times.size();++j) std::cout<<(j ? ";" : "")<<times[j];
    std::cout<<'\n'<<std::flush;
  }
}

int main(int argc,char** argv) {
  @autoreleasepool { try {
    if(argc<3) throw std::runtime_error("usage: optimized base|optimized|cpu|all smoke|full|correctness|stress [samples]");
    std::string mode=argv[1],suite=argv[2]; int samples=argc>3 ? std::stoi(argv[3]) : 7;
    if(samples<1) throw std::runtime_error("invalid samples");
    OptimizedMetal metal; Pool pool(10);
    std::cout<<"case,n,e,queries,method,wall_ms,distance_ms,gpu_ms,edge_visits,rounds,waits,fallback64,upload_ms,samples,correct,wall_samples_ms\n";
    auto run=[&](Case c,int q) { run_case(std::move(c),q,mode,samples,metal,pool); };
    if(suite=="correctness") {
      for(unsigned seed=1;seed<=12;++seed) run(make_case("random",37+seed*3,6,true,seed%2,seed%3==0,seed),8);
      for(Cost weight:{Cost(700000000),Cost(715827883),(Cost(1)<<53)+1})
        for(bool reverse:{false,true}) for(bool blocked:{false,true}) run(boundary_case(weight,reverse,blocked),6);
    } else if(suite=="stress") {
      run(make_case("grid",1024),64); run(make_case("random",1024,8,true,false,true),64);
    } else {
      run(make_case("random",1024),1); run(make_case("random",1024),64);
      run(make_case("clos",200),64); run(make_case("grid",1024),1);
      run(hub_case(1024),1); run(make_case("random",1024,8,true,true,true),8);
      if(suite=="full") {
        run(make_case("random",4096),256);
        run(make_case("random",16384),1); run(make_case("random",16384),64);
        run(make_case("random",65536),1);
        run(make_case("clos",800),1); run(make_case("clos",800),64);
        run(make_case("grid",10000),1); run(make_case("grid",10000),8);
        run(make_case("chain",2048),1); run(hub_case(4096),64);
        run(make_case("random",4096,32),1); run(make_case("random",4096,32),64);
      } else if(suite!="smoke") throw std::runtime_error("bad suite");
    }
  } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;} }
}
