#define main original_harness_main
#include "../harness.cpp"
#undef main

static std::uint64_t hash_result(const Result& r) {
  std::uint64_t h=1469598103934665603ULL;
  auto mix=[&](auto const& v) { for (auto x:v) {h^=static_cast<std::uint64_t>(x);h*=1099511628211ULL;} };
  mix(r.first); mix(r.second.parent_offsets);mix(r.second.parents);mix(r.second.via_edges);
  return h;
}
int main(int argc,char**argv) {
  std::string mode=argc>1?argv[1]:"focus";
  int samples=argc>2?std::stoi(argv[2]):7;
  if(mode=="correctness") return correctness();
  if(mode=="reverse") {
    std::cout<<"case,median_ms,hash,samples_ms\n";
    for(auto kind:{"random","clos","chain","zerocyc"}) {
      auto c=make_case(kind,std::string(kind)=="clos"?800:4096);
      std::vector<double> times;
      std::uint64_t hash=0;
      for(int rep=-2;rep<samples;++rep) {
        std::vector<Result> out;
        auto start=Clock::now();
        for(int q=0;q<16;++q) out.push_back(shortest_paths_to(c.graph,(q*997)%c.graph.num_nodes(),true,EdgeSelection{},c.residual,c.nm(),c.em()));
        auto elapsed=ms(start);
        std::uint64_t h=0;for(auto& r:out) h^=hash_result(r);
        if(rep==-2) hash=h; else if(hash!=h) throw std::runtime_error("reverse nondeterminism");
        if(rep>=0) times.push_back(elapsed);
      }
      auto sorted=times;std::sort(sorted.begin(),sorted.end());
      std::cout<<kind<<','<<sorted[sorted.size()/2]<<','<<hash<<',';
      for(auto t:times) std::cout<<t<<';';std::cout<<std::endl;
    }
    return 0;
  }
  std::vector<Workload> w;
  if(mode=="boundary") {
    for(Cost k:{32767,65534,65535,65536,131071,262143,1048575}) {
      w.push_back({"wide"+std::to_string(k),make_case("wide",16384,8,false,false,11,k),8,"none"});
      w.push_back({"near"+std::to_string(k),make_case("wide",16384,8,false,false,11,k),64,"near"});
    }
  } else {
    w.push_back({"random16",make_case("random",16384),16,"none"});
    w.push_back({"clos800",make_case("clos",800),16,"none"});
    w.push_back({"chain2048",make_case("chain",2048),256,"none"});
    w.push_back({"near16",make_case("random",16384),64,"near"});
    w.push_back({"wide20",make_case("wide",16384,8,false,false,11,1048575),16,"none"});
    w.push_back({"wide25",make_case("wide",16384,8,false,false,11,33554431),16,"none"});
  }
  time_workloads(w,"ref",samples,2097152);
}
