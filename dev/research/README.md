# Research record

Investigations whose outcome informed the code but whose prototypes are not
part of the library. They are kept for the evidence and the dead ends they
document, not as maintained tools.

- [`metal_pathfinding/`](metal_pathfinding/): Apple Metal GPU shortest-path
  prototypes (two phases, with raw measurements and validation audits). Their
  conclusion, that a better CPU queue should come first, led to the work in
  [`../perf/spf_queue/`](../perf/spf_queue/). Objective-C++ and Metal; builds
  only on macOS with Xcode.
