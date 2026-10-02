#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <iostream>

int main() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) return 1;
    std::cout << "device=" << device.name.UTF8String
              << " unified=" << device.hasUnifiedMemory
              << " max_buffer=" << device.maxBufferLength << '\n';
    NSArray<NSString*>* bodies = @[
      @"kernel void test(device long* x [[buffer(0)]], uint i [[thread_position_in_grid]]) { x[i] += 1; }",
      @"kernel void test(device double* x [[buffer(0)]], uint i [[thread_position_in_grid]]) { x[i] += 1.0; }",
      @"kernel void test(device atomic_ulong* x [[buffer(0)]]) { atomic_fetch_min_explicit(x, 1ul, memory_order_relaxed); }",
      @"kernel void test(device atomic_uint* x [[buffer(0)]]) { atomic_fetch_min_explicit(x, 1u, memory_order_relaxed); }",
      @"kernel void test(device atomic_ulong* x [[buffer(0)]]) { atomic_min_explicit(x, 1ul, memory_order_relaxed); }",
      @"kernel void test(device atomic_ulong* x [[buffer(0)]], device ulong* y [[buffer(1)]]) { y[0] = atomic_load_explicit(x, memory_order_relaxed); }"
    ];
    for (NSUInteger i = 0; i < bodies.count; ++i) {
      NSError* error = nil;
      NSString* source = [@"#include <metal_stdlib>\nusing namespace metal;\n" stringByAppendingString:bodies[i]];
      id<MTLLibrary> library = [device newLibraryWithSource:source options:nil error:&error];
      std::cout << "probe=" << i << " compiled=" << (library != nil) << '\n';
      if (!library) std::cout << error.localizedDescription.UTF8String << '\n';
    }
  }
}
