// Benchmark: Fibonacci baseline (kept minimal; real CFD bench is too slow for CI).
test "fib20" (b : @bench.T) {
  b.bench(fn() { b.keep(fib(20)) })
}
