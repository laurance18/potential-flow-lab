// parallel.hpp: a loop split across the processor's cores.
//
// parallelFor(count, grain, body) calls body(begin, end) on disjoint ranges that together cover
// [0, count), each at least `grain` long (except perhaps the last), on the standard library's
// thread pool (std::execution::par). It returns when every range is done. Small loops (count <= grain)
// run on the calling thread, so there is no overhead where there is nothing to gain.
//
// Results are deterministic as long as `body` only writes to what belongs to its own range and
// calls nothing that keeps shared state: the flow field's evaluations are pure const functions, so
// they are safe. Work that is not (a random number generator, appending to a shared list) belongs
// after the loop, done in order on one thread; then the result is the same as a plain loop's, bit for
// bit, however many threads ran.
#pragma once

#include <algorithm>
#include <cstddef>
#include <execution>
#include <numeric>
#include <thread>
#include <vector>

namespace pfl {

template <typename Body>
void parallelFor(std::size_t count, std::size_t grain, const Body& body) {
    if (count == 0) return;
    grain = std::max<std::size_t>(grain, 1);
    // A few ranges per thread, so an uneven range (particles near a body sub-step more) does not
    // keep the others waiting.
    const std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t ranges = std::min((count + grain - 1) / grain, 4 * threads);
    if (ranges <= 1) {
        body(std::size_t{0}, count);
        return;
    }
    std::vector<std::size_t> index(ranges);
    std::iota(index.begin(), index.end(), std::size_t{0});
    std::for_each(std::execution::par, index.begin(), index.end(), [&](std::size_t r) {
        body(count * r / ranges, count * (r + 1) / ranges);
    });
}

} // namespace pfl
