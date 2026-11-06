// verification.hpp: the physics checks behind Verify → Physics checks.
//
// Each check builds a flow with the program's own code (elements, superposition, wall images,
// stagnation points, body surfaces, streamlines, integrated forces) and compares it with a result
// from potential-flow theory. The exact side is written out independently, in polar form, without
// calling the program's complex-potential code: V_r = Λ/2πr for a source, V_θ = -Γ/2πr for a
// vortex, Cp = 1 - 4 sin²θ on a cylinder, L' = ρ V∞ Γ, and so on. A few compare the program with
// itself in a way only correct physics can pass (a turned scene gives a turned flow).
//
// A check passes when its error (a maximum over every point, node or value compared, made
// dimensionless) is at most its tolerance. The tolerances are set by the method, not by what the
// program happens to achieve: about 1e-12 where both sides are closed-form expressions (only
// rounding separates them), looser where the program integrates (RK4 tracing, the trapezoidal rule
// over a traced surface) or the check differentiates numerically.
//
// The checks are pure functions of nothing, safe to run on any thread, and the unit tests run
// every one of them (tests/test_verification.cpp). The one that saves and reopens a scene needs the
// file code, so it is added in io/verification_checks.hpp, which has the list the window shows.
#pragma once

#include "core/scene.hpp"

#include <span>
#include <string>
#include <vector>

namespace pfl {

// One value set side by side with theory, for display: e.g. "V at (0, 1)", "(2, 0)", "(2, 0)".
struct Comparison {
    std::string quantity;
    std::string program;
    std::string exact;
};

struct VerificationResult {
    bool passed = false;
    double error = 0.0;     // dimensionless; +infinity when the check could not be made
    double tolerance = 0.0; // passes when error <= tolerance
    std::string errorLabel; // what `error` measures, e.g. "max |V − V_exact| / |V_exact|"
    int samples = 0;        // points, grid nodes or values compared
    std::vector<Comparison> comparisons;
    std::string note;       // anything else worth saying (nodes skipped, why it failed)
};

struct VerificationCheck {
    const char* group;  // e.g. "Elementary flows"
    const char* title;  // e.g. "Point vortex"
    const char* theory; // the exact result, in a sentence or two (UTF-8)
    const char* method; // how the program is put to the test
    Scene (*scene)();   // the flow checked, to show in the plot
    VerificationResult (*run)();
};

// Every check that needs only core, grouped and in display order.
std::span<const VerificationCheck> verificationChecks();

} // namespace pfl
