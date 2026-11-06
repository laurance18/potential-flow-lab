// verification_checks.hpp: the list of checks that Verify → Physics checks shows.
//
// Every check of core/verification.hpp, then the one that needs the file code: save a scene to a
// file, reopen it, and confirm that the element list, the rakes, the plot settings and the
// streamline plot come back.
#pragma once

#include "core/verification.hpp"

#include <span>

namespace pfl {

std::span<const VerificationCheck> allVerificationChecks();

} // namespace pfl
