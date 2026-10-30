#include "core/scene.hpp"
#include "core/verification.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <set>
#include <string>

using namespace pfl;

// The checks behind Verify → Physics checks must all pass, or the app would show a failure.
TEST_CASE("every physics check in the Verify window passes") {
    for (const VerificationCheck& check : verificationChecks()) {
        const VerificationResult result = check.run();
        INFO(check.group << " / " << check.title << ": error " << result.error << ", tolerance " << result.tolerance
                         << ", " << result.samples << " samples. " << result.note);
        CHECK(result.passed);
        CHECK(std::isfinite(result.error));
        CHECK(result.error <= result.tolerance);
        CHECK(result.samples > 0);
        CHECK_FALSE(result.comparisons.empty());
        CHECK_FALSE(result.errorLabel.empty());
    }
}

TEST_CASE("physics checks: unique titles, and every one shows a valid scene") {
    std::set<std::string> titles;
    for (const VerificationCheck& check : verificationChecks()) {
        INFO(check.title);
        CHECK(titles.insert(check.title).second);
        const Scene scene = check.scene();
        CHECK_FALSE(scene.elements.empty());
        for (const SceneElement& element : scene.elements) CHECK_FALSE(sceneElementError(element).has_value());
    }
}
