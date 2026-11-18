// airfoil_io.hpp: airfoil coordinate files (.dat), read into a source-panel body with an outline.
//
// Two layouts are common (the UIUC Airfoil Coordinates Database uses both):
//   - Selig: a name line, then x y pairs from the trailing edge over the upper surface to the
//     leading edge and back along the lower surface to the trailing edge;
//   - Lednicer: a name line, a line with the number of upper and lower points, then the upper surface
//     from the leading edge to the trailing edge and the lower surface the same way (often with
//     blank lines between).
// The points become the body's outline, as given (a smooth outline, which keeps the trailing edge
// sharp). The body gets a chord of 2 (a = 1) and the airfoil's own thickness, 120 panels, and the name
// from the file. A file without a name line is named after the file.
#pragma once

#include "core/panel_method.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

namespace pfl {

struct AirfoilLoadResult {
    std::optional<PanelBody> body;
    std::string error; // why it could not be read, if body is empty
};

AirfoilLoadResult airfoilFromText(std::string_view text, const std::string& fallbackName);
AirfoilLoadResult loadAirfoilFile(const std::filesystem::path& path);

} // namespace pfl
