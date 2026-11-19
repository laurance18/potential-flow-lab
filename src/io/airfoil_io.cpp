#include "io/airfoil_io.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

namespace pfl {
namespace {

constexpr int k_importedPanels = 120;

// The two numbers on a line, if that is all it holds.
std::optional<Vec2> numberPair(const std::string& line) {
    std::istringstream in(line);
    Vec2 p;
    if (!(in >> p.x >> p.y)) return std::nullopt;
    std::string rest;
    if (in >> rest) return std::nullopt;
    return p;
}

std::string trimmed(std::string_view text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return std::string(text.substr(first, last - first + 1));
}

} // namespace

AirfoilLoadResult airfoilFromText(std::string_view text, const std::string& fallbackName) {
    std::vector<std::string> lines;
    std::istringstream in{std::string(text)};
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }

    std::string name;
    std::vector<Vec2> numbers;
    for (const std::string& line : lines) {
        if (trimmed(line).empty()) continue;
        if (const std::optional<Vec2> p = numberPair(line)) {
            numbers.push_back(*p);
        } else if (numbers.empty() && name.empty()) {
            name = trimmed(line);
        } else {
            return {std::nullopt, "This line is not a pair of numbers: \"" + trimmed(line) + "\"."};
        }
    }
    if (name.empty()) name = fallbackName;
    if (numbers.size() < 3) return {std::nullopt, "The file has fewer than 3 points."};

    std::vector<Vec2> outline;
    const Vec2 first = numbers.front();
    const bool lednicer = first.x > 1.5 && first.y > 1.5 && first.x == std::floor(first.x) && first.y == std::floor(first.y);
    if (lednicer) {
        const auto upper = static_cast<std::size_t>(first.x);
        const auto lower = static_cast<std::size_t>(first.y);
        if (numbers.size() != 1 + upper + lower) {
            return {std::nullopt, "The file says it has " + std::to_string(upper) + " upper and " + std::to_string(lower) +
                                      " lower points (Lednicer layout), but has " + std::to_string(numbers.size() - 1) +
                                      "."};
        }
        // Upper surface from the trailing edge back to the leading edge, then the lower one out again.
        outline.assign(numbers.rbegin() + static_cast<std::ptrdiff_t>(lower), numbers.rend() - 1);
        outline.insert(outline.end(), numbers.begin() + 1 + static_cast<std::ptrdiff_t>(upper), numbers.end());
    } else {
        outline = numbers;
    }

    PanelBody body{.shape = PanelShape::Custom, .panels = k_importedPanels};
    body.outline = std::move(outline);
    body.smooth = true;
    body.name = name;
    double xLo = 1e300, xHi = -1e300, yLo = 1e300, yHi = -1e300;
    for (const Vec2 p : body.outline) {
        xLo = std::min(xLo, p.x);
        xHi = std::max(xHi, p.x);
        yLo = std::min(yLo, p.y);
        yHi = std::max(yHi, p.y);
    }
    if (!(xHi > xLo) || !(yHi > yLo)) return {std::nullopt, "The points do not enclose an area."};
    body.semiAxisA = 1.0;
    body.semiAxisB = (yHi - yLo) / (xHi - xLo); // the airfoil's own proportions, chord 2
    if (const std::optional<std::string> error = panelBodyError(body)) return {std::nullopt, *error};
    return {std::move(body), {}};
}

AirfoilLoadResult loadAirfoilFile(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return {std::nullopt, "Could not open the file."};
    std::ostringstream text;
    text << file.rdbuf();
    return airfoilFromText(text.str(), path.stem().string());
}

} // namespace pfl
