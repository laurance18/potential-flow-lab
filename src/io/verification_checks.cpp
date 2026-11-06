#include "io/verification_checks.hpp"

#include "core/parameters.hpp"
#include "core/streamline.hpp"
#include "io/scene_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <limits>
#include <numbers>
#include <string>
#include <system_error>
#include <vector>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;
constexpr double deg = pi / 180.0;

// Every kind of element, a rake, and plot settings that all differ from the defaults.
Scene everythingScene() {
    Scene scene;
    scene.elements = {UniformFlow{.speed = 1.3, .angle = 12.0 * deg},
                      SourceSink{.strength = 2.5, .position = {-1.5, 0.6}},
                      Vortex{.circulation = -3.0, .position = {0.7, 1.1}},
                      Doublet{.strength = 0.4, .position = {1.5, 1.4}, .angle = 35.0 * deg},
                      SectorFlow{.strength = 0.8, .vertex = {-3.0, -2.0}, .wall1 = 10.0 * deg, .wall2 = 100.0 * deg},
                      SymmetryWall{.point = {0.0, -0.5}, .angle = 7.0 * deg}};
    scene.rakes.push_back({.start = {-2.0, 0.0}, .end = {-1.0, 1.5}, .emitters = 7});
    scene.plot.view = {{-3.5, -2.2}, {4.5, 2.8}};
    scene.plot.streamlineDensity = 31;
    scene.plot.showArrows = false;
    scene.plot.arrowSpacingPx = 150.0;
    scene.plot.showImages = false;
    scene.plot.showEquipotentials = true;
    scene.plot.background = Background::PressureCoefficient;
    scene.plot.showSurfacePressure = true;
    scene.plot.showStagnationPoints = false;
    scene.plot.showParticles = false;
    scene.plot.particleTrails = false;
    return scene;
}

// Largest relative difference between the parameters of two element lists of the same kinds.
double elementDifference(const Scene& a, const Scene& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < a.elements.size(); ++i) {
        const ElementKind kind = kindOf(a.elements[i]);
        if (kind != kindOf(b.elements[i])) return std::numeric_limits<double>::infinity();
        for (const ElementParameter& parameter : elementParameters(kind)) {
            const double x = parameter.get(a.elements[i]);
            const double y = parameter.get(b.elements[i]);
            worst = std::max(worst, std::abs(x - y) / std::max(std::abs(x), 1.0));
        }
    }
    return worst;
}

// Largest distance between corresponding points of two sets of lines; infinite if they differ in shape.
double lineDifference(const StreamlineSet& a, const StreamlineSet& b) {
    if (a.lines.size() != b.lines.size()) return std::numeric_limits<double>::infinity();
    double worst = 0.0;
    for (std::size_t i = 0; i < a.lines.size(); ++i) {
        if (a.lines[i].size() != b.lines[i].size()) return std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < a.lines[i].size(); ++k) worst = std::max(worst, length(a.lines[i][k] - b.lines[i][k]));
    }
    return worst;
}

StreamlineSet linesOf(const Scene& scene) {
    StreamlineSettings settings;
    settings.bounds = scene.plot.view;
    settings.separation = scene.plot.view.size().y / scene.plot.streamlineDensity;
    settings.singularRadius = 0.02;
    return computeStreamlines(buildFlowModel(scene).field, settings);
}

bool samePlotSettings(const PlotSettings& a, const PlotSettings& b) {
    return a.streamlineDensity == b.streamlineDensity && a.showArrows == b.showArrows &&
           a.showImages == b.showImages && a.showEquipotentials == b.showEquipotentials &&
           a.background == b.background && a.showSurfacePressure == b.showSurfacePressure &&
           a.showStagnationPoints == b.showStagnationPoints && a.showParticles == b.showParticles &&
           a.particleTrails == b.particleTrails && std::abs(a.arrowSpacingPx - b.arrowSpacingPx) < 1e-12 &&
           length(a.view.min - b.view.min) < 1e-12 && length(a.view.max - b.view.max) < 1e-12;
}

VerificationResult failed(std::string why) {
    VerificationResult r;
    r.error = std::numeric_limits<double>::infinity();
    r.note = std::move(why);
    return r;
}

VerificationResult saveAndReopen() {
    const Scene saved = everythingScene();
    std::error_code ec;
    const std::filesystem::path folder = std::filesystem::temp_directory_path(ec);
    if (ec) return failed("No temporary folder to save in: " + ec.message());
    const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path path = folder / std::format("potential_flow_lab_check_{}.json", tick);

    std::string error;
    if (!saveSceneFile(saved, path, error)) return failed("Could not save: " + error);
    SceneLoadResult loaded = loadSceneFile(path);
    std::filesystem::remove(path, ec);
    if (!loaded.scene) return failed("Could not reopen: " + loaded.error);
    const Scene& reopened = *loaded.scene;
    if (reopened.elements.size() != saved.elements.size() || reopened.rakes.size() != saved.rakes.size()) {
        return failed(std::format("{} elements and {} rakes came back, expected {} and {}.", reopened.elements.size(),
                                  reopened.rakes.size(), saved.elements.size(), saved.rakes.size()));
    }

    const StreamlineSet before = linesOf(saved);
    const StreamlineSet after = linesOf(reopened);
    const double elements = elementDifference(saved, reopened);
    const double lines = lineDifference(before, after);
    const bool settings = samePlotSettings(saved.plot, reopened.plot);
    const bool rakes = saved.rakes == reopened.rakes;

    VerificationResult r;
    r.error = std::max(elements, lines);
    r.tolerance = 1e-9;
    r.errorLabel = "max relative change of any parameter, and max distance between the points of the streamlines";
    for (const std::vector<Vec2>& line : before.lines) r.samples += static_cast<int>(line.size());
    const auto shown = [](double x) { return x == 0.0 ? std::string("0") : std::format("{:.1e}", x); };
    r.comparisons = {{"Elements reopened", std::to_string(reopened.elements.size()), "6, one of each kind"},
                     {"Largest change of a parameter", shown(elements), "0"},
                     {"Rakes", rakes ? "the same" : "different", "the same"},
                     {"View and plot settings", settings ? "the same" : "different", "the same"},
                     {"Streamlines", std::format("{} lines", after.lines.size()), std::format("{} lines", before.lines.size())},
                     {"Largest distance between the plots", shown(lines), "0"}};
    r.note = "A scene with every kind of element, written to a temporary file and read back. Angles are stored in degrees, "
             "so they may come back a rounding error away.";
    if (!rakes || !settings) r.error = std::numeric_limits<double>::infinity();
    r.passed = std::isfinite(r.error) && r.error <= r.tolerance && r.samples > 0;
    return r;
}

const VerificationCheck k_saveAndReopen{
    "Saving and reopening", "Save a scene and reopen it",
    "A saved configuration must come back complete: the element list and every parameter, the walls, the rakes and the "
    "plot settings, and with them the same streamline plot.",
    "Save a scene with every kind of element to a file, reopen it, and compare everything, then the streamlines point by "
    "point.",
    everythingScene, saveAndReopen};

} // namespace

std::span<const VerificationCheck> allVerificationChecks() {
    static const std::vector<VerificationCheck> checks = [] {
        const std::span<const VerificationCheck> core = verificationChecks();
        std::vector<VerificationCheck> all(core.begin(), core.end());
        all.push_back(k_saveAndReopen);
        return all;
    }();
    return checks;
}

} // namespace pfl
