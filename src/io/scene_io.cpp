#include "io/scene_io.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace pfl {
namespace {

using nlohmann::json;

constexpr const char* k_format = "pfl-scene";
constexpr int k_version = 1;
constexpr double k_degrees = 180.0 / std::numbers::pi;
constexpr std::array k_panelShapes = {PanelShape::Ellipse, PanelShape::Rectangle, PanelShape::Naca, PanelShape::Custom};
constexpr std::array k_camberLines = {CamberLine::Flat, CamberLine::Naca, CamberLine::Arc};

template <typename... Ts>
struct Overloaded : Ts... {
    using Ts::operator()...;
};

json elementToJson(const SceneElement& element) {
    return std::visit(
        Overloaded{
            [](const UniformFlow& e) -> json {
                return {{"type", "uniform"}, {"speed", e.speed}, {"alpha_deg", e.angle * k_degrees}};
            },
            [](const SourceSink& e) -> json {
                return {{"type", "source"}, {"strength", e.strength}, {"x", e.position.x}, {"y", e.position.y}};
            },
            [](const Vortex& e) -> json {
                return {{"type", "vortex"}, {"circulation", e.circulation}, {"x", e.position.x}, {"y", e.position.y}};
            },
            [](const Doublet& e) -> json {
                return {{"type", "doublet"},  {"strength", e.strength},          {"x", e.position.x},
                        {"y", e.position.y}, {"beta_deg", e.angle * k_degrees}};
            },
            [](const SectorFlow& e) -> json {
                return {{"type", "sector"},
                        {"strength", e.strength},
                        {"x", e.vertex.x},
                        {"y", e.vertex.y},
                        {"theta1_deg", e.wall1 * k_degrees},
                        {"theta2_deg", e.wall2 * k_degrees}};
            },
            [](const SymmetryWall& e) -> json {
                return {{"type", "wall"}, {"x", e.point.x}, {"y", e.point.y}, {"angle_deg", e.angle * k_degrees}};
            },
            [](const PanelBody& e) -> json {
                json j = {{"type", "source_panels"},
                          {"shape", shapeName(e.shape)},
                          {"x", e.position.x},
                          {"y", e.position.y},
                          {"a", e.semiAxisA},
                          {"b", e.semiAxisB},
                          {"angle_deg", e.angle * k_degrees},
                          {"panels", e.panels}};
                if (e.shape == PanelShape::Naca) {
                    j["naca_camber"] = e.nacaCamber;
                    j["naca_position"] = e.nacaPosition;
                }
                if (e.shape == PanelShape::Custom) {
                    json points = json::array();
                    for (const Vec2 p : e.outline) points.push_back({p.x, p.y});
                    j["name"] = e.name;
                    j["smooth"] = e.smooth;
                    j["outline"] = points;
                }
                return j;
            },
            [](const ThinAirfoil& e) -> json {
                json j = {{"type", "vortex_sheet"},       {"camber", camberName(e.camber)}, {"x", e.position.x},
                          {"y", e.position.y},            {"chord", e.chord},               {"angle_deg", e.angle * k_degrees},
                          {"panels", e.panels}};
                if (e.camber == CamberLine::Naca) {
                    j["naca_camber"] = e.nacaCamber;
                    j["naca_position"] = e.nacaPosition;
                }
                if (e.camber == CamberLine::Arc) j["arc_camber"] = e.arcCamber;
                return j;
            },
        },
        element);
}

// Throws json exceptions (missing keys, wrong types) or std::runtime_error for an unknown type.
SceneElement elementFromJson(const json& j) {
    const std::string type = j.at("type").get<std::string>();
    const auto num = [&j](const char* key) { return j.at(key).get<double>(); };
    const auto point = [&num] { return Vec2{num("x"), num("y")}; };
    if (type == "uniform") return UniformFlow{.speed = num("speed"), .angle = num("alpha_deg") / k_degrees};
    if (type == "source") return SourceSink{.strength = num("strength"), .position = point()};
    if (type == "vortex") return Vortex{.circulation = num("circulation"), .position = point()};
    if (type == "doublet") {
        return Doublet{.strength = num("strength"), .position = point(), .angle = num("beta_deg") / k_degrees};
    }
    if (type == "sector") {
        return SectorFlow{.strength = num("strength"),
                          .vertex = point(),
                          .wall1 = num("theta1_deg") / k_degrees,
                          .wall2 = num("theta2_deg") / k_degrees};
    }
    if (type == "wall") return SymmetryWall{.point = point(), .angle = num("angle_deg") / k_degrees};
    if (type == "source_panels") {
        const std::string shape = j.at("shape").get<std::string>();
        PanelBody body{.position = point(),
                       .semiAxisA = num("a"),
                       .semiAxisB = num("b"),
                       .angle = j.value("angle_deg", 0.0) / k_degrees,
                       .panels = j.at("panels").get<int>()};
        const auto known = std::ranges::find_if(k_panelShapes, [&](PanelShape s) { return shape == shapeName(s); });
        if (known == k_panelShapes.end()) throw std::runtime_error("unknown panel body shape \"" + shape + "\"");
        body.shape = *known;
        if (body.shape == PanelShape::Naca) {
            body.nacaCamber = j.at("naca_camber").get<int>();
            body.nacaPosition = j.at("naca_position").get<int>();
        }
        if (body.shape == PanelShape::Custom) {
            body.name = j.value("name", std::string());
            body.smooth = j.value("smooth", true);
            for (const json& p : j.at("outline")) body.outline.push_back({p.at(0).get<double>(), p.at(1).get<double>()});
        }
        return body;
    }
    if (type == "vortex_sheet") {
        const std::string camber = j.at("camber").get<std::string>();
        ThinAirfoil airfoil{.position = point(),
                            .chord = num("chord"),
                            .angle = j.value("angle_deg", 0.0) / k_degrees,
                            .panels = j.at("panels").get<int>()};
        const auto known = std::ranges::find_if(k_camberLines, [&](CamberLine c) { return camber == camberName(c); });
        if (known == k_camberLines.end()) throw std::runtime_error("unknown camber line \"" + camber + "\"");
        airfoil.camber = *known;
        if (airfoil.camber == CamberLine::Naca) {
            airfoil.nacaCamber = j.at("naca_camber").get<int>();
            airfoil.nacaPosition = j.at("naca_position").get<int>();
        }
        if (airfoil.camber == CamberLine::Arc) airfoil.arcCamber = num("arc_camber");
        return airfoil;
    }
    throw std::runtime_error("unknown element type \"" + type + "\"");
}

constexpr std::array<std::pair<Background, const char*>, 3> k_backgrounds = {
    {{Background::None, "none"}, {Background::Speed, "speed"}, {Background::PressureCoefficient, "cp"}}};

const char* backgroundName(Background background) {
    for (const auto& [value, name] : k_backgrounds) {
        if (value == background) return name;
    }
    return "none";
}

json plotToJson(const PlotSettings& p) {
    return {{"x_min", p.view.min.x},
            {"x_max", p.view.max.x},
            {"y_min", p.view.min.y},
            {"y_max", p.view.max.y},
            {"streamline_density", p.streamlineDensity},
            {"arrows", p.showArrows},
            {"arrow_spacing_px", p.arrowSpacingPx},
            {"show_images", p.showImages},
            {"equipotentials", p.showEquipotentials},
            {"background", backgroundName(p.background)},
            {"surface_pressure", p.showSurfacePressure},
            {"stagnation_points", p.showStagnationPoints},
            {"panel_strengths", p.showPanelStrengths},
            {"particles", p.showParticles},
            {"particle_trails", p.particleTrails}};
}

// Missing keys keep their defaults, so older or hand-written files still load.
PlotSettings plotFromJson(const json& j) {
    PlotSettings p;
    p.view.min.x = j.value("x_min", p.view.min.x);
    p.view.max.x = j.value("x_max", p.view.max.x);
    p.view.min.y = j.value("y_min", p.view.min.y);
    p.view.max.y = j.value("y_max", p.view.max.y);
    p.streamlineDensity = j.value("streamline_density", p.streamlineDensity);
    p.showArrows = j.value("arrows", p.showArrows);
    p.arrowSpacingPx = j.value("arrow_spacing_px", p.arrowSpacingPx);
    p.showImages = j.value("show_images", p.showImages);
    p.showEquipotentials = j.value("equipotentials", p.showEquipotentials);
    if (j.value("speed_contours", false)) p.background = Background::Speed; // before "background" existed
    if (j.contains("background")) {
        const std::string name = j.at("background").get<std::string>();
        const auto known = std::ranges::find(k_backgrounds, name, [](const auto& entry) { return std::string(entry.second); });
        if (known == k_backgrounds.end()) throw std::runtime_error("unknown background \"" + name + "\"");
        p.background = known->first;
    }
    p.showSurfacePressure = j.value("surface_pressure", p.showSurfacePressure);
    p.showStagnationPoints = j.value("stagnation_points", p.showStagnationPoints);
    p.showPanelStrengths = j.value("panel_strengths", p.showPanelStrengths);
    p.showParticles = j.value("particles", p.showParticles);
    p.particleTrails = j.value("particle_trails", p.particleTrails);
    return p;
}

json rakeToJson(const DyeRake& r) {
    return {{"x1", r.start.x}, {"y1", r.start.y}, {"x2", r.end.x}, {"y2", r.end.y}, {"emitters", r.emitters}};
}

DyeRake rakeFromJson(const json& j) {
    const auto num = [&j](const char* key) { return j.at(key).get<double>(); };
    DyeRake r;
    r.start = {num("x1"), num("y1")};
    r.end = {num("x2"), num("y2")};
    r.emitters = j.value("emitters", r.emitters);
    return r;
}

std::optional<std::string> plotError(const PlotSettings& p) {
    const bool finite = std::isfinite(p.view.min.x) && std::isfinite(p.view.max.x) && std::isfinite(p.view.min.y) &&
                        std::isfinite(p.view.max.y) && std::isfinite(p.arrowSpacingPx);
    if (!finite || p.view.min.x >= p.view.max.x || p.view.min.y >= p.view.max.y) {
        return "plot limits must be finite with x_min < x_max and y_min < y_max";
    }
    if (p.streamlineDensity < 1 || p.arrowSpacingPx <= 0.0) {
        return "streamline_density and arrow_spacing_px must be positive";
    }
    return std::nullopt;
}

} // namespace

std::string sceneToJson(const Scene& scene) {
    json elements = json::array();
    for (const SceneElement& e : scene.elements) elements.push_back(elementToJson(e));
    json doc = {{"format", k_format}, {"version", k_version}, {"elements", elements}};
    if (!scene.rakes.empty()) {
        json rakes = json::array();
        for (const DyeRake& r : scene.rakes) rakes.push_back(rakeToJson(r));
        doc["dye_rakes"] = rakes;
    }
    doc["plot"] = plotToJson(scene.plot);
    return doc.dump(2);
}

SceneLoadResult sceneFromJson(std::string_view text) {
    SceneLoadResult result;
    const json doc = json::parse(text, nullptr, false);
    if (doc.is_discarded() || !doc.is_object()) {
        result.error = "The file is not valid JSON.";
        return result;
    }
    if (doc.value("format", std::string{}) != k_format) {
        result.error = "The file is not a Potential Flow Lab scene (missing \"format\": \"pfl-scene\").";
        return result;
    }
    if (const int version = doc.value("version", 0); version != k_version) {
        result.error = "Unsupported scene version " + std::to_string(version) + " (this program reads version " +
                       std::to_string(k_version) + ").";
        return result;
    }

    Scene scene;
    try {
        const json& elements = doc.at("elements");
        if (!elements.is_array()) throw std::runtime_error("\"elements\" must be an array");
        for (std::size_t i = 0; i < elements.size(); ++i) {
            try {
                scene.elements.push_back(elementFromJson(elements[i]));
            } catch (const std::exception& e) {
                throw std::runtime_error("element #" + std::to_string(i + 1) + ": " + e.what());
            }
            if (auto error = sceneElementError(scene.elements.back())) {
                throw std::runtime_error("element #" + std::to_string(i + 1) + " (" +
                                         kindName(kindOf(scene.elements.back())) + "): " + *error);
            }
        }
        if (doc.contains("dye_rakes")) {
            const json& rakes = doc.at("dye_rakes");
            if (!rakes.is_array()) throw std::runtime_error("\"dye_rakes\" must be an array");
            for (std::size_t i = 0; i < rakes.size(); ++i) {
                try {
                    scene.rakes.push_back(rakeFromJson(rakes[i]));
                } catch (const std::exception& e) {
                    throw std::runtime_error("dye rake #" + std::to_string(i + 1) + ": " + e.what());
                }
                if (auto error = dyeRakeError(scene.rakes.back())) {
                    throw std::runtime_error("dye rake #" + std::to_string(i + 1) + ": " + *error);
                }
            }
        }
        if (doc.contains("plot")) scene.plot = plotFromJson(doc.at("plot"));
        if (auto error = plotError(scene.plot)) throw std::runtime_error("plot settings: " + *error);
    } catch (const std::exception& e) {
        result.error = std::string("The scene could not be read: ") + e.what();
        return result;
    }
    result.scene = std::move(scene);
    return result;
}

bool saveSceneFile(const Scene& scene, const std::filesystem::path& path, std::string& error) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (file) file << sceneToJson(scene) << '\n';
    if (!file) {
        error = "Could not write the file (is the folder writable?).";
        return false;
    }
    return true;
}

SceneLoadResult loadSceneFile(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return {std::nullopt, "Could not open the file."};
    std::ostringstream text;
    text << file.rdbuf();
    return sceneFromJson(text.str());
}

} // namespace pfl
