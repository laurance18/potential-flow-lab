#include "io/data_export.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <iterator>

namespace pfl {
namespace {

// Nodes this close (relative to the grid size) to a singular point get NaN instead of a huge value.
constexpr double k_singularFraction = 1e-6;

void appendNumber(std::string& out, double value) {
    if (std::isfinite(value)) {
        std::format_to(std::back_inserter(out), "{}", value); // shortest form that reads back exactly
    } else {
        out += "NaN";
    }
}

void appendLines(std::string& out, const char* kind, const StreamlineSet& set) {
    for (std::size_t i = 0; i < set.lines.size(); ++i) {
        for (const Vec2 p : set.lines[i]) {
            std::format_to(std::back_inserter(out), "{},{},", kind, i + 1);
            appendNumber(out, p.x);
            out += ',';
            appendNumber(out, p.y);
            out += '\n';
        }
    }
}

} // namespace

std::string fieldCsv(const FlowModel& model, const Rect& bounds, int columns, int rows) {
    const FlowField& field = model.field;
    const double nan = std::nan("");
    const double guard = k_singularFraction * length(bounds.size());
    const Vec2 step{bounds.size().x / std::max(columns - 1, 1), bounds.size().y / std::max(rows - 1, 1)};

    std::string out = "x,y,u,v,speed,phi,psi,cp\n";
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < columns; ++i) {
            const Vec2 p = bounds.min + Vec2{step.x * i, step.y * j};
            const bool defined = field.inDomain(p) && !field.nearSingularity(p, guard);
            const Vec2 v = defined ? field.velocity(p) : Vec2{nan, nan};
            const Complex f = defined ? field.complexPotential(p) : Complex{nan, nan};
            const double speed = length(v);
            const double cp = model.freeStreamSpeed > 0.0
                                  ? 1.0 - (speed * speed) / (model.freeStreamSpeed * model.freeStreamSpeed)
                                  : nan;
            for (const double value : {p.x, p.y, v.x, v.y, speed, f.real(), f.imag(), cp}) {
                appendNumber(out, value);
                out += ',';
            }
            out.back() = '\n';
        }
    }
    return out;
}

std::string linesCsv(const StreamlineSet& streamlines, const StreamlineSet& equipotentials) {
    std::string out = "kind,line,x,y\n";
    appendLines(out, "streamline", streamlines);
    appendLines(out, "equipotential", equipotentials);
    return out;
}

std::string surfaceCsv(const std::vector<BodySurface>& bodies, const FlowModel& model) {
    const double nan = std::nan("");
    std::string out = "body,closed,x,y,s,streamwise,theta_deg,cp,cp_exact\n";
    for (std::size_t b = 0; b < bodies.size(); ++b) {
        const BodySurface& body = bodies[b];
        const std::optional<CylinderReference> cylinder = cylinderReference(body, model);
        for (std::size_t k = 0; k < body.points.size(); ++k) {
            const Vec2 p = body.points[k];
            double exact = nan;
            if (cylinder) exact = cylinder->cp(std::atan2(p.y - cylinder->centre.y, p.x - cylinder->centre.x));
            std::format_to(std::back_inserter(out), "{},{},", b + 1, body.closed ? 1 : 0);
            for (const double value : {p.x, p.y, body.arcLength[k], body.streamwise[k],
                                       body.closed ? body.angle[k] : nan, body.cp[k], exact}) {
                appendNumber(out, value);
                out += ',';
            }
            out.back() = '\n';
        }
    }
    return out;
}

std::string sweepCsv(const std::vector<SweepSample>& samples) {
    const double nan = std::nan("");
    std::string out = "step,value,kind,index,x,y,closed,chord,circulation,lift_q,drag_q,cl,cd,cp_min,lift_q_exact\n";
    const auto row = [&out](std::size_t step, double value, const char* kind, std::size_t index,
                            std::initializer_list<double> values) {
        std::format_to(std::back_inserter(out), "{},", step);
        appendNumber(out, value);
        std::format_to(std::back_inserter(out), ",{},{},", kind, index);
        for (const double v : values) {
            appendNumber(out, v);
            out += ',';
        }
        out.back() = '\n';
    };
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const SweepSample& sample = samples[s];
        for (std::size_t i = 0; i < sample.stagnationPoints.size(); ++i) {
            const Vec2 p = sample.stagnationPoints[i];
            row(s + 1, sample.value, "stagnation", i + 1, {p.x, p.y, nan, nan, nan, nan, nan, nan, nan, nan, nan});
        }
        for (std::size_t i = 0; i < sample.bodies.size(); ++i) {
            const SweepBody& b = sample.bodies[i];
            row(s + 1, sample.value, "body", i + 1,
                {b.centroid.x, b.centroid.y, b.closed ? 1.0 : 0.0, b.chord, b.circulation, b.lift, b.drag,
                 b.lift / b.chord, b.drag / b.chord, b.cpMin, b.exactLift});
        }
    }
    return out;
}

bool writeTextFile(const std::filesystem::path& path, std::string_view text, std::string& error) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (file) file.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!file) {
        error = "Could not write the file (is the folder writable, or is it open in another program?).";
        return false;
    }
    return true;
}

} // namespace pfl
