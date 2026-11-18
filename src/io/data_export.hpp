// data_export.hpp: plain-text CSV exports for post-processing (MATLAB, Excel, Python).
//
// Field data: one row per node of a regular grid (nodes on the rectangle's edges included) with
//     x, y, u, v, speed, phi, psi, cp
// Numbers are written in the shortest form that reads back as the same double. "NaN" marks nodes
// outside the flow domain or at a singular point, and cp where no free stream defines it. phi and psi use the principal branch of log, so phi jumps by Gamma across the line
// running in -x from each vortex, and psi by Lambda across the same line from each source.
//
// Lines: one row per point of every streamline and equipotential line,
//     kind, line, x, y
// with kind "streamline" or "equipotential" and line numbered from 1 within each kind.
//
// Surface pressure: one row per point of every body surface (see core/surface_pressure.hpp),
//     body, closed, x, y, s, streamwise, theta_deg, cp, cp_exact, lambda, lambda_exact
// with body numbered from 1, closed 1 or 0, s the arc length from the stagnation point, streamwise
// x/c for a closed body (distance downstream of the stagnation point for a half-body), theta_deg the
// angle around a closed body from the free-stream direction, and cp_exact the closed-form Cp for a
// circular cylinder, or for a source-panel ellipse or circle alone in a uniform stream (whose points
// are its panels' control points). lambda is the source strength of a panel body's panel at that
// point, lambda_exact that of the exact continuous sheet where it is known. Columns that do not
// apply are NaN.
//
// Parameter sweep (see core/sweep.hpp): rows of three kinds for each step of the sweep,
//     step, value, kind, index, x, y, closed, chord, circulation, lift_q, drag_q, cl, cd, cp_min, lift_q_exact,
//     alpha_deg, cm_c4, cl_thin, cm_c4_thin
// kind "stagnation": one row per stagnation point in the view, at (x, y); the other columns NaN.
// kind "body": one row per body found, (x, y) its centroid (a half-body: its stagnation point),
// closed 1 or 0, and for a closed body its chord, the circulation inside, L'/q_inf and D'/q_inf from
// the surface Cp, c_l = L'/(q_inf c) and c_d, and the Kutta-Joukowski lift 2 Gamma / V_inf when the
// body is isolated. cp_min is the lowest surface Cp.
// kind "airfoil": one row per vortex-sheet airfoil, (x, y) the middle of its chord, its chord and
// circulation Gamma, L'/q_inf and c_l from Delta Cp along the chord, the lowest Cp on either side,
// the Kutta-Joukowski 2 Gamma / V_inf, the angle of attack in degrees, c_m about the quarter chord,
// and thin-airfoil theory's c_l and c_m,c/4 at that angle; closed, drag_q and cd are NaN.
// step and index count from 1. Steps where the value makes the element invalid have no rows. The
// value of an angle is in radians.
#pragma once

#include "core/scene.hpp"
#include "core/streamline.hpp"
#include "core/surface_pressure.hpp"
#include "core/sweep.hpp"

#include <vector>

#include <filesystem>
#include <string>
#include <string_view>

namespace pfl {

std::string fieldCsv(const FlowModel& model, const Rect& bounds, int columns, int rows);
std::string linesCsv(const StreamlineSet& streamlines, const StreamlineSet& equipotentials);
std::string surfaceCsv(const std::vector<BodySurface>& bodies, const FlowModel& model);
std::string sweepCsv(const std::vector<SweepSample>& samples);

// Writes `text` to `path`, replacing the file. On failure returns false and sets `error`.
bool writeTextFile(const std::filesystem::path& path, std::string_view text, std::string& error);

} // namespace pfl
