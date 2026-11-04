// history.hpp: undo / redo of scene edits (Edit → Undo, Redo).
//
// The history holds snapshots of what the user builds: the element list and the dye rakes. The
// view and the plot settings are not part of it, so undo never pans the plot or flips a View
// option. commit() records the scene as it is now; the caller commits once an edit is finished
// (a drag released, a slider let go, a sweep stopped), so a whole gesture is one undo step. A
// commit that changes nothing adds no step, and a new step clears the redo list.
#pragma once

#include "core/dye.hpp"
#include "core/scene.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace pfl {

struct SceneSnapshot {
    std::vector<SceneElement> elements;
    std::vector<DyeRake> rakes;

    friend bool operator==(const SceneSnapshot&, const SceneSnapshot&) = default;
};

SceneSnapshot snapshotOf(const Scene& scene);
// `scene` with the elements and rakes of `snapshot`; its view and plot settings are kept.
void restore(Scene& scene, SceneSnapshot snapshot);

class EditHistory {
public:
    explicit EditHistory(std::size_t limit = 200) : m_limit(limit) {}

    // Forgets every step: `current` is the scene as loaded.
    void reset(SceneSnapshot current);

    // Records `current` as a new step if it differs from the last commit. Returns whether it did.
    bool commit(SceneSnapshot current);

    bool canUndo() const { return !m_undo.empty(); }
    bool canRedo() const { return !m_redo.empty(); }

    // `current` is committed first (an edit not yet committed is the step undone). Returns the
    // snapshot to show, or nullopt if there is nothing to undo.
    std::optional<SceneSnapshot> undo(SceneSnapshot current);
    // Returns the snapshot to show, or nullopt if there is nothing to redo (or `current` holds an
    // edit made after the last undo, which ends the redo list).
    std::optional<SceneSnapshot> redo(SceneSnapshot current);

private:
    std::size_t m_limit;
    SceneSnapshot m_committed;
    std::vector<SceneSnapshot> m_undo; // oldest first
    std::vector<SceneSnapshot> m_redo; // the next redo last
};

} // namespace pfl
