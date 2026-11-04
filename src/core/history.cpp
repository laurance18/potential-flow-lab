#include "core/history.hpp"

#include <utility>

namespace pfl {

SceneSnapshot snapshotOf(const Scene& scene) { return {scene.elements, scene.rakes}; }

void restore(Scene& scene, SceneSnapshot snapshot) {
    scene.elements = std::move(snapshot.elements);
    scene.rakes = std::move(snapshot.rakes);
}

void EditHistory::reset(SceneSnapshot current) {
    m_committed = std::move(current);
    m_undo.clear();
    m_redo.clear();
}

bool EditHistory::commit(SceneSnapshot current) {
    if (current == m_committed) return false;
    m_undo.push_back(std::exchange(m_committed, std::move(current)));
    if (m_undo.size() > m_limit) m_undo.erase(m_undo.begin());
    m_redo.clear();
    return true;
}

std::optional<SceneSnapshot> EditHistory::undo(SceneSnapshot current) {
    commit(std::move(current));
    if (m_undo.empty()) return std::nullopt;
    m_redo.push_back(std::exchange(m_committed, std::move(m_undo.back())));
    m_undo.pop_back();
    return m_committed;
}

std::optional<SceneSnapshot> EditHistory::redo(SceneSnapshot current) {
    commit(std::move(current));
    if (m_redo.empty()) return std::nullopt;
    m_undo.push_back(std::exchange(m_committed, std::move(m_redo.back())));
    m_redo.pop_back();
    return m_committed;
}

} // namespace pfl
