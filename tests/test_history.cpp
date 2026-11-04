#include "core/history.hpp"

#include <doctest/doctest.h>

#include <cstddef>
#include <optional>
#include <variant>

using namespace pfl;

namespace {

SceneSnapshot withVortices(std::size_t count) {
    SceneSnapshot snapshot;
    for (std::size_t i = 0; i < count; ++i) {
        snapshot.elements.push_back(Vortex{.circulation = 1.0, .position = {static_cast<double>(i), 0.0}});
    }
    return snapshot;
}

} // namespace

TEST_CASE("undo and redo step through the committed edits") {
    EditHistory history;
    history.reset(withVortices(0));
    CHECK_FALSE(history.canUndo());
    CHECK(history.commit(withVortices(1)));
    CHECK(history.commit(withVortices(2)));

    CHECK(history.undo(withVortices(2)) == withVortices(1));
    CHECK(history.undo(withVortices(1)) == withVortices(0));
    CHECK_FALSE(history.undo(withVortices(0)));
    CHECK(history.canRedo());
    CHECK(history.redo(withVortices(0)) == withVortices(1));
    CHECK(history.redo(withVortices(1)) == withVortices(2));
    CHECK_FALSE(history.redo(withVortices(2)));
}

TEST_CASE("a commit that changes nothing adds no step") {
    EditHistory history;
    history.reset(withVortices(1));
    CHECK_FALSE(history.commit(withVortices(1)));
    CHECK_FALSE(history.canUndo());

    // A parameter dragged away and back again is no edit either.
    SceneSnapshot moved = withVortices(1);
    std::get<Vortex>(moved.elements[0]).position = {5.0, 5.0};
    CHECK(history.commit(moved));
    CHECK(history.commit(withVortices(1)));
    CHECK(history.undo(withVortices(1)) == moved);
}

TEST_CASE("an edit not yet committed is what undo takes back") {
    EditHistory history;
    history.reset(withVortices(0));
    history.commit(withVortices(1));
    // A drag still in progress: the scene differs from the last commit.
    SceneSnapshot dragged = withVortices(1);
    std::get<Vortex>(dragged.elements[0]).position = {2.0, 0.0};
    CHECK(history.undo(dragged) == withVortices(1));
    CHECK(history.redo(withVortices(1)) == dragged);
}

TEST_CASE("a new edit after an undo ends the redo list") {
    EditHistory history;
    history.reset(withVortices(0));
    history.commit(withVortices(1));
    history.commit(withVortices(2));
    history.undo(withVortices(2));
    CHECK(history.canRedo());

    SceneSnapshot other = withVortices(1);
    other.rakes.push_back(DyeRake{});
    CHECK(history.commit(other));
    CHECK_FALSE(history.canRedo());
    CHECK_FALSE(history.redo(other));
    CHECK(history.undo(other) == withVortices(1));

    // An uncommitted edit made after an undo ends it too.
    history.undo(withVortices(1));
    CHECK(history.canRedo());
    CHECK_FALSE(history.redo(withVortices(3)));
    CHECK(history.undo(withVortices(3)) == withVortices(0));
}

TEST_CASE("the history keeps at most its limit of steps") {
    EditHistory history(3);
    history.reset(withVortices(0));
    for (std::size_t i = 1; i <= 5; ++i) history.commit(withVortices(i));
    std::size_t undone = 0;
    SceneSnapshot current = withVortices(5);
    while (std::optional<SceneSnapshot> previous = history.undo(current)) {
        current = *previous;
        ++undone;
    }
    CHECK(undone == 3);
    CHECK(current == withVortices(2));
}

TEST_CASE("reset forgets the steps; restore keeps the view and plot settings") {
    EditHistory history;
    history.reset(withVortices(0));
    history.commit(withVortices(1));
    history.reset(withVortices(4));
    CHECK_FALSE(history.canUndo());
    CHECK_FALSE(history.canRedo());

    Scene scene;
    scene.plot.view = {{-1.0, -1.0}, {1.0, 1.0}};
    scene.plot.showEquipotentials = true;
    restore(scene, withVortices(2));
    CHECK(scene.elements.size() == 2);
    CHECK(scene.plot.view == Rect{{-1.0, -1.0}, {1.0, 1.0}});
    CHECK(scene.plot.showEquipotentials);
    CHECK(snapshotOf(scene) == withVortices(2));
}
