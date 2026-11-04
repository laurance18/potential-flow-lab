# pfl_core: pure math/physics. Must never depend on OpenGL, GLFW or ImGui.
add_library(pfl_core STATIC
    vec2.hpp
    flow_element.hpp flow_element.cpp
    flow_field.hpp flow_field.cpp
    image_system.hpp image_system.cpp
    scene.hpp scene.cpp
    history.hpp history.cpp
    streamline.hpp streamline.cpp
    surface_pressure.hpp surface_pressure.cpp
    parameters.hpp parameters.cpp
    sweep.hpp sweep.cpp
    verification.hpp verification.cpp
    particles.hpp particles.cpp
    dye.hpp dye.cpp)
target_include_directories(pfl_core PUBLIC "${PROJECT_SOURCE_DIR}/src")
pfl_set_warnings(pfl_core)
