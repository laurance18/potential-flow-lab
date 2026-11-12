# pfl_ui: Dear ImGui / ImPlot panels. Talks to core types, never to OpenGL directly.
# The user guide is built into the exe for Help → User guide (help_window.hpp).
pfl_embed_resources(pfl_docs
    BASE_DIR "${PROJECT_SOURCE_DIR}"
    NAMESPACE docs
    FILES docs/user_guide.md)

add_library(pfl_ui STATIC
    style.hpp style.cpp
    main_menu.hpp main_menu.cpp
    elements_menu.hpp elements_menu.cpp
    toolbar.hpp toolbar.cpp
    parameter_dialog.hpp parameter_dialog.cpp
    plot_view.hpp plot_view.cpp
    status_bar.hpp status_bar.cpp
    surface_panel.hpp surface_panel.cpp
    parameter_panel.hpp parameter_panel.cpp
    record_dialog.hpp record_dialog.cpp
    verification_window.hpp verification_window.cpp
    help_window.hpp help_window.cpp)
target_include_directories(pfl_ui PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(pfl_ui PUBLIC pfl_core pfl_imgui PRIVATE pfl_io pfl_docs)
pfl_set_warnings(pfl_ui)
