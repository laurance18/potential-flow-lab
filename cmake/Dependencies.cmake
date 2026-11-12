# Third-party dependencies. Everything is pinned and built from source as static libraries.
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# --- GLFW 3.5.1 (commit 70a9bb3881fe80fd483236e2b203cb451c6ecf40) ---------------------------
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
set(USE_MSVC_RUNTIME_LIBRARY_DLL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.5.1
    GIT_SHALLOW    TRUE
    SYSTEM)

# --- Dear ImGui v1.93.0-docking (commit 4e02e84cbb85d48881388e7090c2dfa1708c02a1) -----------
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.93.0-docking
    GIT_SHALLOW    TRUE)

# --- ImPlot master @ 2026-10-10 (v1.0 + fixes for ImGui 1.92.9 / 1.93.0) -----------------------
# The v1.0 tag predates ImGui's line-texture change (1.92.9) and draws invisible plot lines with
# ImGui 1.93.0; commit 7eeb916 "use ImDrawList::_SelectLineTexture()" fixes that.
FetchContent_Declare(implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        851cdbdd60d1dfb64c23e987969694de0f547a6b)

# --- nlohmann/json v3.12.0 (release tarball) -------------------------------------------------
set(PFL_JSON_VERSION 3.12.0) # also shown in Help → About
FetchContent_Declare(nlohmann_json
    URL      https://github.com/nlohmann/json/releases/download/v${PFL_JSON_VERSION}/json.tar.xz
    URL_HASH SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SYSTEM)

# --- doctest v2.5.3 (commit 2d0a9359a60c51affe2a9bebb1be1dca47868151) ------------------------
FetchContent_Declare(doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.5.3
    GIT_SHALLOW    TRUE
    SYSTEM)

FetchContent_MakeAvailable(glfw imgui implot nlohmann_json doctest)

# --- glad2 loader (GL 4.6 core + KHR_debug), generated once and committed --------------------
add_library(glad STATIC "${PROJECT_SOURCE_DIR}/third_party/glad/src/gl.c")
target_include_directories(glad SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/third_party/glad/include")

# --- Dear ImGui + backends + ImPlot as one static library ------------------------------------
add_library(pfl_imgui STATIC
    "${imgui_SOURCE_DIR}/imgui.cpp"
    "${imgui_SOURCE_DIR}/imgui_draw.cpp"
    "${imgui_SOURCE_DIR}/imgui_tables.cpp"
    "${imgui_SOURCE_DIR}/imgui_widgets.cpp"
    "${imgui_SOURCE_DIR}/imgui_demo.cpp"
    "${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp"
    "${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp"
    "${implot_SOURCE_DIR}/implot.cpp"
    "${implot_SOURCE_DIR}/implot_items.cpp")
target_include_directories(pfl_imgui SYSTEM PUBLIC
    "${imgui_SOURCE_DIR}"
    "${imgui_SOURCE_DIR}/backends"
    "${implot_SOURCE_DIR}")
target_compile_options(pfl_imgui PRIVATE /utf-8)
target_link_libraries(pfl_imgui PUBLIC glfw)
