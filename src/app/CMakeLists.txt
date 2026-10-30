# PotentialFlowLab: the executable. Window, main loop, and wiring of core/render/ui/io.
set(PFL_ICON_PATH "${PROJECT_SOURCE_DIR}/resources/app.ico")
configure_file("${PROJECT_SOURCE_DIR}/resources/app.rc.in" "${CMAKE_CURRENT_BINARY_DIR}/app.rc" @ONLY)

add_executable(PotentialFlowLab
    main.cpp
    app.hpp app.cpp
    platform_win32.hpp platform_win32.cpp
    mp4_writer_win32.hpp mp4_writer_win32.cpp
    video_recorder.hpp video_recorder.cpp
    "${CMAKE_CURRENT_BINARY_DIR}/app.rc")
target_link_libraries(PotentialFlowLab PRIVATE pfl_core pfl_render pfl_ui pfl_io pfl_imgui glfw comdlg32 ole32 windowscodecs
                                               mfplat mfreadwrite mfuuid delayimp)
pfl_set_warnings(PotentialFlowLab)

# Media Foundation (MP4 export) is missing on Windows "N" editions: load it only when a video is
# recorded, so the program still starts there (mp4_writer_win32.cpp checks for it first).
target_link_options(PotentialFlowLab PRIVATE /DELAYLOAD:mfplat.dll /DELAYLOAD:mfreadwrite.dll)

# Release: GUI subsystem (no console window). Debug keeps the console for logs.
# mainCRTStartup lets us keep a plain int main() in both.
target_link_options(PotentialFlowLab PRIVATE
    "$<$<NOT:$<CONFIG:Debug>>:/SUBSYSTEM:WINDOWS;/ENTRY:mainCRTStartup>")

set_property(DIRECTORY "${PROJECT_SOURCE_DIR}" PROPERTY VS_STARTUP_PROJECT PotentialFlowLab)
set_target_properties(PotentialFlowLab PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
