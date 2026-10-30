// record_dialog.hpp: the "Record animation" dialog.
//
// Chooses the file format (MP4 video or animated GIF), the length, frame rate and size of the
// recording, whether it starts from fresh particles and dye, and whether the parameter sweep
// plays once over it. The app then records the plot, axes included, as it is drawn: the scene can
// still be edited (elements dragged, sliders moved) while it records.
#pragma once

#include <imgui.h>

#include <optional>
#include <string>

namespace pfl {

enum class VideoFormat { Mp4, Gif };

struct RecordSettings {
    VideoFormat format = VideoFormat::Mp4;
    double seconds = 6.0;
    int fps = 30;
    int scalePercent = 100;
    bool restart = false; // respawn the particles and clear the dye first
    bool sweep = false;   // play the parameter sweep once, From to To, over the recording
};

// Size of the video for a plot of width x height pixels: scaled, and even for H.264.
struct VideoSize {
    int width = 0;
    int height = 0;
};
VideoSize videoSize(int width, int height, const RecordSettings& settings);

int frameCount(const RecordSettings& settings);

class RecordDialog {
public:
    // `sweep` describes the sweep the recording can play (e.g. "#3 Point vortex: Γ from 0 to 12.6"),
    // empty if none is set up in the parameter panel.
    void open(std::string sweep);

    // Draws the dialog while it is open; returns the settings when the user presses Record.
    // `plotPixels` is the size of the area that will be recorded.
    std::optional<RecordSettings> draw(ImVec2 plotPixels);

private:
    RecordSettings m_settings;
    std::string m_sweep;
    bool m_openRequested = false;
};

} // namespace pfl
