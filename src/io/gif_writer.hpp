// gif_writer.hpp: animated GIF encoder (GIF89a, looping forever), written from scratch so it
// needs no library and runs anywhere.
//
// GIF frames hold at most 256 colours, so every frame gets its own palette: a median cut of the
// colours of the pixels that changed, down to 255 colours (exact when there are fewer), with each
// colour mapped to the nearest palette entry. A pixel has changed when its colour differs from what
// the GIF already shows there by more than a small tolerance (8 of 255 in some channel). The other
// pixels are written as transparent (palette entry 255) so the previous frame shows through, and
// only the rectangle that holds the changes is stored. In a flow animation most of the plot is
// still or nearly so between frames, which keeps the files small and the still parts free of
// palette flicker.
//
// GIF delays are in hundredths of a second; they are rounded so that the total stays exact
// (30 fps: 3, 3, 4, 3, 3, 4, ...). Browsers slow down delays under 2/100 s, so keep fps <= 50.
#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <vector>

namespace pfl {

class GifWriter {
public:
    GifWriter() = default;
    ~GifWriter(); // closes the file if it is still open (errors are lost; call close() to see them)
    GifWriter(const GifWriter&) = delete;
    GifWriter& operator=(const GifWriter&) = delete;

    // Creates `path` for frames of width x height pixels. On failure returns false and sets `error`.
    bool open(const std::filesystem::path& path, int width, int height, int fps, std::string& error);

    // Appends a frame: width * height 8-bit RGBA pixels, rows top to bottom (alpha is ignored).
    bool addFrame(const std::uint8_t* rgba, std::string& error);

    // Writes the trailer and closes the file.
    bool close(std::string& error);

    int frameCount() const { return m_frames; }

private:
    bool check(std::string& error);

    std::ofstream m_file;
    int m_width = 0;
    int m_height = 0;
    int m_fps = 0;
    int m_frames = 0;
    std::vector<std::uint8_t> m_shown; // RGB of what the GIF shows after the last frame
};

// The GIF LZW data of `indices` (8-bit codes, before splitting into sub-blocks): a clear code
// first, then variable-width codes (9 to 12 bits, least significant bit first), restarting with a
// clear code when the 4096-entry table is full. Exposed for the unit tests.
std::vector<std::uint8_t> gifLzwEncode(std::span<const std::uint8_t> indices);

} // namespace pfl
