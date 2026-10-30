// mp4_writer_win32.hpp: H.264 video in an MP4 file, encoded by Windows Media Foundation.
//
// Media Foundation comes with Windows 10/11, except the "N" editions without the Media Feature
// Pack. Its DLLs are delay-loaded, so the program still starts there; open() then says that MP4
// is unavailable. Use one writer from one thread: open() initialises COM and Media Foundation on
// the calling thread and close() shuts them down again.
#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>

namespace pfl {

class Mp4Writer {
public:
    Mp4Writer();
    ~Mp4Writer(); // closes the file if it is still open (errors are lost; call close() to see them)
    Mp4Writer(const Mp4Writer&) = delete;
    Mp4Writer& operator=(const Mp4Writer&) = delete;

    // Creates `path` for frames of width x height pixels (both even, as H.264 needs).
    bool open(const std::filesystem::path& path, int width, int height, int fps, std::string& error);

    // Appends a frame: width * height 8-bit RGBA pixels, rows top to bottom (alpha is ignored).
    bool addFrame(const std::uint8_t* rgba, std::string& error);

    // Finishes the file (writes the MP4 index) and closes it.
    bool close(std::string& error);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace pfl
