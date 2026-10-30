// video_recorder.hpp: writes recorded frames to an MP4 or animated GIF file on a background thread.
//
// The main thread hands over each captured frame (8-bit RGBA, rows top to bottom, at the capture
// size) and goes on drawing; the encoder thread scales it to the video size (averaging the pixels
// each output pixel covers) and encodes it. A few frames can wait in between. When the encoder
// falls behind, addFrame() waits for room, so no frame is ever dropped: the recording slows the
// screen down rather than lose frames.
#pragma once

#include "ui/record_dialog.hpp"

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace pfl {

class VideoRecorder {
public:
    VideoRecorder() = default;
    ~VideoRecorder(); // finishes the file
    VideoRecorder(const VideoRecorder&) = delete;
    VideoRecorder& operator=(const VideoRecorder&) = delete;

    // Creates `path` and starts the encoder thread. Frames arrive at `capture` size and are stored
    // at `video` size. On failure returns false and sets `error`.
    bool start(const std::filesystem::path& path, VideoFormat format, VideoSize capture, VideoSize video, int fps,
               std::string& error);

    // Queues one frame. False once the encoder has failed (finish() says why).
    bool addFrame(std::vector<std::uint8_t> rgba);

    // Encodes the queued frames and closes the file. False and `error` if anything failed.
    bool finish(std::string& error);

private:
    void encode(VideoFormat format, std::filesystem::path path, int fps, std::string& openError);

    VideoSize m_capture;
    VideoSize m_video;
    std::thread m_thread;
    std::mutex m_mutex;
    std::condition_variable m_changed; // a frame was queued or taken, or the recording is closing
    std::deque<std::vector<std::uint8_t>> m_queue;
    bool m_opened = false;  // the encoder thread has finished opening the file
    bool m_closing = false; // no more frames will come
    std::string m_error;    // the first error of the encoder thread
};

} // namespace pfl
