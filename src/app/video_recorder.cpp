#include "app/video_recorder.hpp"

#include "app/mp4_writer_win32.hpp"
#include "io/gif_writer.hpp"

#include <algorithm>
#include <array>

namespace pfl {
namespace {

constexpr std::size_t k_maxQueued = 4;

std::size_t pixelCount(VideoSize size) {
    return static_cast<std::size_t>(size.width) * static_cast<std::size_t>(size.height);
}

// `rgba` scaled from `from` to `to`: each output pixel is the mean of the input pixels it covers.
std::vector<std::uint8_t> resize(std::vector<std::uint8_t> rgba, VideoSize from, VideoSize to) {
    if (from.width == to.width && from.height == to.height) return rgba;
    std::vector<std::uint8_t> out(4 * pixelCount(to));
    for (int y = 0; y < to.height; ++y) {
        const int y0 = y * from.height / to.height;
        const int y1 = std::max(y0 + 1, (y + 1) * from.height / to.height);
        for (int x = 0; x < to.width; ++x) {
            const int x0 = x * from.width / to.width;
            const int x1 = std::max(x0 + 1, (x + 1) * from.width / to.width);
            std::array<unsigned, 4> sum{};
            for (int yy = y0; yy < y1; ++yy) {
                const std::uint8_t* row = rgba.data() + 4 * (static_cast<std::size_t>(yy) * static_cast<std::size_t>(from.width));
                for (int xx = x0; xx < x1; ++xx) {
                    for (std::size_t c = 0; c < 4; ++c) sum[c] += row[4 * static_cast<std::size_t>(xx) + c];
                }
            }
            const unsigned n = static_cast<unsigned>((y1 - y0) * (x1 - x0));
            std::uint8_t* pixel = out.data() + 4 * (static_cast<std::size_t>(y) * static_cast<std::size_t>(to.width) + static_cast<std::size_t>(x));
            for (std::size_t c = 0; c < 4; ++c) pixel[c] = static_cast<std::uint8_t>((sum[c] + n / 2) / n);
        }
    }
    return out;
}

} // namespace

VideoRecorder::~VideoRecorder() {
    std::string ignored;
    finish(ignored);
}

bool VideoRecorder::start(const std::filesystem::path& path, VideoFormat format, VideoSize capture, VideoSize video,
                          int fps, std::string& error) {
    m_capture = capture;
    m_video = video;
    std::string openError;
    m_thread = std::thread(&VideoRecorder::encode, this, format, path, fps, std::ref(openError));
    std::unique_lock lock(m_mutex);
    m_changed.wait(lock, [this] { return m_opened; });
    if (openError.empty()) return true;
    lock.unlock();
    m_thread.join();
    error = openError;
    return false;
}

bool VideoRecorder::addFrame(std::vector<std::uint8_t> rgba) {
    {
        std::unique_lock lock(m_mutex);
        m_changed.wait(lock, [this] { return m_queue.size() < k_maxQueued || !m_error.empty(); });
        if (!m_error.empty()) return false;
        m_queue.push_back(std::move(rgba));
    }
    m_changed.notify_all();
    return true;
}

bool VideoRecorder::finish(std::string& error) {
    if (m_thread.joinable()) {
        {
            std::lock_guard lock(m_mutex);
            m_closing = true;
        }
        m_changed.notify_all();
        m_thread.join();
    }
    if (m_error.empty()) return true;
    error = m_error;
    return false;
}

// The encoder thread: opens the file (Media Foundation needs that on this thread), then encodes
// frames until finish() closes the recording or a frame fails.
void VideoRecorder::encode(VideoFormat format, std::filesystem::path path, int fps, std::string& openError) {
    GifWriter gif;
    Mp4Writer mp4;
    const bool isGif = format == VideoFormat::Gif;
    std::string error;
    const bool opened = isGif ? gif.open(path, m_video.width, m_video.height, fps, error)
                              : mp4.open(path, m_video.width, m_video.height, fps, error);
    {
        std::lock_guard lock(m_mutex);
        if (!opened) openError = error;
        m_opened = true;
    }
    m_changed.notify_all();
    if (!opened) return;

    for (;;) {
        std::vector<std::uint8_t> frame;
        {
            std::unique_lock lock(m_mutex);
            m_changed.wait(lock, [this] { return !m_queue.empty() || m_closing; });
            if (m_queue.empty()) break;
            frame = std::move(m_queue.front());
            m_queue.pop_front();
        }
        m_changed.notify_all();
        const std::vector<std::uint8_t> scaled = resize(std::move(frame), m_capture, m_video);
        if (!(isGif ? gif.addFrame(scaled.data(), error) : mp4.addFrame(scaled.data(), error))) {
            {
                std::lock_guard lock(m_mutex);
                m_error = error;
            }
            m_changed.notify_all();
            break;
        }
    }
    if (!(isGif ? gif.close(error) : mp4.close(error))) {
        std::lock_guard lock(m_mutex);
        if (m_error.empty()) m_error = error;
    }
}

} // namespace pfl
