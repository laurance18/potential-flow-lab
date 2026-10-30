#include "io/gif_writer.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>

namespace pfl {
namespace {

constexpr std::size_t k_paletteColors = 255; // entry 255 is the transparent one
constexpr std::uint8_t k_transparent = 255;
constexpr int k_bins = 1 << 15; // colours are counted on 5 bits per channel

// LZW with 8-bit pixels: codes 0-255 are the pixels themselves, then the clear and end codes.
constexpr int k_minCodeSize = 8;
constexpr int k_clearCode = 1 << k_minCodeSize;
constexpr int k_endCode = k_clearCode + 1;
constexpr int k_maxCode = 4095;

int binOf(const std::uint8_t* rgba) { return (rgba[0] >> 3) << 10 | (rgba[1] >> 3) << 5 | (rgba[2] >> 3); }

// A pixel whose colour is this close (in every channel) to what the GIF already shows there is
// left as it is. Fading particle trails change nearly every pixel a little in every frame; without
// the tolerance each frame would store the whole plot. Comparing with what is shown, not with the
// previous frame, lets slow fades through once they add up.
constexpr int k_tolerance = 8;

bool closeColor(const std::uint8_t* rgba, const std::uint8_t* shown) {
    for (int a = 0; a < 3; ++a) {
        if (std::abs(rgba[a] - shown[a]) > k_tolerance) return false;
    }
    return true;
}

// Codes of `bits` bits, packed least significant bit first.
class BitWriter {
public:
    void write(int code, int bits) {
        m_buffer |= static_cast<std::uint32_t>(code) << m_count;
        m_count += bits;
        while (m_count >= 8) {
            m_bytes.push_back(static_cast<std::uint8_t>(m_buffer & 0xFFu));
            m_buffer >>= 8;
            m_count -= 8;
        }
    }

    std::vector<std::uint8_t> finish() {
        if (m_count > 0) m_bytes.push_back(static_cast<std::uint8_t>(m_buffer & 0xFFu));
        m_buffer = 0;
        m_count = 0;
        return std::move(m_bytes);
    }

private:
    std::vector<std::uint8_t> m_bytes;
    std::uint32_t m_buffer = 0;
    int m_count = 0;
};

// One occupied colour bin: how many pixels fell in it and their mean colour.
struct ColorBin {
    int bin = 0;
    std::uint32_t count = 0;
    std::array<double, 3> color{};
};

// A range of bins in the median cut, with the axis of its largest spread and the summed squared
// deviation along it (0: all one colour, nothing to split).
struct Box {
    std::size_t begin = 0;
    std::size_t end = 0;
    int axis = 0;
    double spread = 0.0;
};

void measure(Box& box, const std::vector<ColorBin>& bins) {
    double n = 0.0;
    std::array<double, 3> sum{};
    std::array<double, 3> squares{};
    for (std::size_t i = box.begin; i < box.end; ++i) {
        const double w = bins[i].count;
        n += w;
        for (int a = 0; a < 3; ++a) {
            sum[a] += w * bins[i].color[a];
            squares[a] += w * bins[i].color[a] * bins[i].color[a];
        }
    }
    box.spread = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double spread = squares[a] - sum[a] * sum[a] / n;
        if (spread > box.spread) {
            box.spread = spread;
            box.axis = a;
        }
    }
}

std::array<double, 3> meanColor(const Box& box, const std::vector<ColorBin>& bins) {
    double n = 0.0;
    std::array<double, 3> sum{};
    for (std::size_t i = box.begin; i < box.end; ++i) {
        n += bins[i].count;
        for (int a = 0; a < 3; ++a) sum[a] += bins[i].count * bins[i].color[a];
    }
    for (double& c : sum) c /= n;
    return sum;
}

// Median cut of `bins` into at most k_paletteColors boxes: the box with the largest spread is split
// at the pixel-weighted median along its widest axis until every colour has its own box or the
// palette is full.
std::vector<std::array<double, 3>> medianCut(std::vector<ColorBin>& bins) {
    std::vector<Box> boxes{{0, bins.size()}};
    measure(boxes.front(), bins);
    while (boxes.size() < k_paletteColors) {
        const auto widest = std::ranges::max_element(boxes, {}, &Box::spread);
        if (widest->spread <= 0.0) break;
        const Box box = *widest;
        const auto first = bins.begin() + static_cast<std::ptrdiff_t>(box.begin);
        const auto last = bins.begin() + static_cast<std::ptrdiff_t>(box.end);
        std::sort(first, last, [axis = box.axis](const ColorBin& a, const ColorBin& b) { return a.color[axis] < b.color[axis]; });

        std::uint64_t total = 0;
        for (auto it = first; it != last; ++it) total += it->count;
        std::size_t split = box.begin + 1;
        std::uint64_t below = 0;
        for (std::size_t i = box.begin; i + 1 < box.end; ++i) {
            below += bins[i].count;
            split = i + 1;
            if (2 * below >= total) break;
        }
        Box lower{box.begin, split};
        Box upper{split, box.end};
        measure(lower, bins);
        measure(upper, bins);
        *widest = lower;
        boxes.push_back(upper);
    }

    std::vector<std::array<double, 3>> palette;
    for (const Box& box : boxes) palette.push_back(meanColor(box, bins));
    return palette;
}

void put16(std::ofstream& file, int value) {
    const char bytes[2] = {static_cast<char>(value & 0xFF), static_cast<char>((value >> 8) & 0xFF)};
    file.write(bytes, 2);
}

void put8(std::ofstream& file, int value) { file.put(static_cast<char>(value & 0xFF)); }

} // namespace

std::vector<std::uint8_t> gifLzwEncode(std::span<const std::uint8_t> indices) {
    // next[code][pixel]: the code of the string `code` followed by `pixel`, 0 if not in the table yet.
    std::vector<std::array<std::uint16_t, 256>> next(k_maxCode + 1);
    BitWriter out;
    int codeSize = k_minCodeSize + 1;
    int lastCode = k_endCode; // the highest code in the table
    out.write(k_clearCode, codeSize);

    int current = -1; // the code of the string matched so far
    for (const std::uint8_t pixel : indices) {
        if (current < 0) {
            current = pixel;
            continue;
        }
        if (const std::uint16_t code = next[static_cast<std::size_t>(current)][pixel]; code != 0) {
            current = code;
            continue;
        }
        out.write(current, codeSize);
        next[static_cast<std::size_t>(current)][pixel] = static_cast<std::uint16_t>(++lastCode);
        if (lastCode >= (1 << codeSize)) ++codeSize;
        if (lastCode == k_maxCode) { // table full: start over
            out.write(k_clearCode, codeSize);
            for (auto& entry : next) entry.fill(0);
            codeSize = k_minCodeSize + 1;
            lastCode = k_endCode;
        }
        current = pixel;
    }
    if (current >= 0) {
        out.write(current, codeSize);
        // The decoder adds one more table entry for this last code (the encoder never needs it),
        // and widens its codes if that entry fills the current width.
        if (lastCode > k_endCode && lastCode + 1 == (1 << codeSize) && codeSize < 12) ++codeSize;
    }
    out.write(k_endCode, codeSize);
    return out.finish();
}

GifWriter::~GifWriter() {
    std::string ignored;
    if (m_file.is_open()) close(ignored);
}

bool GifWriter::open(const std::filesystem::path& path, int width, int height, int fps, std::string& error) {
    if (width < 1 || height < 1 || width > 0xFFFF || height > 0xFFFF || fps < 1) {
        error = "Invalid GIF size or frame rate.";
        return false;
    }
    m_width = width;
    m_height = height;
    m_fps = fps;
    m_frames = 0;
    m_shown.clear();
    m_file.open(path, std::ios::binary | std::ios::trunc);

    m_file.write("GIF89a", 6);
    put16(m_file, width);
    put16(m_file, height);
    put8(m_file, 0x70); // no global colour table (every frame has its own), 8-bit colour resolution
    put8(m_file, 0);    // background colour index
    put8(m_file, 0);    // square pixels
    // Netscape application extension: loop forever.
    m_file.write("\x21\xFF\x0B" "NETSCAPE2.0" "\x03\x01\x00\x00\x00", 19);
    return check(error);
}

bool GifWriter::addFrame(const std::uint8_t* rgba, std::string& error) {
    const std::size_t width = static_cast<std::size_t>(m_width);
    const std::size_t pixels = width * static_cast<std::size_t>(m_height);
    const bool first = m_shown.empty();
    if (first) m_shown.assign(3 * pixels, 0);
    std::vector<std::uint8_t> changedMask(pixels);
    for (std::size_t i = 0; i < pixels; ++i) changedMask[i] = first || !closeColor(rgba + 4 * i, m_shown.data() + 3 * i);
    const auto changed = [&changedMask](std::size_t i) { return changedMask[i] != 0; };

    // The rectangle holding every pixel that changed (all of the first frame). Nothing changed:
    // a single transparent pixel, to keep the frame's delay.
    std::size_t left = width, top = static_cast<std::size_t>(m_height), right = 0, bottom = 0;
    for (std::size_t y = 0; y < static_cast<std::size_t>(m_height); ++y) {
        for (std::size_t x = 0; x < width; ++x) {
            if (!changed(y * width + x)) continue;
            left = std::min(left, x);
            right = std::max(right, x);
            top = std::min(top, y);
            bottom = std::max(bottom, y);
        }
    }
    const bool still = left > right;
    if (still) left = right = top = bottom = 0;

    // Palette from the colours of the changed pixels.
    std::vector<std::uint32_t> counts(k_bins, 0);
    std::vector<std::array<std::uint64_t, 3>> sums(k_bins, {0, 0, 0});
    for (std::size_t y = top; y <= bottom && !still; ++y) {
        for (std::size_t x = left; x <= right; ++x) {
            const std::size_t i = y * width + x;
            if (!changed(i)) continue;
            const std::uint8_t* p = rgba + 4 * i;
            const int bin = binOf(p);
            ++counts[static_cast<std::size_t>(bin)];
            for (int a = 0; a < 3; ++a) sums[static_cast<std::size_t>(bin)][static_cast<std::size_t>(a)] += p[a];
        }
    }
    std::vector<ColorBin> bins;
    for (int bin = 0; bin < k_bins; ++bin) {
        const std::uint32_t n = counts[static_cast<std::size_t>(bin)];
        if (n == 0) continue;
        const auto& s = sums[static_cast<std::size_t>(bin)];
        bins.push_back({bin, n, {static_cast<double>(s[0]) / n, static_cast<double>(s[1]) / n, static_cast<double>(s[2]) / n}});
    }
    const std::vector<std::array<double, 3>> palette = bins.empty() ? std::vector<std::array<double, 3>>{} : medianCut(bins);

    // Every bin maps to its nearest palette colour.
    std::vector<std::uint8_t> indexOfBin(k_bins, 0);
    for (const ColorBin& b : bins) {
        double best = 1e300;
        for (std::size_t k = 0; k < palette.size(); ++k) {
            double d = 0.0;
            for (int a = 0; a < 3; ++a) d += (b.color[a] - palette[k][a]) * (b.color[a] - palette[k][a]);
            if (d < best) {
                best = d;
                indexOfBin[static_cast<std::size_t>(b.bin)] = static_cast<std::uint8_t>(k);
            }
        }
    }

    const std::size_t rectWidth = right - left + 1;
    const std::size_t rectHeight = bottom - top + 1;
    std::array<std::uint8_t, 256 * 3> table{};
    for (std::size_t k = 0; k < palette.size(); ++k) {
        for (std::size_t a = 0; a < 3; ++a) {
            table[3 * k + a] = static_cast<std::uint8_t>(std::clamp(static_cast<int>(palette[k][a] + 0.5), 0, 255));
        }
    }
    std::vector<std::uint8_t> indices;
    indices.reserve(rectWidth * rectHeight);
    for (std::size_t y = top; y <= bottom; ++y) {
        for (std::size_t x = left; x <= right; ++x) {
            const std::size_t i = y * width + x;
            if (still || !changed(i)) {
                indices.push_back(k_transparent);
                continue;
            }
            const std::uint8_t index = indexOfBin[static_cast<std::size_t>(binOf(rgba + 4 * i))];
            indices.push_back(index);
            std::copy_n(table.begin() + 3 * index, 3, m_shown.begin() + static_cast<std::ptrdiff_t>(3 * i));
        }
    }

    // Graphic control extension: delay, keep the previous frame under this one, transparency.
    const auto hundredths = [this](int frame) { return (frame * 100 + m_fps / 2) / m_fps; };
    const int delay = hundredths(m_frames + 1) - hundredths(m_frames);
    m_file.write("\x21\xF9\x04", 3);
    put8(m_file, (1 << 2) | (first ? 0 : 1));
    put16(m_file, delay);
    put8(m_file, k_transparent);
    put8(m_file, 0);

    // Image descriptor with a 256-entry local colour table.
    put8(m_file, 0x2C);
    put16(m_file, static_cast<int>(left));
    put16(m_file, static_cast<int>(top));
    put16(m_file, static_cast<int>(rectWidth));
    put16(m_file, static_cast<int>(rectHeight));
    put8(m_file, 0x80 | 0x07);
    m_file.write(reinterpret_cast<const char*>(table.data()), static_cast<std::streamsize>(table.size()));

    // Image data in sub-blocks of at most 255 bytes.
    put8(m_file, k_minCodeSize);
    const std::vector<std::uint8_t> data = gifLzwEncode(indices);
    for (std::size_t offset = 0; offset < data.size(); offset += 255) {
        const std::size_t length = std::min<std::size_t>(255, data.size() - offset);
        put8(m_file, static_cast<int>(length));
        m_file.write(reinterpret_cast<const char*>(data.data() + offset), static_cast<std::streamsize>(length));
    }
    put8(m_file, 0);

    ++m_frames;
    return check(error);
}

bool GifWriter::close(std::string& error) {
    if (!m_file.is_open()) return true;
    put8(m_file, 0x3B);
    m_file.close();
    m_shown.clear();
    return check(error);
}

bool GifWriter::check(std::string& error) {
    if (!m_file.fail()) return true;
    error = "Could not write the GIF file (is the folder writable and the disk not full?).";
    return false;
}

} // namespace pfl
