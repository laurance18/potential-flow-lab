#include "io/gif_writer.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <vector>

using namespace pfl;

namespace {

// GIF LZW decoder (8-bit codes), written from the GIF89a specification, independently of the encoder.
std::vector<std::uint8_t> lzwDecode(const std::vector<std::uint8_t>& data) {
    std::vector<std::vector<std::uint8_t>> table;
    const auto reset = [&table] {
        table.assign(258, {}); // 256 = clear, 257 = end
        for (int i = 0; i < 256; ++i) table[static_cast<std::size_t>(i)] = {static_cast<std::uint8_t>(i)};
    };
    reset();
    int codeSize = 9;
    std::size_t bit = 0;
    int previous = -1;
    std::vector<std::uint8_t> out;
    for (;;) {
        REQUIRE(bit <= 8 * data.size() + 16); // ran past the data without an end code
        int code = 0;
        for (int b = 0; b < codeSize; ++b, ++bit) {
            const std::size_t byte = bit / 8;
            const int value = byte < data.size() ? (data[byte] >> (bit % 8)) & 1 : 0;
            code |= value << b;
        }
        if (code == 256) {
            reset();
            codeSize = 9;
            previous = -1;
            continue;
        }
        if (code == 257) break;
        std::vector<std::uint8_t> entry;
        if (static_cast<std::size_t>(code) < table.size()) {
            REQUIRE((code < 256 || code > 257));
            entry = table[static_cast<std::size_t>(code)];
        } else {
            REQUIRE(static_cast<std::size_t>(code) == table.size()); // the string being defined
            REQUIRE(previous >= 0);
            entry = table[static_cast<std::size_t>(previous)];
            entry.push_back(entry.front());
        }
        out.insert(out.end(), entry.begin(), entry.end());
        if (previous >= 0 && table.size() < 4096) {
            std::vector<std::uint8_t> added = table[static_cast<std::size_t>(previous)];
            added.push_back(entry.front());
            table.push_back(std::move(added));
        }
        previous = code;
        if (table.size() == (std::size_t{1} << codeSize) && codeSize < 12) ++codeSize;
    }
    return out;
}

struct DecodedGif {
    int width = 0;
    int height = 0;
    bool loops = false;
    std::vector<int> delays;                         // hundredths of a second
    std::vector<std::array<int, 4>> rects;           // left, top, width, height
    std::vector<std::vector<std::uint8_t>> frames;   // RGBA after composing each frame
};

int read16(const std::vector<std::uint8_t>& f, std::size_t at) { return f[at] | f[at + 1] << 8; }

// Reads the subset of GIF that GifWriter writes: no global table, local 256-colour tables,
// disposal "keep", optional transparency.
DecodedGif decodeGif(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    const std::vector<std::uint8_t> f((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    REQUIRE(f.size() > 13);
    REQUIRE(std::string(f.begin(), f.begin() + 6) == "GIF89a");
    DecodedGif gif;
    gif.width = read16(f, 6);
    gif.height = read16(f, 8);
    REQUIRE((f[10] & 0x80) == 0);
    std::vector<std::uint8_t> canvas(4 * static_cast<std::size_t>(gif.width * gif.height), 0);
    std::size_t at = 13;
    int transparent = -1;
    for (;;) {
        REQUIRE(at < f.size());
        const std::uint8_t block = f[at++];
        if (block == 0x3B) break;
        if (block == 0x21) {
            const std::uint8_t label = f[at++];
            if (label == 0xF9) {
                REQUIRE(f[at] == 4);
                REQUIRE(((f[at + 1] >> 2) & 7) == 1); // keep the previous frame
                gif.delays.push_back(read16(f, at + 2));
                transparent = (f[at + 1] & 1) ? f[at + 4] : -1;
            }
            if (label == 0xFF && std::string(f.begin() + static_cast<std::ptrdiff_t>(at + 1), f.begin() + static_cast<std::ptrdiff_t>(at + 12)) == "NETSCAPE2.0") {
                gif.loops = f[at + 12] == 3 && f[at + 13] == 1 && read16(f, at + 14) == 0; // loop count 0: forever
            }
            while (f[at] != 0) at += f[at] + 1u; // skip the sub-blocks
            ++at;
            continue;
        }
        REQUIRE(block == 0x2C);
        const int left = read16(f, at), top = read16(f, at + 2), w = read16(f, at + 4), h = read16(f, at + 6);
        REQUIRE((f[at + 8] & 0x87) == 0x87); // a 256-entry local table
        at += 9;
        const std::size_t table = at;
        at += 3 * 256;
        REQUIRE(f[at++] == 8);
        std::vector<std::uint8_t> data;
        while (f[at] != 0) {
            data.insert(data.end(), f.begin() + static_cast<std::ptrdiff_t>(at + 1), f.begin() + static_cast<std::ptrdiff_t>(at + 1 + f[at]));
            at += f[at] + 1u;
        }
        ++at;
        const std::vector<std::uint8_t> indices = lzwDecode(data);
        REQUIRE(indices.size() == static_cast<std::size_t>(w * h));
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                const int index = indices[static_cast<std::size_t>(y * w + x)];
                if (index == transparent) continue;
                std::uint8_t* p = canvas.data() + 4 * static_cast<std::size_t>((top + y) * gif.width + left + x);
                for (int c = 0; c < 3; ++c) p[c] = f[table + 3 * static_cast<std::size_t>(index) + static_cast<std::size_t>(c)];
                p[3] = 255;
            }
        }
        gif.rects.push_back({left, top, w, h});
        gif.frames.push_back(canvas);
    }
    return gif;
}

std::vector<std::uint8_t> image(int width, int height, auto color) {
    std::vector<std::uint8_t> rgba;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const std::array<std::uint8_t, 3> c = color(x, y);
            rgba.insert(rgba.end(), {c[0], c[1], c[2], 255});
        }
    }
    return rgba;
}

std::filesystem::path tempFile(const char* name) { return std::filesystem::temp_directory_path() / name; }

} // namespace

TEST_CASE("GIF LZW round trip, through table resets and every code width") {
    std::mt19937 random(7);
    std::vector<std::vector<std::uint8_t>> inputs = {{}, {42}, std::vector<std::uint8_t>(70'000, 255)};
    for (const int alphabet : {2, 16, 256}) {
        std::uniform_int_distribution<int> pick(0, alphabet - 1);
        std::vector<std::uint8_t> data(50'000);
        for (std::uint8_t& v : data) v = static_cast<std::uint8_t>(pick(random));
        inputs.push_back(data);
    }
    // Every length up to a few thousand codes, so the input ends at each code-width boundary.
    std::uniform_int_distribution<int> byte(0, 255);
    std::vector<std::uint8_t> prefix(3'000);
    for (std::uint8_t& v : prefix) v = static_cast<std::uint8_t>(byte(random));
    for (std::size_t n = 1; n <= prefix.size(); n += 7) inputs.emplace_back(prefix.begin(), prefix.begin() + static_cast<std::ptrdiff_t>(n));

    for (const std::vector<std::uint8_t>& input : inputs) {
        CAPTURE(input.size());
        CHECK(lzwDecode(gifLzwEncode(input)) == input);
    }
}

TEST_CASE("animated GIF: exact colours, changed rectangles only, exact total delay, looping") {
    constexpr int w = 24;
    constexpr int h = 10;
    // 240 distinct colours: fewer than 255, so the palette holds them exactly.
    const auto gradient = [](int x, int y) {
        return std::array<std::uint8_t, 3>{static_cast<std::uint8_t>(10 * x), static_cast<std::uint8_t>(25 * y), 77};
    };
    const std::vector<std::uint8_t> first = image(w, h, gradient);
    const std::vector<std::uint8_t> second = image(w, h, [&](int x, int y) {
        if (x >= 5 && x < 8 && y >= 2 && y < 4) return std::array<std::uint8_t, 3>{255, 0, 255};
        return gradient(x, y);
    });

    const std::filesystem::path path = tempFile("pfl_test_animation.gif");
    {
        GifWriter writer;
        std::string error;
        REQUIRE(writer.open(path, w, h, 30, error));
        REQUIRE(writer.addFrame(first.data(), error));
        REQUIRE(writer.addFrame(second.data(), error));
        REQUIRE(writer.addFrame(second.data(), error)); // nothing changed
        REQUIRE(writer.close(error));
        CHECK(writer.frameCount() == 3);
    }
    const DecodedGif gif = decodeGif(path);
    std::filesystem::remove(path);

    CHECK(gif.width == w);
    CHECK(gif.height == h);
    CHECK(gif.loops);
    REQUIRE(gif.frames.size() == 3);
    CHECK(gif.delays == std::vector<int>{3, 4, 3}); // 0.10 s for 3 frames at 30 fps
    CHECK(gif.rects[0] == std::array<int, 4>{0, 0, w, h});
    CHECK(gif.rects[1] == std::array<int, 4>{5, 2, 3, 2});
    CHECK(gif.rects[2] == std::array<int, 4>{0, 0, 1, 1});
    CHECK(gif.frames[0] == first);
    CHECK(gif.frames[1] == second);
    CHECK(gif.frames[2] == second);
}

TEST_CASE("GIF with more than 255 colours stays close to the original") {
    constexpr int w = 64;
    constexpr int h = 64;
    std::mt19937 random(3);
    std::uniform_int_distribution<int> byte(0, 255);
    const std::vector<std::uint8_t> noise = image(w, h, [&](int, int) {
        return std::array<std::uint8_t, 3>{static_cast<std::uint8_t>(byte(random)), static_cast<std::uint8_t>(byte(random)),
                                           static_cast<std::uint8_t>(byte(random))};
    });
    const std::filesystem::path path = tempFile("pfl_test_noise.gif");
    {
        GifWriter writer;
        std::string error;
        REQUIRE(writer.open(path, w, h, 20, error));
        REQUIRE(writer.addFrame(noise.data(), error));
        REQUIRE(writer.close(error));
    }
    const DecodedGif gif = decodeGif(path);
    std::filesystem::remove(path);
    REQUIRE(gif.frames.size() == 1);
    double error = 0.0;
    for (std::size_t i = 0; i < noise.size(); ++i) {
        if (i % 4 != 3) error += std::abs(gif.frames[0][i] - noise[i]);
    }
    CHECK(error / (3.0 * w * h) < 16.0); // mean error per channel
}
