#include "app/mp4_writer_win32.hpp"

#include <windows.h>

#include <codecapi.h>
#include <mfapi.h>
#include <mferror.h>
#include <mfidl.h>
#include <mfreadwrite.h>
#include <wrl/client.h>

#include <algorithm>
#include <format>

namespace pfl {
namespace {

using Microsoft::WRL::ComPtr;

constexpr LONGLONG k_ticksPerSecond = 10'000'000; // Media Foundation time unit: 100 ns
constexpr double k_bitsPerPixel = 0.25;           // per frame; particles and dye need a generous rate
constexpr double k_minBitrate = 2e6;
constexpr double k_maxBitrate = 50e6;

std::string failure(const char* what, HRESULT hr) {
    return std::format("Could not write the MP4 file: {} failed (error 0x{:08X}).", what, static_cast<unsigned long>(hr));
}

// Frame size, rate and square pixels: shared by the encoder's input and output types.
HRESULT setVideoFormat(IMFMediaType* type, const GUID& subtype, int width, int height, int fps) {
    HRESULT hr = type->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
    if (SUCCEEDED(hr)) hr = type->SetGUID(MF_MT_SUBTYPE, subtype);
    if (SUCCEEDED(hr)) hr = type->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);
    if (SUCCEEDED(hr)) hr = MFSetAttributeSize(type, MF_MT_FRAME_SIZE, static_cast<UINT32>(width), static_cast<UINT32>(height));
    if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(type, MF_MT_FRAME_RATE, static_cast<UINT32>(fps), 1);
    if (SUCCEEDED(hr)) hr = MFSetAttributeRatio(type, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
    return hr;
}

} // namespace

struct Mp4Writer::Impl {
    ComPtr<IMFSinkWriter> writer;
    DWORD stream = 0;
    int width = 0;
    int height = 0;
    int fps = 0;
    LONGLONG frames = 0;
    bool comInitialized = false;
    bool mfStarted = false;

    void shutdown() {
        writer.Reset();
        if (mfStarted) MFShutdown();
        if (comInitialized) CoUninitialize();
        mfStarted = comInitialized = false;
    }
};

Mp4Writer::Mp4Writer() = default;

Mp4Writer::~Mp4Writer() {
    std::string ignored;
    close(ignored);
}

bool Mp4Writer::open(const std::filesystem::path& path, int width, int height, int fps, std::string& error) {
    std::string ignored;
    close(ignored);
    if (width < 2 || height < 2 || width % 2 != 0 || height % 2 != 0 || fps < 1) {
        error = "Invalid MP4 size or frame rate.";
        return false;
    }
    // The DLLs are delay-loaded: check that they exist before the first call into them.
    if (LoadLibraryExW(L"mfplat.dll", nullptr, LOAD_LIBRARY_SEARCH_SYSTEM32) == nullptr ||
        LoadLibraryExW(L"mfreadwrite.dll", nullptr, LOAD_LIBRARY_SEARCH_SYSTEM32) == nullptr) {
        error = "MP4 export needs Windows Media Foundation, which this edition of Windows does not include "
                "(on Windows N editions, install the Media Feature Pack). Record an animated GIF instead.";
        return false;
    }

    auto impl = std::make_unique<Impl>();
    impl->width = width;
    impl->height = height;
    impl->fps = fps;
    impl->comInitialized = SUCCEEDED(CoInitializeEx(nullptr, COINIT_MULTITHREADED));
    HRESULT hr = MFStartup(MF_VERSION);
    if (FAILED(hr)) {
        error = failure("starting Media Foundation", hr);
        impl->shutdown();
        return false;
    }
    impl->mfStarted = true;

    const double bitrate = std::clamp(k_bitsPerPixel * width * height * fps, k_minBitrate, k_maxBitrate);
    ComPtr<IMFAttributes> attributes;
    ComPtr<IMFMediaType> output;
    ComPtr<IMFMediaType> input;
    const char* step = "creating the attributes";
    hr = MFCreateAttributes(&attributes, 1);
    if (SUCCEEDED(hr)) hr = attributes->SetGUID(MF_TRANSCODE_CONTAINERTYPE, MFTranscodeContainerType_MPEG4);
    if (SUCCEEDED(hr)) {
        step = "creating the file";
        hr = MFCreateSinkWriterFromURL(path.c_str(), nullptr, attributes.Get(), &impl->writer);
    }
    if (SUCCEEDED(hr)) {
        step = "setting up the H.264 stream";
        hr = MFCreateMediaType(&output);
    }
    if (SUCCEEDED(hr)) hr = setVideoFormat(output.Get(), MFVideoFormat_H264, width, height, fps);
    if (SUCCEEDED(hr)) hr = output->SetUINT32(MF_MT_AVG_BITRATE, static_cast<UINT32>(bitrate));
    if (SUCCEEDED(hr)) hr = output->SetUINT32(MF_MT_MPEG2_PROFILE, eAVEncH264VProfile_High);
    if (SUCCEEDED(hr)) hr = impl->writer->AddStream(output.Get(), &impl->stream);
    if (SUCCEEDED(hr)) {
        step = "setting up the H.264 encoder";
        hr = MFCreateMediaType(&input);
    }
    if (SUCCEEDED(hr)) hr = setVideoFormat(input.Get(), MFVideoFormat_RGB32, width, height, fps);
    // A positive stride: rows top to bottom (uncompressed RGB is otherwise taken as bottom-up).
    if (SUCCEEDED(hr)) hr = input->SetUINT32(MF_MT_DEFAULT_STRIDE, static_cast<UINT32>(4 * width));
    if (SUCCEEDED(hr)) hr = impl->writer->SetInputMediaType(impl->stream, input.Get(), nullptr);
    if (SUCCEEDED(hr)) {
        step = "starting the encoder";
        hr = impl->writer->BeginWriting();
    }
    if (FAILED(hr)) {
        error = failure(step, hr);
        impl->shutdown();
        return false;
    }
    m_impl = std::move(impl);
    return true;
}

bool Mp4Writer::addFrame(const std::uint8_t* rgba, std::string& error) {
    if (!m_impl) {
        error = "The MP4 file is not open.";
        return false;
    }
    Impl& impl = *m_impl;
    const DWORD size = static_cast<DWORD>(4 * impl.width * impl.height);
    ComPtr<IMFMediaBuffer> buffer;
    HRESULT hr = MFCreateMemoryBuffer(size, &buffer);
    BYTE* data = nullptr;
    if (SUCCEEDED(hr)) hr = buffer->Lock(&data, nullptr, nullptr);
    if (SUCCEEDED(hr)) {
        // MFVideoFormat_RGB32 is B, G, R, X in memory.
        for (DWORD i = 0; i < size; i += 4) {
            data[i + 0] = rgba[i + 2];
            data[i + 1] = rgba[i + 1];
            data[i + 2] = rgba[i + 0];
            data[i + 3] = 255;
        }
        hr = buffer->Unlock();
    }
    if (SUCCEEDED(hr)) hr = buffer->SetCurrentLength(size);

    ComPtr<IMFSample> sample;
    if (SUCCEEDED(hr)) hr = MFCreateSample(&sample);
    if (SUCCEEDED(hr)) hr = sample->AddBuffer(buffer.Get());
    // Exact frame times: frame k starts at k / fps seconds.
    const LONGLONG start = impl.frames * k_ticksPerSecond / impl.fps;
    const LONGLONG end = (impl.frames + 1) * k_ticksPerSecond / impl.fps;
    if (SUCCEEDED(hr)) hr = sample->SetSampleTime(start);
    if (SUCCEEDED(hr)) hr = sample->SetSampleDuration(end - start);
    if (SUCCEEDED(hr)) hr = impl.writer->WriteSample(impl.stream, sample.Get());
    if (FAILED(hr)) {
        error = failure("encoding a frame", hr);
        return false;
    }
    ++impl.frames;
    return true;
}

bool Mp4Writer::close(std::string& error) {
    if (!m_impl) return true;
    const HRESULT hr = m_impl->frames > 0 ? m_impl->writer->Finalize() : S_OK;
    m_impl->shutdown();
    m_impl.reset();
    if (FAILED(hr)) {
        error = failure("finishing the file", hr);
        return false;
    }
    return true;
}

} // namespace pfl
