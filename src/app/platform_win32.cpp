#include "app/platform_win32.hpp"

#include <windows.h>

#include <commdlg.h>
#include <shellapi.h>
#include <wincodec.h>
#include <wrl/client.h>

#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cwchar>
#include <format>
#include <string>
#include <vector>

namespace pfl {
namespace {

constexpr wchar_t k_sceneFilter[] = L"Flow scenes (*.json)\0*.json\0All files (*.*)\0*.*\0";
constexpr wchar_t k_pngFilter[] = L"PNG images (*.png)\0*.png\0All files (*.*)\0*.*\0";
constexpr wchar_t k_csvFilter[] = L"CSV files (*.csv)\0*.csv\0All files (*.*)\0*.*\0";
constexpr wchar_t k_mp4Filter[] = L"MP4 videos (*.mp4)\0*.mp4\0All files (*.*)\0*.*\0";
constexpr wchar_t k_gifFilter[] = L"Animated GIFs (*.gif)\0*.gif\0All files (*.*)\0*.*\0";
constexpr const char* k_usage =
    "Usage: PotentialFlowLab.exe [scene.json] [--export-image out.png] [--export-field out.csv] "
    "[--export-lines out.csv] [--export-surface out.csv] [--export-video out.mp4|out.gif "
    "[--video-seconds 6] [--video-fps 30]]";

using Microsoft::WRL::ComPtr;

std::wstring widen(std::string_view utf8) {
    if (utf8.empty()) return {};
    const int length = MultiByteToWideChar(CP_UTF8, 0, utf8.data(), static_cast<int>(utf8.size()), nullptr, 0);
    std::wstring wide(static_cast<std::size_t>(length), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8.data(), static_cast<int>(utf8.size()), wide.data(), length);
    return wide;
}

HWND ownerWindow(GLFWwindow* window) { return window != nullptr ? glfwGetWin32Window(window) : nullptr; }

// Runs GetOpenFileNameW / GetSaveFileNameW with `filter`. `initial` pre-fills the name.
std::optional<std::filesystem::path> fileDialog(GLFWwindow* owner, bool save, const std::filesystem::path& initial,
                                                const wchar_t* filter = k_sceneFilter,
                                                const wchar_t* extension = L"json") {
    std::array<wchar_t, 4096> buffer{};
    initial.native().copy(buffer.data(), buffer.size() - 1);

    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = ownerWindow(owner);
    ofn.lpstrFilter = filter;
    ofn.lpstrFile = buffer.data();
    ofn.nMaxFile = static_cast<DWORD>(buffer.size());
    ofn.lpstrDefExt = extension;
    ofn.Flags = OFN_NOCHANGEDIR | OFN_PATHMUSTEXIST |
                (save ? OFN_OVERWRITEPROMPT : static_cast<DWORD>(OFN_FILEMUSTEXIST));
    const BOOL ok = save ? GetSaveFileNameW(&ofn) : GetOpenFileNameW(&ofn);
    if (ok == FALSE) return std::nullopt;
    return std::filesystem::path(buffer.data());
}

} // namespace

void logMessage(std::string_view message) {
    std::string line(message);
    line += '\n';
    std::fputs(line.c_str(), stderr);
    std::fflush(stderr);
    OutputDebugStringW(widen(line).c_str());
}

void showFatalError(std::string_view title, std::string_view message) {
    logMessage(std::string(title) + ": " + std::string(message));
    MessageBoxW(nullptr, widen(message).c_str(), widen(title).c_str(), MB_OK | MB_ICONERROR | MB_SETFOREGROUND);
}

void showErrorMessage(GLFWwindow* owner, std::string_view title, std::string_view message) {
    logMessage(std::string(title) + ": " + std::string(message));
    MessageBoxW(ownerWindow(owner), widen(message).c_str(), widen(title).c_str(), MB_OK | MB_ICONWARNING);
}

CommandLine parseCommandLine() {
    // CommandLineToArgvW gives Unicode arguments, unlike main()'s argv.
    CommandLine result;
    int count = 0;
    LPWSTR* args = CommandLineToArgvW(GetCommandLineW(), &count);
    if (args == nullptr) return result;
    for (int i = 1; i < count && result.error.empty(); ++i) {
        const std::wstring arg = args[i];
        const std::string name = std::filesystem::path(arg).string();
        std::optional<std::filesystem::path>* target = nullptr;
        if (arg == L"--export-image") target = &result.exportImage;
        if (arg == L"--export-field") target = &result.exportField;
        if (arg == L"--export-lines") target = &result.exportLines;
        if (arg == L"--export-surface") target = &result.exportSurface;
        if (arg == L"--export-video") target = &result.exportVideo;
        const bool seconds = arg == L"--video-seconds";
        const bool fps = arg == L"--video-fps";

        if (target != nullptr || seconds || fps) {
            if (i + 1 >= count) {
                result.error = name + (target != nullptr ? " needs a file name." : " needs a number.");
                continue;
            }
            const std::wstring value = args[++i];
            if (target != nullptr) {
                *target = std::filesystem::path(value);
                continue;
            }
            wchar_t* end = nullptr;
            const double number = std::wcstod(value.c_str(), &end);
            const bool whole = end != nullptr && *end == 0 && !value.empty();
            if (seconds && whole && number >= 0.1 && number <= 3600.0) {
                result.videoSeconds = number;
            } else if (fps && whole && number >= 1.0 && number <= 120.0 && number == std::floor(number)) {
                result.videoFps = static_cast<int>(number);
            } else {
                result.error = name + (seconds ? " needs a number of seconds from 0.1 to 3600."
                                               : " needs a whole number of frames per second from 1 to 120.");
            }
        } else if (arg.starts_with(L"--")) {
            result.error = "Unknown option " + name + ".\n\n" + k_usage;
        } else if (!result.scene) {
            result.scene = std::filesystem::path(arg);
        }
    }
    LocalFree(args);
    if (result.error.empty() && !result.exportVideo && (result.videoSeconds || result.videoFps)) {
        result.error = std::string("--video-seconds and --video-fps go with --export-video.\n\n") + k_usage;
    }
    return result;
}

std::optional<std::filesystem::path> showOpenSceneDialog(GLFWwindow* owner) { return fileDialog(owner, false, {}); }

std::optional<std::filesystem::path> showSaveSceneDialog(GLFWwindow* owner, const std::filesystem::path& suggested) {
    return fileDialog(owner, true, suggested);
}

std::optional<std::filesystem::path> showSavePngDialog(GLFWwindow* owner, const std::filesystem::path& suggested) {
    return fileDialog(owner, true, suggested, k_pngFilter, L"png");
}

std::optional<std::filesystem::path> showSaveCsvDialog(GLFWwindow* owner, const std::filesystem::path& suggested) {
    return fileDialog(owner, true, suggested, k_csvFilter, L"csv");
}

std::optional<std::filesystem::path> showSaveVideoDialog(GLFWwindow* owner, const std::filesystem::path& suggested,
                                                         bool gif) {
    return fileDialog(owner, true, suggested, gif ? k_gifFilter : k_mp4Filter, gif ? L"gif" : L"mp4");
}

bool writePng(const std::filesystem::path& path, int width, int height, const std::uint8_t* rgba, std::string& error) {
    // COM may already be initialised on this thread (e.g. by a file dialog); balance only our own call.
    const HRESULT init = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
    const bool uninitialize = SUCCEEDED(init);

    // WIC's PNG encoder takes 24-bit BGR rows.
    const UINT stride = static_cast<UINT>(width) * 3;
    const std::size_t pixels = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    std::vector<BYTE> bgr(3 * pixels);
    for (std::size_t i = 0; i < pixels; ++i) {
        bgr[3 * i + 0] = rgba[4 * i + 2];
        bgr[3 * i + 1] = rgba[4 * i + 1];
        bgr[3 * i + 2] = rgba[4 * i + 0];
    }

    HRESULT hr = S_OK;
    {
        ComPtr<IWICImagingFactory> factory;
        ComPtr<IWICStream> stream;
        ComPtr<IWICBitmapEncoder> encoder;
        ComPtr<IWICBitmapFrameEncode> frame;
        WICPixelFormatGUID format = GUID_WICPixelFormat24bppBGR;
        hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&factory));
        if (SUCCEEDED(hr)) hr = factory->CreateStream(&stream);
        if (SUCCEEDED(hr)) hr = stream->InitializeFromFilename(path.c_str(), GENERIC_WRITE);
        if (SUCCEEDED(hr)) hr = factory->CreateEncoder(GUID_ContainerFormatPng, nullptr, &encoder);
        if (SUCCEEDED(hr)) hr = encoder->Initialize(stream.Get(), WICBitmapEncoderNoCache);
        if (SUCCEEDED(hr)) hr = encoder->CreateNewFrame(&frame, nullptr);
        if (SUCCEEDED(hr)) hr = frame->Initialize(nullptr);
        if (SUCCEEDED(hr)) hr = frame->SetSize(static_cast<UINT>(width), static_cast<UINT>(height));
        if (SUCCEEDED(hr)) hr = frame->SetPixelFormat(&format);
        if (SUCCEEDED(hr) && format != GUID_WICPixelFormat24bppBGR) hr = E_FAIL; // encoder wants another layout
        if (SUCCEEDED(hr)) hr = frame->WritePixels(static_cast<UINT>(height), stride, static_cast<UINT>(bgr.size()), bgr.data());
        if (SUCCEEDED(hr)) hr = frame->Commit();
        if (SUCCEEDED(hr)) hr = encoder->Commit();
    } // release the COM objects before CoUninitialize
    if (uninitialize) CoUninitialize();

    if (FAILED(hr)) {
        error = std::format("Could not write the PNG file (error 0x{:08X}). Is the folder writable, or is the file "
                            "open in another program?",
                            static_cast<unsigned long>(hr));
        return false;
    }
    return true;
}

} // namespace pfl
