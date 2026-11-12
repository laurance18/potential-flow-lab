// help_window.hpp: the Help menu's windows.
//
// User guide (F1): docs/user_guide.md, embedded in the exe at build time and rendered in a window
// with a chapter list, a search box and the Markdown turned into ImGui text (headings, paragraphs,
// bullets, tables, code blocks, bold and inline code). Sign conventions: the complex potential of
// every element, what the signs of Λ, Γ, μ and the angles mean, the wall images and Cp. Keyboard
// shortcuts: every key and mouse gesture the program understands. About: version, author and the
// libraries and OpenGL renderer in use.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pfl {

// What the About window shows: the version, and rows such as ("OpenGL", "4.6 on ..."), ("GLFW", "3.5.1").
struct AboutInfo {
    std::string version;
    std::vector<std::pair<std::string, std::string>> system;
};

class HelpWindow {
public:
    HelpWindow();
    ~HelpWindow(); // defined where Chapter is complete
    HelpWindow(const HelpWindow&) = delete;
    HelpWindow& operator=(const HelpWindow&) = delete;

    void setAboutInfo(AboutInfo info) { m_about = std::move(info); }

    // Each opens its window, or brings it to the front if it is open already.
    void openGuide();
    void openConventions();
    void openShortcuts();
    void openAbout();

    // Draws whichever windows are open. Call every frame.
    void draw();

    // For tests and tools: the chapter titles of the embedded user guide (empty if it is not built in).
    std::vector<std::string> chapterTitles() const;

private:
    struct Block;
    struct Chapter;

    void parseGuide(std::string_view markdown);
    void drawGuide();
    void drawChapter(const Chapter& chapter);
    void drawBlock(const Block& block);
    void drawConventions();
    void drawShortcuts();
    void drawAbout();

    std::vector<Chapter> m_chapters;
    std::string m_guideTitle;
    std::size_t m_chapter = 0;
    char m_search[128] = {};
    bool m_scrollToTop = false;

    bool m_guideOpen = false;
    bool m_conventionsOpen = false;
    bool m_shortcutsOpen = false;
    bool m_aboutOpen = false;
    bool m_focusGuide = false;
    bool m_focusConventions = false;
    bool m_focusShortcuts = false;
    bool m_focusAbout = false;
    AboutInfo m_about;
};

} // namespace pfl
