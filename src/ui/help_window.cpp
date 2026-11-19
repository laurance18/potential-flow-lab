#define IMGUI_DEFINE_MATH_OPERATORS
#include "ui/help_window.hpp"

#include "embedded_resources.hpp" // pfl::docs (src/ui/CMakeLists.txt embeds docs/user_guide.md)
#include "ui/style.hpp"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <format>
#include <iterator>
#include <optional>
#include <span>

namespace pfl {

// ------------------------------------------------------------------------------------------------
// The user guide's Markdown, reduced to what docs/user_guide.md uses.

struct HelpWindow::Block {
    enum class Kind { Heading, Paragraph, Bullet, Table, Code };
    Kind kind = Kind::Paragraph;
    int level = 0;                              // Bullet: indent level (0 = top)
    std::string text;                           // Heading, Paragraph, Bullet: inline Markdown
    std::vector<std::vector<std::string>> rows; // Table: the header row first
    std::vector<std::string> lines;             // Code
};

struct HelpWindow::Chapter {
    std::string title;
    std::string searchText; // everything in the chapter, lower case, for the search box
    std::vector<Block> blocks;
};

namespace {

constexpr const char* k_guideTitle = "User guide";
constexpr const char* k_conventionsTitle = "Sign conventions";
constexpr const char* k_shortcutsTitle = "Keyboard shortcuts";
constexpr const char* k_aboutTitle = "About Potential Flow Lab";

constexpr ImVec4 k_bodyText{0.88f, 0.88f, 0.90f, 1.0f};
constexpr ImVec4 k_boldText{1.0f, 1.0f, 1.0f, 1.0f};
constexpr ImVec4 k_codeText{1.0f, 0.82f, 0.56f, 1.0f};
constexpr ImVec4 k_headingText{0.55f, 0.76f, 1.0f, 1.0f};

std::string_view trimView(std::string_view s) {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t'))
        s.remove_prefix(1);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r'))
        s.remove_suffix(1);
    return s;
}

std::string lowerAscii(std::string_view s) {
    std::string out(s);
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

// "| a | b |" -> {"a", "b"}; "\|" is a literal bar.
std::vector<std::string> splitCells(std::string_view row) {
    std::vector<std::string> cells;
    std::string cell;
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (row[i] == '\\' && i + 1 < row.size() && row[i + 1] == '|') {
            cell += '|';
            ++i;
        } else if (row[i] == '|') {
            cells.emplace_back(trimView(cell));
            cell.clear();
        } else {
            cell += row[i];
        }
    }
    cells.emplace_back(trimView(cell));
    // The row starts and ends with a bar, so the first and last pieces are empty.
    if (!cells.empty() && cells.front().empty()) cells.erase(cells.begin());
    if (!cells.empty() && cells.back().empty()) cells.pop_back();
    return cells;
}

bool isTableSeparator(std::string_view row) {
    bool dash = false;
    for (const char c : row) {
        if (c == '-')
            dash = true;
        else if (c != '|' && c != ':' && c != ' ')
            return false;
    }
    return dash;
}

// Inline Markdown: **bold**, `code`, *italic* (shown plain) and [text](url) (shown as text), split
// into words so that they can be wrapped.
struct Fragment {
    std::string text;
    bool bold = false;
    bool code = false;
};
using Word = std::vector<Fragment>;

std::vector<Word> splitWords(std::string_view md) {
    std::vector<Word> words;
    Word word;
    Fragment fragment;
    bool bold = false;
    bool code = false;
    const auto flushFragment = [&] {
        if (!fragment.text.empty()) word.push_back(fragment);
        fragment = Fragment{.text = {}, .bold = bold, .code = code};
    };
    const auto flushWord = [&] {
        flushFragment();
        if (!word.empty()) words.push_back(word);
        word.clear();
    };
    for (std::size_t i = 0; i < md.size(); ++i) {
        const char c = md[i];
        if (c == '`') {
            flushFragment();
            code = !code;
            fragment.code = code;
            continue;
        }
        if (!code && c == '*') {
            if (i + 1 < md.size() && md[i + 1] == '*') {
                flushFragment();
                bold = !bold;
                fragment.bold = bold;
                ++i;
            }
            continue; // a single '*' marks italics, which the font cannot show
        }
        if (!code && c == '[') {
            const std::size_t close = md.find("](", i);
            const std::size_t end = close == std::string_view::npos ? close : md.find(')', close);
            if (end != std::string_view::npos) {
                for (const char t : md.substr(i + 1, close - i - 1)) {
                    if (t == ' ')
                        flushWord();
                    else
                        fragment.text += t;
                }
                i = end;
                continue;
            }
        }
        if (c == ' ' || c == '\n' || c == '\t')
            flushWord();
        else
            fragment.text += c;
    }
    flushWord();
    return words;
}

std::string plainText(std::string_view md) {
    std::string out;
    for (const Word& word : splitWords(md)) {
        for (const Fragment& f : word)
            out += f.text;
        out += ' ';
    }
    return out;
}

// Draws one fragment at the cursor.
void drawFragment(const Fragment& f) {
    const float s = uiScale();
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    if (f.code) {
        const ImVec2 size = ImGui::CalcTextSize(f.text.c_str());
        ImGui::GetWindowDrawList()->AddRectFilled(pos - ImVec2(2.0f * s, 1.0f * s),
                                                  pos + size + ImVec2(2.0f * s, 1.0f * s), IM_COL32(255, 255, 255, 20),
                                                  3.0f * s);
        ImGui::PushStyleColor(ImGuiCol_Text, k_codeText);
        ImGui::TextUnformatted(f.text.c_str());
        ImGui::PopStyleColor();
        return;
    }
    if (f.bold) {
        ImGui::PushStyleColor(ImGuiCol_Text, k_boldText);
        ImGui::TextUnformatted(f.text.c_str());
        // Drawn twice, half a pixel apart: the nearest thing to bold with one font.
        ImGui::GetWindowDrawList()->AddText(nullptr, 0.0f, pos + ImVec2(0.5f * s, 0.0f),
                                            ImGui::GetColorU32(ImGuiCol_Text), f.text.c_str());
        ImGui::PopStyleColor();
        return;
    }
    ImGui::TextUnformatted(f.text.c_str());
}

// Draws inline Markdown from the cursor, wrapping at the screen x `right`; wrapped lines start
// under the first one (so bullet text stays clear of its bullet).
void drawInline(std::string_view md, float right, bool allBold = false) {
    const std::vector<Word> words = splitWords(md);
    if (words.empty()) return;
    const float startX = ImGui::GetCursorPosX();
    const float space = ImGui::CalcTextSize(" ").x;
    bool first = true;
    for (Word word : words) {
        float width = 0.0f;
        for (Fragment& f : word) {
            if (allBold) f.bold = true;
            width += ImGui::CalcTextSize(f.text.c_str()).x;
        }
        if (!first) {
            ImGui::SameLine(0.0f, space);
            if (ImGui::GetCursorScreenPos().x + width > right) {
                ImGui::NewLine();
                ImGui::SetCursorPosX(startX);
            }
        }
        for (std::size_t i = 0; i < word.size(); ++i) {
            if (i > 0) ImGui::SameLine(0.0f, 0.0f);
            drawFragment(word[i]);
        }
        first = false;
    }
}

// Right edge of the text area, from a cursor at the start of a line.
float lineRight() {
    return ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
}

void blockGap() {
    ImGui::Dummy(ImVec2(0.0f, 0.35f * ImGui::GetFontSize()));
}

// A title in a larger font.
void drawTitle(const char* text, float scale) {
    ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * scale);
    ImGui::TextUnformatted(text);
    ImGui::PopFont();
}

// A two-column table of (name, description) rows, the names in bold.
void drawPairTable(const char* id, std::span<const std::pair<const char*, const char*>> rows, float nameWeight = 1.4f,
                   float textWeight = 4.0f) {
    if (!ImGui::BeginTable(id, 2,
                           ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        return;
    }
    ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, nameWeight);
    ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, textWeight);
    for (const auto& [name, text] : rows) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        drawInline(name, lineRight(), true);
        ImGui::TableSetColumnIndex(1);
        drawInline(text, lineRight());
    }
    ImGui::EndTable();
}

// Begins a floating window centred in the work area, sized for its content on first use.
bool beginHelpWindow(const char* title, bool& open, bool& focus, ImVec2 size, ImVec2 minSize) {
    const float s = uiScale();
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImVec2 work = viewport->WorkSize;
    ImGui::SetNextWindowSize(ImVec2(std::min(size.x * s, 0.92f * work.x), std::min(size.y * s, work.y - 24.0f * s)),
                             ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(viewport->WorkPos + work * 0.5f, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSizeConstraints(minSize * s, ImVec2(FLT_MAX, FLT_MAX));
    if (focus) {
        ImGui::SetNextWindowFocus();
        focus = false;
    }
    const bool visible = ImGui::Begin(title, &open, ImGuiWindowFlags_NoCollapse);
    // Esc closes the window, unless a text field is using it.
    if (visible && ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) && !ImGui::IsAnyItemActive() &&
        ImGui::IsKeyPressed(ImGuiKey_Escape)) {
        open = false;
    }
    return visible;
}

} // namespace

// ------------------------------------------------------------------------------------------------

HelpWindow::HelpWindow() {
    if (const auto text = docs::findText("docs/user_guide.md")) parseGuide(*text);
}

HelpWindow::~HelpWindow() = default;

void HelpWindow::openGuide() {
    m_guideOpen = true;
    m_focusGuide = true;
}
void HelpWindow::openConventions() {
    m_conventionsOpen = true;
    m_focusConventions = true;
}
void HelpWindow::openShortcuts() {
    m_shortcutsOpen = true;
    m_focusShortcuts = true;
}
void HelpWindow::openAbout() {
    m_aboutOpen = true;
    m_focusAbout = true;
}

std::vector<std::string> HelpWindow::chapterTitles() const {
    std::vector<std::string> titles;
    for (const Chapter& c : m_chapters)
        titles.push_back(c.title);
    return titles;
}

void HelpWindow::parseGuide(std::string_view markdown) {
    m_chapters.clear();
    Chapter current;
    std::vector<Block> intro;        // text before the first "## " heading goes into the first chapter
    std::optional<std::size_t> open; // the paragraph, bullet or table being continued
    bool inCode = false;
    Block code;

    const auto endChapter = [&] {
        if (current.title.empty()) {
            intro = std::move(current.blocks);
        } else {
            std::string search = lowerAscii(current.title);
            for (const Block& b : current.blocks) {
                search += ' ';
                search += lowerAscii(plainText(b.text));
                for (const auto& row : b.rows) {
                    for (const std::string& cell : row) {
                        search += ' ';
                        search += lowerAscii(plainText(cell));
                    }
                }
                for (const std::string& line : b.lines) {
                    search += ' ';
                    search += lowerAscii(line);
                }
            }
            current.searchText = std::move(search);
            m_chapters.push_back(std::move(current));
        }
        current = Chapter{};
        open.reset();
    };

    std::size_t pos = 0;
    while (pos <= markdown.size()) {
        const std::size_t end = markdown.find('\n', pos);
        const std::string_view raw =
            markdown.substr(pos, end == std::string_view::npos ? std::string_view::npos : end - pos);
        pos = end == std::string_view::npos ? markdown.size() + 1 : end + 1;

        if (inCode) {
            if (trimView(raw).starts_with("```")) {
                inCode = false;
                current.blocks.push_back(std::move(code));
            } else {
                std::string_view line = raw;
                if (line.ends_with('\r')) line.remove_suffix(1);
                code.lines.emplace_back(line);
            }
            continue;
        }
        const std::string_view line = trimView(raw);
        const std::size_t firstText = raw.find_first_not_of(" \t");
        const std::size_t indent = firstText == std::string_view::npos ? 0 : firstText;
        if (line.empty()) {
            open.reset();
            continue;
        }
        if (line.starts_with("```")) {
            inCode = true;
            code = Block{.kind = Block::Kind::Code};
            open.reset();
            continue;
        }
        if (line.starts_with("# ")) {
            m_guideTitle = std::string(trimView(line.substr(2)));
            open.reset();
            continue;
        }
        if (line.starts_with("## ")) {
            endChapter();
            current.title = std::string(trimView(line.substr(3)));
            current.blocks = std::move(intro);
            intro.clear();
            continue;
        }
        if (line.starts_with("### ")) {
            current.blocks.push_back(
                Block{.kind = Block::Kind::Heading, .text = std::string(trimView(line.substr(4)))});
            open.reset();
            continue;
        }
        if (line.starts_with('|')) {
            if (isTableSeparator(line)) continue;
            if (!open || current.blocks[*open].kind != Block::Kind::Table) {
                current.blocks.push_back(Block{.kind = Block::Kind::Table});
                open = current.blocks.size() - 1;
            }
            current.blocks[*open].rows.push_back(splitCells(line));
            continue;
        }
        if (line.starts_with("- ") || line.starts_with("* ")) {
            current.blocks.push_back(Block{.kind = Block::Kind::Bullet,
                                           .level = static_cast<int>(indent / 2),
                                           .text = std::string(trimView(line.substr(2)))});
            open = current.blocks.size() - 1;
            continue;
        }
        if (open && (current.blocks[*open].kind == Block::Kind::Paragraph ||
                     current.blocks[*open].kind == Block::Kind::Bullet)) {
            current.blocks[*open].text += ' ';
            current.blocks[*open].text += line;
            continue;
        }
        current.blocks.push_back(Block{.kind = Block::Kind::Paragraph, .text = std::string(line)});
        open = current.blocks.size() - 1;
    }
    if (inCode) current.blocks.push_back(std::move(code));
    endChapter();
    if (m_chapters.empty() && !intro.empty()) {
        current.title = m_guideTitle.empty() ? "User guide" : m_guideTitle;
        current.blocks = std::move(intro);
        endChapter();
    }
}

void HelpWindow::draw() {
    if (m_guideOpen) drawGuide();
    if (m_conventionsOpen) drawConventions();
    if (m_shortcutsOpen) drawShortcuts();
    if (m_aboutOpen) drawAbout();
}

// ------------------------------------------------------------------------------------------------
// User guide

void HelpWindow::drawGuide() {
    const float s = uiScale();
    if (!beginHelpWindow(k_guideTitle, m_guideOpen, m_focusGuide, ImVec2(920.0f, 700.0f), ImVec2(560.0f, 380.0f))) {
        ImGui::End();
        return;
    }
    if (m_chapters.empty()) {
        ImGui::TextWrapped(
            "The user guide is not built into this exe. It is docs/user_guide.md in the source repository.");
        ImGui::End();
        return;
    }

    // Left: the chapters, filtered by the search box.
    const std::string query = lowerAscii(trimView(m_search));
    std::vector<bool> matches(m_chapters.size(), true);
    bool any = false;
    for (std::size_t i = 0; i < m_chapters.size(); ++i) {
        matches[i] = query.empty() || m_chapters[i].searchText.find(query) != std::string::npos;
        any = any || matches[i];
    }
    if (!query.empty() && any && !matches[m_chapter]) {
        m_chapter = static_cast<std::size_t>(std::find(matches.begin(), matches.end(), true) - matches.begin());
        m_scrollToTop = true;
    }

    ImGui::BeginChild("chapters", ImVec2(250.0f * s, 0.0f), ImGuiChildFlags_Borders);
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputTextWithHint("##search", "Search...", m_search, sizeof(m_search));
    ImGui::SetItemTooltip("Show only the chapters that mention this.");
    ImGui::Spacing();
    for (std::size_t i = 0; i < m_chapters.size(); ++i) {
        if (!matches[i]) continue;
        if (ImGui::Selectable(m_chapters[i].title.c_str(), i == m_chapter) && i != m_chapter) {
            m_chapter = i;
            m_scrollToTop = true;
        }
    }
    if (!any) ImGui::TextDisabled("No chapter mentions that.");
    ImGui::EndChild();

    ImGui::SameLine();
    ImGui::BeginChild("chapter", ImVec2(0.0f, 0.0f), ImGuiChildFlags_AlwaysUseWindowPadding);
    if (m_scrollToTop) {
        ImGui::SetScrollY(0.0f);
        m_scrollToTop = false;
    }
    ImGui::PushStyleColor(ImGuiCol_Text, k_bodyText);
    drawChapter(m_chapters[m_chapter]);
    if (query.empty() && m_chapter + 1 < m_chapters.size()) {
        blockGap();
        ImGui::Separator();
        const std::string next = std::format("Next: {} \xE2\x80\xBA", m_chapters[m_chapter + 1].title);
        if (ImGui::SmallButton(next.c_str())) {
            ++m_chapter;
            m_scrollToTop = true;
        }
    }
    ImGui::PopStyleColor();
    ImGui::EndChild();
    ImGui::End();
}

void HelpWindow::drawChapter(const Chapter& chapter) {
    ImGui::PushStyleColor(ImGuiCol_Text, k_boldText);
    drawTitle(chapter.title.c_str(), 1.45f);
    ImGui::PopStyleColor();
    blockGap();
    for (std::size_t i = 0; i < chapter.blocks.size(); ++i) {
        ImGui::PushID(static_cast<int>(i));
        drawBlock(chapter.blocks[i]);
        ImGui::PopID();
    }
}

void HelpWindow::drawBlock(const Block& block) {
    const float s = uiScale();
    switch (block.kind) {
    case Block::Kind::Heading:
        blockGap();
        ImGui::PushStyleColor(ImGuiCol_Text, k_headingText);
        drawTitle(block.text.c_str(), 1.15f);
        ImGui::PopStyleColor();
        ImGui::Spacing();
        break;
    case Block::Kind::Paragraph:
        drawInline(block.text, lineRight());
        blockGap();
        break;
    case Block::Kind::Bullet: {
        const float right = lineRight();
        const float indent = block.level * 18.0f * s; // Indent(0) would apply ImGui's default indent
        if (indent > 0.0f) ImGui::Indent(indent);
        ImGui::Bullet();
        ImGui::SameLine();
        drawInline(block.text, right);
        if (indent > 0.0f) ImGui::Unindent(indent);
        ImGui::Spacing();
        break;
    }
    case Block::Kind::Table: {
        std::size_t columns = 0;
        for (const auto& row : block.rows)
            columns = std::max(columns, row.size());
        if (columns == 0) break;
        std::vector<float> weights(columns, 1.0f);
        for (const auto& row : block.rows) {
            for (std::size_t c = 0; c < row.size(); ++c) {
                weights[c] = std::max(weights[c],
                                      std::clamp(std::sqrt(static_cast<float>(plainText(row[c]).size())), 1.5f, 6.0f));
            }
        }
        constexpr ImGuiTableFlags flags =
            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp;
        if (ImGui::BeginTable("table", static_cast<int>(columns), flags)) {
            for (std::size_t c = 0; c < columns; ++c)
                ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, weights[c]);
            for (std::size_t r = 0; r < block.rows.size(); ++r) {
                ImGui::TableNextRow(r == 0 ? ImGuiTableRowFlags_Headers : ImGuiTableRowFlags_None);
                for (std::size_t c = 0; c < block.rows[r].size() && c < columns; ++c) {
                    ImGui::TableSetColumnIndex(static_cast<int>(c));
                    drawInline(block.rows[r][c], lineRight(), r == 0);
                }
            }
            ImGui::EndTable();
        }
        blockGap();
        break;
    }
    case Block::Kind::Code:
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.0f, 0.0f, 0.0f, 0.35f));
        ImGui::PushStyleColor(ImGuiCol_Text, k_codeText);
        if (ImGui::BeginChild("code", ImVec2(0.0f, 0.0f),
                              ImGuiChildFlags_AutoResizeY | ImGuiChildFlags_Borders |
                                  ImGuiChildFlags_AlwaysUseWindowPadding,
                              ImGuiWindowFlags_HorizontalScrollbar)) {
            for (const std::string& line : block.lines)
                ImGui::TextUnformatted(line.c_str());
        }
        ImGui::EndChild();
        ImGui::PopStyleColor(2);
        blockGap();
        break;
    }
}

// ------------------------------------------------------------------------------------------------
// Sign conventions

void HelpWindow::drawConventions() {
    if (!beginHelpWindow(k_conventionsTitle, m_conventionsOpen, m_focusConventions, ImVec2(820.0f, 640.0f),
                         ImVec2(520.0f, 320.0f))) {
        ImGui::End();
        return;
    }
    ImGui::PushStyleColor(ImGuiCol_Text, k_bodyText);

    ImGui::PushStyleColor(ImGuiCol_Text, k_headingText);
    drawTitle("Coordinates and angles", 1.15f);
    ImGui::PopStyleColor();
    ImGui::Spacing();
    static constexpr std::pair<const char*, const char*> k_coordinates[] = {
        {"Axes", "x to the right, y up, both to the same scale. Lengths, speeds and strengths are in any consistent "
                 "units of your choice; the program carries no units of its own."},
        {"Angles", "In degrees, counter-clockwise from the +x axis: the uniform-flow direction α, the doublet axis β, "
                   "the sector walls θ₁ and θ₂ and the wall angle φ."},
        {"Complex notation", "z = x + iy, F(z) = φ + iψ (velocity potential and stream function), W = dF/dz = u − iv. "
                             "The velocity at a point is the sum of every element's W, with the wall images added."},
    };
    drawPairTable("coordinates", k_coordinates);
    blockGap();

    ImGui::PushStyleColor(ImGuiCol_Text, k_headingText);
    drawTitle("Elements (ζ = z − z₀, from the element's position)", 1.15f);
    ImGui::PopStyleColor();
    ImGui::Spacing();
    constexpr ImGuiTableFlags flags =
        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp;
    if (ImGui::BeginTable("elements", 4, flags)) {
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, 1.3f);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, 3.6f);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, 1.6f);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, 1.8f);
        static constexpr const char* k_rows[][4] = {
            {"Element", "Parameters and their signs", "F(z)", "W = u − iv"},
            {"Uniform flow", "Speed V∞ ≥ 0 and direction α: u = V∞ cos α, v = V∞ sin α.", "V∞ e^(−iα) z", "V∞ e^(−iα)"},
            {"Source / sink",
             "Strength Λ = volume flow rate per unit depth. Λ > 0 is a source (flow outwards), Λ < 0 a sink. "
             "Radial velocity V_r = Λ/(2πr).",
             "(Λ/2π) ln ζ", "Λ / (2π ζ)"},
            {"Point vortex",
             "Circulation Γ. **Γ > 0 turns clockwise** (Anderson's convention), Γ < 0 counter-clockwise. "
             "Tangential velocity V_θ = −Γ/(2πr).",
             "(iΓ/2π) ln ζ", "iΓ / (2π ζ)"},
            {"Doublet",
             "Strength μ = κ/2π and axis angle β. With a uniform flow at α, β = α gives the cylinder of radius "
             "a = √(μ/V∞). A negative μ is the same doublet turned by 180°.",
             "μ e^(iβ) / ζ", "−μ e^(iβ) / ζ²"},
            {"Sector flow",
             "Strength V, vertex z₀, walls θ₁ < θ₂; the opening β_s = θ₂ − θ₁ ≤ 360° sets n = 180°/β_s. "
             "The flow exists only inside the wedge. With V > 0 the fluid comes in along the second wall (θ₂) "
             "and leaves along the first (θ₁); V < 0 reverses it.",
             "V (ζ e^(−iθ₁))ⁿ", "V n e^(−iθ₁) (ζ e^(−iθ₁))ⁿ⁻¹"},
            {"Symmetry wall",
             "A point on the wall and its angle φ. The wall is enforced with images: a source keeps Λ, a vortex "
             "gets −Γ, a doublet gets the axis 2φ − β, and a uniform flow keeps only its component along the "
             "wall. Sector flows are not mirrored. Walls at 180°/k to each other and parallel walls (a channel) "
             "are exact; other layouts are truncated and the status bar says so.",
             "F + F_image", "W + W_image"},
            {"Source-panel body",
             "An ellipse (a circle when a = b), rectangle, NACA 4-digit airfoil or drawn / imported outline, half-sizes "
             "a, b along and across its own axis, centre z₀ and rotation θ (an airfoil's angle of attack is −θ), cut "
             "into N straight panels numbered clockwise. Panel j from z₁ to z₂ in direction e^(iθ_j) is a source sheet of "
             "strength λ_j (> 0 pushes fluid out); the λ_j are solved so that v·n = 0 at every panel midpoint.",
             "Σ (λ_j/2π) ∫ ln(z − z′) ds′ along panel j", "Σ (λ_j/2π) e^(−iθ_j) ln((z − z₁)/(z − z₂))"},
            {"Vortex-sheet airfoil",
             "A camber line (flat plate, NACA 4-digit mean line or circular arc) of chord c, placed by the middle of its "
             "chord z₀ and turned by θ (in a stream along +x, α = −θ). It carries a vortex sheet γ(s), circulation per "
             "unit length, **γ > 0 clockwise** like Γ: the flow over the top is faster, u_upper − u_lower = γ. "
             "Γ = ∫γ ds > 0 gives lift L′ = ρV∞Γ. N straight panels from the leading to the trailing edge, γ linear "
             "along each; the Kutta condition sets γ = 0 at the trailing edge. φ jumps by Γ across the ray running on "
             "from the trailing edge along the chord.",
             "(i/2π) ∫ γ(s) ln(z − z(s)) ds", "(i/2π) ∫ γ(s) ds / (z − z(s))"},
        };
        for (std::size_t r = 0; r < std::size(k_rows); ++r) {
            ImGui::TableNextRow(r == 0 ? ImGuiTableRowFlags_Headers : ImGuiTableRowFlags_None);
            for (int c = 0; c < 4; ++c) {
                ImGui::TableSetColumnIndex(c);
                drawInline(k_rows[r][c], lineRight(), r == 0 || c == 0);
            }
        }
        ImGui::EndTable();
    }
    blockGap();

    ImGui::PushStyleColor(ImGuiCol_Text, k_headingText);
    drawTitle("Pressure and forces", 1.15f);
    ImGui::PopStyleColor();
    ImGui::Spacing();
    static constexpr std::pair<const char*, const char*> k_pressure[] = {
        {"Pressure coefficient",
         "Cp = 1 − |V|²/V∞² (Bernoulli), with V∞ the speed of the uniform flow (the sum of the "
         "uniform flows in the scene). Without a uniform flow there is no reference speed and no Cp. "
         "Cp = 1 at a stagnation point; Cp < 0 where the flow is faster than the free stream."},
        {"Stagnation points",
         "Where W = 0 (Cp = 1). The dividing streamlines through them bound the bodies: the fluid of a "
         "source or doublet that the free stream never enters."},
        {"Lift and drag",
         "Per unit span, from integrating −Cp n around a body surface: D′ along the free stream, L′ at "
         "90° counter-clockwise from it. Shown as D′/q∞ and L′/q∞ with q∞ = ½ρV∞², and as c_d and c_l "
         "with the chord c as the reference length."},
        {"Exact values",
         "For a body alone in a uniform stream: D′ = 0 (d'Alembert) and L′ = ρ V∞ Γ with Γ > 0 clockwise "
         "(Kutta–Joukowski), so L′/q∞ = 2Γ/V∞. A clockwise vortex on a cylinder in a flow along +x "
         "lifts it upwards."},
        {"Thin airfoils",
         "ΔCp = Cp_lower − Cp_upper = 2γ/V∞, > 0 where the airfoil is pushed up. c_l = 2Γ/(V∞c) (Kutta–Joukowski), "
         "c_m,c/4 > 0 nose up. Thin-airfoil theory: γ(θ) = 2V∞[A₀(1 + cos θ)/sin θ + Σ Aₙ sin nθ] with "
         "x = (c/2)(1 − cos θ), c_l = 2π(A₀ + A₁/2), c_m,c/4 = (π/4)(A₂ − A₁). The angle of attack α is the "
         "stream's direction relative to the chord, positive nose up."},
        {"φ and ψ",
         "φ is multivalued around a vortex (by Γ) and ψ around a source (by Λ). The exports use the principal "
         "branch, which jumps across the line running in −x from the element; around a thin airfoil φ jumps "
         "across the ray from its trailing edge along the chord."},
    };
    drawPairTable("pressure", k_pressure);

    ImGui::PopStyleColor();
    ImGui::End();
}

// ------------------------------------------------------------------------------------------------
// Keyboard shortcuts

void HelpWindow::drawShortcuts() {
    if (!beginHelpWindow(k_shortcutsTitle, m_shortcutsOpen, m_focusShortcuts, ImVec2(640.0f, 640.0f),
                         ImVec2(440.0f, 300.0f))) {
        ImGui::End();
        return;
    }
    ImGui::PushStyleColor(ImGuiCol_Text, k_bodyText);

    struct Group {
        const char* title;
        std::span<const std::pair<const char*, const char*>> rows;
    };
    static constexpr std::pair<const char*, const char*> k_files[] = {
        {"Ctrl+N", "New scene"},
        {"Ctrl+O", "Open a scene file"},
        {"Ctrl+S", "Save"},
        {"Ctrl+Shift+S", "Save As"},
        {"Alt+F4", "Quit (asks about unsaved changes)"},
    };
    static constexpr std::pair<const char*, const char*> k_edit[] = {
        {"Ctrl+Z", "Undo the last change to the elements or dye rakes"},
        {"Ctrl+Y", "Redo (also Ctrl+Shift+Z)"},
    };
    static constexpr std::pair<const char*, const char*> k_windows[] = {
        {"Ctrl+P", "Parameters panel (live sliders and the sweep)"},
        {"Ctrl+Shift+F", "Performance stats (frame rate, GPU load) in the menu bar"},
        {"F5", "Physics checks"},
        {"F1", "This user guide; Esc closes a help window"},
    };
    static constexpr std::pair<const char*, const char*> k_plot[] = {
        {"Scroll", "Zoom at the mouse"},
        {"Drag", "Pan the plot (on empty space)"},
        {"Hover", "Element parameters, stagnation point coordinates, velocity and Cp in the status bar"},
        {"Double-click", "Edit the element under the mouse"},
        {"Drag an element", "Move it; the flow follows live"},
        {"Shift while dragging", "Snap the position to a grid"},
        {"Esc while dragging", "Put the element back"},
        {"Click a stagnation point", "Show or hide its dividing streamlines"},
        {"Click to place on, element tool, click", "Place a new element there (Shift snaps, Esc cancels)"},
    };
    static constexpr std::pair<const char*, const char*> k_rakes[] = {
        {"Rake tool, then drag", "Draw a dye rake (a click places a vertical one)"},
        {"Drag a rake", "Move it; drag an end handle to turn or stretch it"},
        {"Right-click a rake", "Number of emitters, or delete it"},
        {"Delete", "Remove the rake under the mouse"},
    };
    static constexpr std::pair<const char*, const char*> k_dialogs[] = {
        {"Esc", "Cancel a dialog"},
        {"Ctrl+click a slider", "Type an exact value (Parameters panel and dialogs)"},
        {"Right-click a slider", "Sweep this parameter (Parameters panel)"},
    };
    static constexpr std::pair<const char*, const char*> k_debug[] = {
        {"Ctrl+R", "Reload the shaders from shaders/ (Debug builds only)"},
    };
    static constexpr Group k_groups[] = {
        {"Files", k_files},       {"Edit", k_edit},       {"Panels and windows", k_windows},
        {"Plot", k_plot},         {"Dye rakes", k_rakes}, {"Dialogs and sliders", k_dialogs},
        {"Development", k_debug},
    };
    for (const Group& group : k_groups) {
        ImGui::PushStyleColor(ImGuiCol_Text, k_headingText);
        drawTitle(group.title, 1.15f);
        ImGui::PopStyleColor();
        ImGui::Spacing();
        drawPairTable(group.title, group.rows, 1.6f, 3.4f);
        blockGap();
    }

    ImGui::PopStyleColor();
    ImGui::End();
}

// ------------------------------------------------------------------------------------------------
// About

void HelpWindow::drawAbout() {
    if (!beginHelpWindow(k_aboutTitle, m_aboutOpen, m_focusAbout, ImVec2(560.0f, 480.0f), ImVec2(420.0f, 300.0f))) {
        ImGui::End();
        return;
    }
    ImGui::PushStyleColor(ImGuiCol_Text, k_bodyText);

    ImGui::PushStyleColor(ImGuiCol_Text, k_boldText);
    drawTitle("Potential Flow Lab", 1.6f);
    ImGui::PopStyleColor();
    const std::string version = m_about.version.empty() ? std::string("Version unknown") : "Version " + m_about.version;
    ImGui::TextUnformatted(version.c_str());
    blockGap();
    ImGui::TextWrapped(
        "An interactive two-dimensional potential-flow simulator: elementary flows (uniform flow, source "
        "and sink, vortex, doublet, sector flow) and symmetry walls superposed, with evenly spaced "
        "streamlines, equipotential lines, pressure contours, surface pressure and forces, tracer particles "
        "and dye, parameter sweeps, recordings and built-in physics checks.");
    blockGap();
    ImGui::TextWrapped("Senior mechanical-engineering semester project, Fall 2026. Written in C++20 by Deniz Işık.");
    blockGap();

    ImGui::PushStyleColor(ImGuiCol_Text, k_headingText);
    drawTitle("Built with", 1.15f);
    ImGui::PopStyleColor();
    ImGui::Spacing();
    std::vector<std::pair<const char*, const char*>> rows;
    rows.emplace_back("Dear ImGui", IMGUI_VERSION);
    rows.emplace_back("ImPlot", IMPLOT_VERSION);
    for (const auto& [name, value] : m_about.system)
        rows.emplace_back(name.c_str(), value.c_str());
    drawPairTable("about", rows, 1.4f, 3.6f);
    blockGap();

    if (ImGui::Button("Copy")) {
        std::string text = "Potential Flow Lab " + m_about.version + "\n";
        for (const auto& [name, value] : rows)
            text += std::format("{}: {}\n", name, value);
        ImGui::SetClipboardText(text.c_str());
    }
    ImGui::SetItemTooltip("Copy the version and the list above, for example to report a problem.");

    ImGui::PopStyleColor();
    ImGui::End();
}

} // namespace pfl
