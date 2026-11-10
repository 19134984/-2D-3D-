/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef IO_TEXT_TABLE_H
#define IO_TEXT_TABLE_H

#include <string>
#include <vector>
#include <sstream>
#include <ostream>
#include <optional>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <utility>
#include <type_traits>

namespace olb {

/// Consistent house-style text rendering for console / introspection output
/**
 * A tiny formatting toolbox shared by mesh, statistics and field-storage
 * summaries so that they no longer hand-roll setw/`|`/`---+---`. Everything
 * renders through a plain std::ostream so that passing an OstreamManager keeps
 * the `[Tag]` prefix intact.
 **/
namespace io {

/// Reusable text-formatting primitives (tables, key/value lists, human units)
namespace text {

/// Table rendering layout
enum struct TableLayout {
  columns,   ///< bare aligned columns
  ruled,     ///< columns with a rule under the header
  grid,      ///< inner column separators and header rule (no outer box)
  framed,    ///< outer box, no inner separators
  markdown   ///< GitHub-flavored pipe table
};

/// Parse a table layout name
inline std::optional<TableLayout> parseTableLayout(const std::string& name)
{
  if (name == "columns")  return TableLayout::columns;
  if (name == "ruled")    return TableLayout::ruled;
  if (name == "grid")     return TableLayout::grid;
  if (name == "framed")   return TableLayout::framed;
  if (name == "markdown") return TableLayout::markdown;
  return std::nullopt;
}

/// Process-global text styling, consulted by every io::text formatter
struct Style {
  /// Status-line grammar: `key=value`, `key value`, or `key=value;`-separated
  enum struct Fields { KeyValue, Labeled, Semicolon };
  Fields      fields  = Fields::Semicolon;
  /// Table layout, overridable per Table
  TableLayout table   = TableLayout::ruled;
  /// Box-drawing glyphs (rules, separators) vs plain ASCII
  bool        unicode = true;
  /// ANSI color escapes
  bool        color   = false;
};

/// The process-global style instance
inline Style& style() {
  static Style instance;
  return instance;
}

/// Scientific notation with the given number of decimals, e.g. 0.00014 -> "1.4e-04"
inline std::string sci(double value, int precision=1)
{
  std::ostringstream s;
  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(precision);
  s << value;
  return s.str();
}

/// Fixed-point with the given number of decimals, e.g. 2.0 -> "2.000"
inline std::string fixed(double value, int precision=3)
{
  std::ostringstream s;
  s.setf(std::ios::fixed, std::ios::floatfield);
  s.precision(precision);
  s << value;
  return s.str();
}

/// Repeat a (possibly multi-byte) unit string n times
inline std::string repeat(const std::string& unit, std::size_t n)
{
  std::string out;
  out.reserve(unit.size() * n);
  for (std::size_t i=0; i < n; ++i) {
    out += unit;
  }
  return out;
}

/// A titled section rule, e.g. "━━━ Mesh ━━━━━━━━━━━━━━━━━━"
inline std::string section(const std::string& title)
{
  const std::string dash = style().unicode ? "\xe2\x94\x81" : "="; // ━
  constexpr std::size_t width = 52;
  const std::size_t used = 4 + title.size() + 1;
  const std::size_t tail = width > used ? width - used : 3;
  std::string line = repeat(dash, 3) + " " + title + " " + repeat(dash, tail);
  return style().color ? "\x1b[1m" + line + "\x1b[0m" : line;
}

/// Format a byte count in binary units, e.g. 737388 -> "720.1 KiB"
inline std::string humanBytes(std::size_t bytes)
{
  constexpr const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
  if (bytes < 1024) {
    return std::to_string(bytes) + " B";
  }
  double value = static_cast<double>(bytes);
  int unit = 0;
  while (value >= 1024.0 && unit < 5) {
    value /= 1024.0;
    ++unit;
  }
  std::ostringstream s;
  s.setf(std::ios::fixed);
  s.precision(1);
  s << value << " " << units[unit];
  return s.str();
}

/// Format a fraction as a percentage, e.g. 0.847 -> "85%" (default 0 decimals)
inline std::string percent(double frac, int precision=0)
{
  std::ostringstream s;
  s.setf(std::ios::fixed);
  s.precision(precision);
  s << (frac * 100.0) << "%";
  return s.str();
}

/// Format a large count with thousands separators, e.g. 1234567 -> "1,234,567"
inline std::string humanCount(std::size_t count)
{
  std::string digits = std::to_string(count);
  std::string out;
  int c = 0;
  for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
    if (c != 0 && c % 3 == 0) {
      out.push_back(',');
    }
    out.push_back(*it);
    ++c;
  }
  return std::string(out.rbegin(), out.rend());
}

/// Compact SI-style count for dense one-liners, e.g. 12345 -> "12.3k", 1.7e12 -> "1.7T"
inline std::string shortCount(std::size_t count)
{
  constexpr const char* units[] = {"", "k", "M", "G", "T", "P", "E"};
  if (count < 1000) {
    return std::to_string(count);
  }
  double value = static_cast<double>(count);
  int unit = 0;
  while (value >= 1000.0 && unit < 6) {
    value /= 1000.0;
    ++unit;
  }
  std::ostringstream s;
  s.setf(std::ios::fixed);
  s.precision(1);
  s << value << units[unit];
  return s.str();
}

/// SI-scaled magnitude of a rate, e.g. lups=1.36e7 -> "13.6 MLUPs", 2.0e10 -> "20.0 GLUPs"
inline std::string lups(double lupsPerSecond)
{
  constexpr const char* units[] = {"", "k", "M", "G", "T", "P", "E"};
  double value = lupsPerSecond;
  int unit = 0;
  while (value >= 1000.0 && unit < 6) {
    value /= 1000.0;
    ++unit;
  }
  std::ostringstream s;
  s.setf(std::ios::fixed);
  s.precision(1);
  s << value << " " << units[unit] << "LUPs";
  return s.str();
}

/// Detail level of the case-wide statistics report
enum struct Detail {
  glance,   ///< one condensed line over all components (cfl / energy health-first)
  compact,  ///< one classic release line per lattice (step=…; maxCFL=…; avgEnergy=…)
  summary,  ///< one table row per component, levels folded
  full      ///< per-component dashboard with per-level sub-table and min/max/mean
};

/// Parse a detail name
inline std::optional<Detail> parseDetail(const std::string& name)
{
  if (name == "glance")  return Detail::glance;
  if (name == "compact") return Detail::compact;
  if (name == "summary") return Detail::summary;
  if (name == "full")    return Detail::full;
  return std::nullopt;
}

/// Name of a detail level
inline std::string detailName(Detail detail)
{
  switch (detail) {
    case Detail::glance:  return "glance";
    case Detail::compact: return "compact";
    case Detail::summary: return "summary";
    case Detail::full:    return "full";
  }
  return "compact";
}

/// Make a (possibly namespace-qualified) name safe as a filename component,
/// e.g. "names::NavierStokes" -> "names_NavierStokes"
inline std::string filenameSafe(std::string name)
{
  for (auto pos = name.find("::"); pos != std::string::npos; pos = name.find("::")) {
    name.replace(pos, 2, "_");
  }
  return name;
}

/// Convert an arbitrary streamable value into a string
template <typename T>
inline std::string toCell(const T& value)
{
  if constexpr (std::is_convertible_v<T, std::string>) {
    return std::string(value);
  } else {
    std::ostringstream s;
    s << value;
    return s.str();
  }
}

/// Auto-width, column-aligned table with optional horizontal rule rows
class Table {
public:
  enum struct Align { Left, Right };

  struct Column {
    std::string title;
    Align       align    = Align::Left;
    /// Cells longer than this are truncated with an ellipsis (0 = no limit)
    std::size_t maxWidth = 0;
  };

private:
  std::vector<Column>                   _columns;
  std::vector<std::vector<std::string>> _rows;
  std::vector<bool>                     _isRule; // parallel to _rows
  std::optional<TableLayout>            _layout;

public:
  Table(std::initializer_list<Column> columns)
    : _columns(columns)
  { }

  /// Override the process-global table layout for this table
  Table& layout(TableLayout chosen)
  {
    _layout = chosen;
    return *this;
  }

  /// Append a data row; each cell is stringified via operator<<
  template <typename... CELLS>
  Table& row(CELLS&&... cells)
  {
    std::vector<std::string> cellStrings{toCell(std::forward<CELLS>(cells))...};
    for (std::size_t i=0; i < cellStrings.size() && i < _columns.size(); ++i) {
      cellStrings[i] = clip(cellStrings[i], _columns[i].maxWidth);
    }
    _rows.push_back(std::move(cellStrings));
    _isRule.push_back(false);
    return *this;
  }

  /// Append a full-width horizontal rule row
  Table& rule()
  {
    _rows.push_back({});
    _isRule.push_back(true);
    return *this;
  }

  /// Render the table to out, indented by the given prefix
  void render(std::ostream& out, const std::string& indent="  ") const
  {
    std::vector<std::size_t> width(_columns.size(), 0);
    for (std::size_t i=0; i < _columns.size(); ++i) {
      width[i] = _columns[i].title.size();
    }
    for (const auto& row : _rows) {
      for (std::size_t i=0; i < row.size() && i < width.size(); ++i) {
        width[i] = std::max(width[i], row[i].size());
      }
    }

    const TableLayout layout = _layout.value_or(style().table);
    const bool color = style().color;
    const bool uni   = style().unicode;
    const std::string dash = uni ? "\xe2\x94\x80" : "-"; // ─

    std::string colGap      = "  ";
    std::size_t colGapWidth = 2;
    std::string rowLead;
    std::string rowTail;
    switch (layout) {
      case TableLayout::grid:
        colGap = uni ? " \xe2\x94\x82 " : " | "; // │
        colGapWidth = 3;
        break;
      case TableLayout::markdown:
        colGap = " | "; colGapWidth = 3;
        rowLead = "| "; rowTail = " |";
        break;
      case TableLayout::framed:
        rowLead = uni ? "\xe2\x94\x82 " : "| "; // │
        rowTail = uni ? " \xe2\x94\x82" : " |";
        break;
      default:
        break;
    }

    std::size_t innerWidth = 0;
    for (std::size_t i=0; i < width.size(); ++i) {
      innerWidth += width[i] + (i+1 < width.size() ? colGapWidth : 0);
    }

    auto emitCells = [&](const std::vector<std::string>& cells, bool header) {
      out << indent << rowLead;
      for (std::size_t i=0; i < _columns.size(); ++i) {
        const std::string cell = pad(i < cells.size() ? cells[i] : std::string{},
                                     width[i], _columns[i].align);
        if (header && color) {
          out << "\x1b[1m" << cell << "\x1b[0m";
        } else {
          out << cell;
        }
        if (i+1 < _columns.size()) {
          out << colGap;
        }
      }
      out << rowTail << std::endl;
    };

    auto emitRule = [&]() {
      out << indent;
      switch (layout) {
        case TableLayout::grid: {
          const std::string cross = uni ? "\xe2\x94\x80\xe2\x94\xbc\xe2\x94\x80" : "-+-"; // ─┼─
          for (std::size_t i=0; i < width.size(); ++i) {
            out << repeat(dash, width[i]);
            if (i+1 < width.size()) {
              out << cross;
            }
          }
          break;
        }
        case TableLayout::framed: {
          const std::string left  = uni ? "\xe2\x94\x9c" : "+"; // ├
          const std::string right = uni ? "\xe2\x94\xa4" : "+"; // ┤
          out << left << repeat(dash, innerWidth + 2) << right;
          break;
        }
        case TableLayout::markdown: {
          out << "|";
          for (std::size_t i=0; i < width.size(); ++i) {
            std::string cell(std::max<std::size_t>(3, width[i]), '-');
            if (_columns[i].align == Align::Right) {
              cell.back() = ':';
            }
            out << " " << cell << " |";
          }
          break;
        }
        default: {
          out << repeat(dash, innerWidth);
          break;
        }
      }
      out << std::endl;
    };

    if (layout == TableLayout::framed) {
      const std::string tl = uni ? "\xe2\x94\x8c" : "+"; // ┌
      const std::string tr = uni ? "\xe2\x94\x90" : "+"; // ┐
      out << indent << tl << repeat(dash, innerWidth + 2) << tr << std::endl;
    }

    std::vector<std::string> titles(_columns.size());
    for (std::size_t i=0; i < _columns.size(); ++i) {
      titles[i] = _columns[i].title;
    }
    emitCells(titles, true);

    if (layout != TableLayout::columns) {
      emitRule();
    }

    for (std::size_t r=0; r < _rows.size(); ++r) {
      if (_isRule[r]) {
        emitRule();
      } else {
        emitCells(_rows[r], false);
      }
    }

    if (layout == TableLayout::framed) {
      const std::string bl = uni ? "\xe2\x94\x94" : "+"; // └
      const std::string br = uni ? "\xe2\x94\x98" : "+"; // ┘
      out << indent << bl << repeat(dash, innerWidth + 2) << br << std::endl;
    }
  }

private:
  /// Truncate value to maxWidth code units, appending an ellipsis (0 = no limit)
  static std::string clip(const std::string& value, std::size_t maxWidth)
  {
    if (maxWidth == 0 || value.size() <= maxWidth) {
      return value;
    }
    if (maxWidth <= 3) {
      return value.substr(0, maxWidth);
    }
    return value.substr(0, maxWidth - 3) + "...";
  }

  static std::string pad(const std::string& value, std::size_t width, Align align)
  {
    if (value.size() >= width) {
      return value;
    }
    const std::string fill(width - value.size(), ' ');
    return align == Align::Right ? fill + value : value + fill;
  }
};

/// Aligned "key : value" block
class KeyValueList {
private:
  std::vector<std::pair<std::string,std::string>> _items;

public:
  KeyValueList& item(std::string key, std::string value)
  {
    _items.emplace_back(std::move(key), std::move(value));
    return *this;
  }

  template <typename T>
  KeyValueList& item(std::string key, const T& value)
  {
    _items.emplace_back(std::move(key), toCell(value));
    return *this;
  }

  bool empty() const {
    return _items.empty();
  }

  void render(std::ostream& out, const std::string& indent="  ") const
  {
    std::size_t width = 0;
    for (const auto& [key, value] : _items) {
      width = std::max(width, key.size());
    }
    for (const auto& [key, value] : _items) {
      out << indent << key << std::string(width - key.size(), ' ')
          << " : " << value << std::endl;
    }
  }
};

/// Bullet separator honouring the unicode style knob
inline const char* bullet()
{
  return style().unicode ? "\xc2\xb7" : "|"; // ·
}

/// Compact one-line status of fields, rendered as `key=value` or `key value` per style
class StatusLine {
private:
  struct Field { std::string key; std::string value; bool isToken; bool isSeparator; };
  std::vector<Field> _fields;

public:
  /// A standalone data token (e.g. "done", a component name, an MLUPs value)
  StatusLine& token(std::string text)
  {
    _fields.push_back({std::move(text), {}, true, false});
    return *this;
  }

  /// A structural divider (e.g. "|" or "·") — never semicolon-wrapped
  StatusLine& separator(std::string text)
  {
    _fields.push_back({std::move(text), {}, true, true});
    return *this;
  }

  StatusLine& field(std::string key, std::string value)
  {
    _fields.push_back({std::move(key), std::move(value), false, false});
    return *this;
  }

  template <typename T>
  StatusLine& field(std::string key, const T& value)
  {
    return field(std::move(key), toCell(value));
  }

  std::string str() const
  {
    const auto grammar = style().fields;
    const bool kv        = grammar != Style::Fields::Labeled;
    const bool semicolon = grammar == Style::Fields::Semicolon;
    const bool color     = style().color;
    std::ostringstream out;
    for (std::size_t i=0; i < _fields.size(); ++i) {
      if (i != 0) {
        const bool sepAdjacent = _fields[i].isSeparator || _fields[i-1].isSeparator;
        out << (semicolon && !sepAdjacent ? "; " : "  ");
      }
      const auto& f = _fields[i];
      if (f.isToken) {
        out << f.key;
        continue;
      }
      if (color) {
        out << "\x1b[2m" << f.key << (kv ? "=" : " ") << "\x1b[0m" << f.value;
      } else {
        out << f.key << (kv ? "=" : " ") << f.value;
      }
    }
    return out.str();
  }

  void render(std::ostream& out) const
  {
    out << str() << std::endl;
  }
};

}

}

}

#endif
