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

#ifndef CASE_STATISTICS_REPORTER_H
#define CASE_STATISTICS_REPORTER_H

#include "reporter.h"
#include "refinement/lattice.h"
#include "io/textTable.h"
#include "io/input.h"
#include "io/ostreamManager.h"

#include <string>
#include <vector>
#include <algorithm>

namespace olb {

namespace reporter {

/// min/max/mean of one quantity
struct QuantityStat {
  std::string         name;
  unsigned            components = 1;
  std::vector<double> min;
  std::vector<double> max;
  std::vector<double> mean;
};

/// Per-level statistics snapshot
struct LevelStat {
  std::size_t cells     = 0;
  double      maxCFL    = 0;
  double      avgRho    = 0;
  double      avgEnergy = 0;
  /// rho, u and registered fields
  std::vector<QuantityStat> quantities;
};

/// Per-component snapshot
struct ComponentStat {
  std::string            name;
  bool                   refined = false;
  std::vector<LevelStat> levels;
};

namespace detail {

inline double worstCFL(const ComponentStat& component) {
  double cfl = 0;
  for (const auto& level : component.levels) {
    cfl = std::max(cfl, level.maxCFL);
  }
  return cfl;
}

inline double meanEnergy(const ComponentStat& component) {
  double sumEnergy = 0;
  std::size_t cells = 0;
  for (const auto& level : component.levels) {
    sumEnergy += level.avgEnergy * static_cast<double>(level.cells);
    cells     += level.cells;
  }
  return cells > 0 ? sumEnergy / static_cast<double>(cells) : 0.0;
}

inline io::text::StatusLine header(std::size_t iT, double physT) {
  return io::text::StatusLine().field("step", iT).field("t", io::text::fixed(physT) + "s");
}

inline void renderGlance(std::ostream& out, std::size_t iT, double physT,
                         const std::vector<ComponentStat>& components) {
  auto line = header(iT, physT);
  line.separator("|");
  for (std::size_t i=0; i < components.size(); ++i) {
    const auto& component = components[i];
    if (i != 0) {
      line.separator(io::text::bullet());
    }
    line.token(component.name)
        .field("cfl", io::text::sci(worstCFL(component)))
        .field("E", io::text::sci(meanEnergy(component)));
    if (component.refined) {
      line.field("lvl", component.levels.size());
    }
  }
  line.render(out);
}

inline void renderCompactLine(std::ostream& out, const std::string& label,
                              std::size_t iT, double physT, const LevelStat& level) {
  io::text::StatusLine line;
  if (!label.empty()) {
    line.token(label);
  }
  line.field("step", iT)
      .field("physT", io::text::fixed(physT) + "s")
      .field("maxCFL", io::text::sci(level.maxCFL))
      .field("avgEnergy", io::text::sci(level.avgEnergy))
      .field("avgRho", io::text::fixed(level.avgRho));
  for (const auto& q : level.quantities) {
    std::string value = "[";
    for (unsigned iC=0; iC < q.components; ++iC) {
      if (iC > 0) {
        value += ", ";
      }
      value += io::text::sci(q.min[iC]) + "/" + io::text::sci(q.max[iC])
             + "/" + io::text::sci(q.mean[iC]);
    }
    value += "]";
    line.field(q.name, value);
  }
  line.render(out);
}

inline void renderCompact(std::ostream& out, std::size_t iT, double physT,
                          const std::vector<ComponentStat>& components) {
  std::size_t lines = 0;
  for (const auto& component : components) {
    lines += component.levels.size();
  }
  const bool label = lines > 1;
  for (const auto& component : components) {
    for (std::size_t iL=0; iL < component.levels.size(); ++iL) {
      std::string name;
      if (label) {
        name = component.name;
        if (component.refined) {
          name += "/L" + std::to_string(iL);
        }
      }
      renderCompactLine(out, name, iT, physT, component.levels[iL]);
    }
  }
}

inline void renderSummary(std::ostream& out, std::size_t iT, double physT,
                          const std::vector<ComponentStat>& components) {
  header(iT, physT).render(out);
  io::text::Table table({
    {"component"}, {"cells", io::text::Table::Align::Right},
    {"avgRho", io::text::Table::Align::Right},
    {"maxCFL", io::text::Table::Align::Right},
    {"avgEnergy", io::text::Table::Align::Right},
    {"levels", io::text::Table::Align::Right},
  });
  for (const auto& component : components) {
    double sumRho = 0, sumEnergy = 0;
    std::size_t cells = 0;
    for (const auto& level : component.levels) {
      sumRho    += level.avgRho * static_cast<double>(level.cells);
      sumEnergy += level.avgEnergy * static_cast<double>(level.cells);
      cells     += level.cells;
    }
    const double invCells = cells > 0 ? 1.0 / static_cast<double>(cells) : 0.0;
    table.row(component.name, io::text::shortCount(cells),
              io::text::fixed(sumRho * invCells), io::text::sci(worstCFL(component)),
              io::text::sci(sumEnergy * invCells),
              component.refined ? std::to_string(component.levels.size()) : std::string{"-"});
  }
  table.render(out);
}

inline void renderFull(std::ostream& out, std::size_t iT, double physT,
                       const std::vector<ComponentStat>& components) {
  header(iT, physT).render(out);
  for (const auto& component : components) {
    out << std::endl << component.name;
    if (component.refined) {
      out << " (" << component.levels.size() << " levels)";
    }
    out << std::endl;

    io::text::Table levels({
      {"level"}, {"cells", io::text::Table::Align::Right},
      {"maxCFL", io::text::Table::Align::Right},
      {"avgRho", io::text::Table::Align::Right},
      {"avgEnergy", io::text::Table::Align::Right},
    });
    for (std::size_t iL=0; iL < component.levels.size(); ++iL) {
      const auto& level = component.levels[iL];
      levels.row(component.refined ? "L" + std::to_string(iL) : std::string{"-"},
                 io::text::shortCount(level.cells), io::text::sci(level.maxCFL),
                 io::text::fixed(level.avgRho), io::text::sci(level.avgEnergy));
    }
    levels.render(out);

    io::text::Table quantities({
      {"quantity"}, {"min", io::text::Table::Align::Right},
      {"max", io::text::Table::Align::Right},
      {"mean", io::text::Table::Align::Right},
    });
    const auto& q0 = component.levels.empty() ? std::vector<QuantityStat>{}
                                              : component.levels.front().quantities;
    for (const auto& q : q0) {
      for (unsigned iC=0; iC < q.components; ++iC) {
        const std::string label = q.components > 1 ? q.name + "[" + std::to_string(iC) + "]"
                                                   : q.name;
        quantities.row(label, io::text::sci(q.min[iC]), io::text::sci(q.max[iC]),
                       io::text::sci(q.mean[iC]));
      }
    }
    quantities.render(out);
  }
}

}

/// Render the statistics report
inline void render(std::ostream& out, std::size_t iT, double physT,
                   io::text::Detail detail, const std::vector<ComponentStat>& components) {
  if (components.empty()) {
    return;
  }
  switch (detail) {
    case io::text::Detail::glance:  detail::renderGlance(out, iT, physT, components);  break;
    case io::text::Detail::compact: detail::renderCompact(out, iT, physT, components); break;
    case io::text::Detail::summary: detail::renderSummary(out, iT, physT, components); break;
    case io::text::Detail::full:    detail::renderFull(out, iT, physT, components);    break;
  }
}

/// Case-wide statistics report
class Statistics {
private:
  /// The converter defining iT
  template <typename CASE>
  static decltype(auto) referenceConverter(CASE& myCase) {
    auto& component = myCase.get(typename CASE::map_t::keys_t::template get<0>{});
    if constexpr (CASE::is_refined) {
      return component.getLevel(0).getUnitConverter();
    } else {
      return component.getUnitConverter();
    }
  }

  static LevelStat gatherLevel(auto& lattice, bool wantQuantities) {
    auto& stats = lattice.getStatistics();
    LevelStat levelStat;
    levelStat.cells     = lattice.getCuboidDecomposition().getLatticeVolume();
    levelStat.maxCFL    = stats.getMaxU();
    levelStat.avgRho    = stats.getAverageRho();
    levelStat.avgEnergy = stats.getAverageEnergy();
    if (wantQuantities) {
      for (const auto& entry : stats.getEntries()) {
        QuantityStat quantity;
        quantity.name       = entry.name;
        quantity.components = entry.components;
        for (unsigned iC=0; iC < entry.components; ++iC) {
          quantity.min.push_back(entry.min[iC]);
          quantity.max.push_back(entry.max[iC]);
          quantity.mean.push_back(entry.mean[iC]);
        }
        levelStat.quantities.push_back(std::move(quantity));
      }
    }
    return levelStat;
  }

  template <typename CASE>
  static std::vector<ComponentStat> gather(CASE& myCase, std::size_t iT, io::text::Detail forDetail) {
    const bool wantQuantities = forDetail == io::text::Detail::full
                             || forDetail == io::text::Detail::compact;
    std::vector<ComponentStat> snapshots;
    CASE::map_t::keys_t::for_each([&](auto name) {
      using NAME = typename decltype(name)::type;
      using COMPONENT = typename CASE::map_t::template value<NAME>;
      auto& component = myCase.get(NAME{});
      ComponentStat componentStat;
      std::string componentName{fields::name<NAME>()};
      if (componentName.rfind("names::", 0) == 0) {
        componentName = componentName.substr(7);
      }
      componentStat.name = componentName;
      if constexpr (is_refined_lattice<COMPONENT>::value) {
        componentStat.refined = true;
        for (std::size_t iLevel=0; iLevel < component.getLevelCount(); ++iLevel) {
          auto& level = component.getLevel(iLevel);
          componentStat.levels.push_back(gatherLevel(level, wantQuantities));
          level.getStatistics().writeCSV(static_cast<int>(iT),
            level.getUnitConverter().getPhysTime(iT));
        }
        snapshots.push_back(std::move(componentStat));
      } else if constexpr (requires (COMPONENT& c) { c.getStatistics(); c.getUnitConverter(); }) {
        componentStat.levels.push_back(gatherLevel(component, wantQuantities));
        component.getStatistics().writeCSV(static_cast<int>(iT),
          component.getUnitConverter().getPhysTime(iT));
        snapshots.push_back(std::move(componentStat));
      }
    });
    return snapshots;
  }

  template <typename CASE>
  void emit(CASE& myCase, std::size_t iT, io::text::Detail forDetail) {
    auto snapshots = gather(myCase, iT, forDetail);
    if (snapshots.empty()) {
      return;
    }
    const double physT = referenceConverter(myCase).getPhysTime(iT);
    OstreamManager clout(std::cout, "Statistics");
    render(clout, iT, physT, forDetail, snapshots);
  }

  /// Lattice steps between reports (0 disables)
  template <typename CASE, typename PARAMETERS>
  static std::size_t interval(CASE& myCase, PARAMETERS& params) {
    const double physInterval = params.template get<olb::parameters::PHYS_ITER_T>();
    if (physInterval <= 0) {
      return 0;
    }
    return static_cast<std::size_t>(referenceConverter(myCase).getLatticeTime(physInterval));
  }

  /// Verbosity from the top-level --detail option (compact / summary / full)
  static io::text::Detail detail() {
    if (auto value = singleton::input().resolve(OverrideSource::Query{{}, "detail", true, false})) {
      return io::text::parseDetail(*value).value_or(io::text::Detail::compact);
    }
    return io::text::Detail::compact;
  }

public:
  using parameters = meta::list<olb::parameters::PHYS_ITER_T>;

  template <typename CASE, typename PARAMETERS>
  void report(CASE& myCase, std::size_t iT, PARAMETERS& params) {
    const std::size_t every = interval(myCase, params);
    if (every == 0 || iT % every != 0) {
      return;
    }
    emit(myCase, iT, detail());
  }

  template <typename CASE, typename PARAMETERS>
  void finalize(CASE& myCase, std::size_t iT, PARAMETERS& params) {
    if (interval(myCase, params) == 0) {
      return;
    }
    OstreamManager clout(std::cout, "Statistics");
    clout << std::endl << "-- Final report --" << std::endl;
    emit(myCase, iT, io::text::Detail::full);
  }
};

}

}

#endif
