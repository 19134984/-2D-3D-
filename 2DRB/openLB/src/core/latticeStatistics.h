/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * Interface for post-processing steps -- header file.
 */
#ifndef LATTICE_STATISTICS_H
#define LATTICE_STATISTICS_H

#include <vector>
#include <string>
#include <memory>
#include <limits>
#include <concepts>
#include <map>
#include <typeindex>
#include <functional>
#include "io/ostreamManager.h"
#include "io/gnuplotWriter.h"

namespace olb {

template<typename T, typename DESCRIPTOR> class SuperLattice;
template<typename T, typename DESCRIPTOR> class LatticeStatistics;

/////////////////// Statistics Postprocessing ////////////////////////////

namespace concepts {

/// A parameter whose value is sourced from lattice statistics via compute(stats)
template <typename FIELD, typename T, typename DESCRIPTOR>
concept StatisticsField = requires (const LatticeStatistics<T,DESCRIPTOR>& stats) {
  { FIELD::compute(stats) } -> std::convertible_to<T>;
};

}

namespace statistics {

/// Lattice-wide average density, consumed by ConstRhoBGK
struct AVERAGE_RHO : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static T compute(const LatticeStatistics<T,DESCRIPTOR>& stats) {
    return stats.getAverageRho();
  }
};

/// Lattice-wide average kinetic energy
struct AVERAGE_ENERGY : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static T compute(const LatticeStatistics<T,DESCRIPTOR>& stats) {
    return stats.getAverageEnergy();
  }
};

/// Maximum lattice velocity, reported as maxCFL
struct MAX_U : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static T compute(const LatticeStatistics<T,DESCRIPTOR>& stats) {
    return stats.getMaxU();
  }
};

/// Parameters sourced live from lattice statistics at every collide
using live_parameters = meta::list<AVERAGE_RHO>;

}

/// Min/max/mean of one scalar or vector quantity over a lattice
/**
 * One value per component. Populated by the built-in collision-time reduction
 * (rho, u) and by the registered-field reductions.
 **/
template<typename T>
struct StatisticEntry {
  std::string    name;
  unsigned       components = 1;
  std::vector<T> min;
  std::vector<T> max;
  std::vector<T> mean;
};

/// Collects statistics (min/max/mean of rho, u and registered fields) of a lattice
/**
 * rho and the velocity magnitude u are built-in and reduced at collision time.
 * Arbitrary descriptor fields are added via registerField<FIELD>() and reduced
 * at statistics time (see SuperLattice::collectStatistics). The scalar getters
 * getAverageRho / getAverageEnergy / getMaxU / getNumCells expose the built-in
 * aggregates, and print() emits the avgEnergy / maxCFL Gnuplot output.
 *
 * \ingroup core
 **/
template<typename T, typename DESCRIPTOR>
class LatticeStatistics {
public:
  /// Per-block partial reduction of the built-in quantities (rho, u)
  /**
   * Kept an aggregate (no user constructors) so the OpenMP reduction initializer
   * `omp_priv={}` seeds each thread with the min/max identities below. avRho /
   * avEnergy are running sums (averaged in reset); maxU keeps max(uSqr) (turned
   * into max|u| in reset). minRho / maxRho / minUSqr / sumU add the remaining
   * pieces for min/max/mean of rho and |u|. The Expr guard mirrors the original
   * maxU handling for the CSE expression value type.
   **/
  struct Aggregatable {
    std::size_t nCells = 0;
    T avRho    { };
    T avEnergy { };
    T maxU     { };
    T minRho  = std::numeric_limits<T>::max();
    T maxRho  = std::numeric_limits<T>::lowest();
    T minUSqr = std::numeric_limits<T>::max();
    T sumU     { };

    void increment(T rho, T uSqr) {
      nCells   += 1;
      avRho    += rho;
      avEnergy += uSqr;
      if constexpr (!std::is_same_v<Expr,T>) {
        maxU    = std::max(uSqr, maxU);
        minRho  = std::min(rho, minRho);
        maxRho  = std::max(rho, maxRho);
        minUSqr = std::min(uSqr, minUSqr);
        sumU   += util::sqrt(uSqr);
      }
    }

    Aggregatable& operator+=(const Aggregatable& rhs) {
      nCells   += rhs.nCells;
      avRho    += rhs.avRho;
      avEnergy += rhs.avEnergy;
      if constexpr (!std::is_same_v<Expr,T>) {
        maxU    = std::max(maxU, rhs.maxU);
        minRho  = std::min(minRho, rhs.minRho);
        maxRho  = std::max(maxRho, rhs.maxRho);
        minUSqr = std::min(minUSqr, rhs.minUSqr);
        sumU   += rhs.sumU;
      }
      return *this;
    }
  };

  /// Raw per-block/instance reduction of one registered field (count + per-component min/max/sum)
  struct FieldResult {
    unsigned       components = 0;
    std::size_t    count      = 0;
    std::vector<T> min;
    std::vector<T> max;
    std::vector<T> sum;
  };

  /// A field registered for reduction: its identity, name and the operator installer
  struct RegisteredField {
    std::type_index id;
    std::string     name;
    unsigned        components;
    std::function<void(SuperLattice<T,DESCRIPTOR>&)> addOperator;
  };

public:
  LatticeStatistics();
  ~LatticeStatistics() = default;

  /// Finalized built-in results, exchanged during super<->block aggregation
  struct Result {
    T avRho    { };
    T avEnergy { };
    T maxU     { };
    T minRho   { };
    T maxRho   { };
    T minU     { };
    T meanU    { };
    std::size_t nCells = 0;
  };

  /// Finalize the partial accumulators into results for this interval
  void reset();
  /// Read the finalized built-in results
  Result getResult() const;
  /// Set finalized built-in results directly (super aggregation and broadcast)
  void setResult(const Result& result);
  void initialize();

  /// Accumulate one cell into the built-in (rho, u) partials
  void incrementStats(T rho, T uSqr);
  /// Accumulate a per-block partial into the running built-in totals
  void incrementStats(Aggregatable& aggregatable);
  /// Increment cell count only
  void incrementStats();

  /// Register descriptor FIELD for on-device min/max/mean reduction at statistics time
  /** Defined in statisticsFieldReducer.h; installs a FieldStatisticO<FIELD> post-processor. **/
  template<typename FIELD>
  void registerField(std::string name = "");
  /// Registered fields (identity/name/installer), driven by SuperLattice
  const std::vector<RegisteredField>& getRegisteredFields() const { return _registered; }
  std::size_t getFieldCount() const { return _registered.size(); }
  /// Store a field's raw reduction result (called by FieldStatisticO on each block)
  template<typename FIELD>
  void setFieldResult(const FieldResult& result) { _fieldResults[typeid(FIELD)] = result; }
  /// Raw per-field results by field type (written by the reduction operator)
  std::map<std::type_index,FieldResult>& fieldResults() { return _fieldResults; }
  const std::map<std::type_index,FieldResult>& fieldResults() const { return _fieldResults; }

  /// Finalized built-in scalar aggregates
  T getAverageRho() const    { return _avRho; }
  T getAverageEnergy() const { return _avEnergy; }
  T getMaxU() const          { return _maxU; }
  std::size_t getNumCells() const { return _numCells; }

  /// Unified min/max/mean access by quantity name ("rho", "u", or a field name)
  T getMin(const std::string& name, unsigned iComponent=0) const;
  T getMax(const std::string& name, unsigned iComponent=0) const;
  T getMean(const std::string& name, unsigned iComponent=0) const;
  /// All quantities (rho, u, registered) as unified entries for print/CSV
  std::vector<StatisticEntry<T>> getEntries() const;

  void incrementTime();
  void resetTime(size_t value=0);
  std::size_t getTime() const;

  /// Stable identity used to name CSV / plot output when no explicit name is given
  void setName(const std::string& name) { _name = name; }
  const std::string& getName() const { return _name; }

  /// Single-line console output with compact registered-field tokens and the avgEnergy/maxCFL Gnuplot plots
  void print(int iterationStep, T physicalTime=-1, const std::string& postfix = "");
  void printWithoutPlot(int iterationStep, T physicalTime=-1) const;
  /// Formatted min/max/mean table over all quantities
  void printSummary(int iterationStep, T physicalTime=-1) const;
  /// Append the current statistics row to ./tmp/statisticsData/<name>_statistics.csv
  void writeCSV(int iterationStep, T physicalTime=-1, const std::string& name="") const;

private:
  const RegisteredField* findRegistered(const std::string& name) const;
  StatisticEntry<T> entryFor(const RegisteredField& rf) const;

  mutable OstreamManager clout;

  Aggregatable _partial;

  T _avRho    { };
  T _avEnergy { };
  T _maxU     { };
  T _minRho   { };
  T _maxRho   { };
  T _minU     { };
  T _meanU    { };
  std::size_t _numCells    = 0;
  std::size_t _latticeTime = 0;
  bool _firstCall = true;

  std::vector<RegisteredField> _registered;
  std::map<std::type_index,FieldResult> _fieldResults;

  static std::size_t _instanceCounter;
  std::size_t _id;
  std::string _name;

  std::unique_ptr<Gnuplot<T>> gplotAvgEnergy;
  std::unique_ptr<Gnuplot<T>> gplotMaxCFL;
};

}  // namespace olb

#endif
