/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Florian Kaiser
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

#ifndef STRUCTURAL_MECHANICS_UNIT_CONVERTER_HH
#define STRUCTURAL_MECHANICS_UNIT_CONVERTER_HH

#include "structuralMechanicsUnitConverter.h"
#include "io/textTable.h"

namespace olb {

  template <typename T, typename DESCRIPTOR>
  void StructuralMechanicsUnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const {
    using Align = io::text::Table::Align;

    clout << io::text::section("Unit Converter") << std::endl;
    io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    parameters.row("Resolution",                    "N",      this->getResolution());
    parameters.row("Voxel length [m]",              "dx",     this->getConversionFactorLength());
    parameters.row("Time step [s]",                 "dt",     this->getConversionFactorTime());
    parameters.row("Poisson ratio",                 "nu",     this->getPoissonRatio());
    parameters.row("Young's modulus [N/m^2]",       "E",      this->getYoungsModulus());
    parameters.row("Lame shear modulus [N/m^2]",    "mu",     this->getLameMu());
    parameters.row("Lame lambda [N/m^2]",           "lambda", this->getLameLambda());
    parameters.row("Bulk modulus [N/m^2]",          "K",      this->getBulkModulus());
    parameters.row("Density [kg/m^d]",              "rho",    this->getPhysDensity());
    parameters.render(clout);

    clout << std::endl;
    clout << io::text::section("Characteristic Units") << std::endl;
    io::text::Table characteristic({{"Quantity"}, {"Symbol"}, {"Value", Align::Right}});
    characteristic.row("Char. time [s]",   "charT", this->getCharPhysTime());
    characteristic.row("Char. length [m]", "charL", this->getCharPhysLength());
    characteristic.render(clout);
  }

  template <typename T, class DESCRIPTOR>
  void StructuralMechanicsUnitConverter<T, DESCRIPTOR>::print() const {
    OstreamManager clout(std::cout, "StructuralMechanicsUnitConverter");
    print(clout);
  }

  template <typename T, typename DESCRIPTOR>
  void StructuralMechanicsUnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const {
    OstreamManager clout(std::cout, "StructuralMechanicsUnitConverter");
    std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";
    if (singleton::mpi().isMainProcessor()) {
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        clout << "error write() function: can not open std::ofstream" << std::endl;
      }
      else {
        print( fout );
        fout.close();
      }
    }
  }

  template <typename T, typename DESCRIPTOR>
  void ElastodynamicsVectorLBMUnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const {
    using Align = io::text::Table::Align;

    clout << io::text::section("Unit Converter") << std::endl;
    io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    parameters.row("Resolution",                    "N",      this->getResolution());
    parameters.row("Voxel length [m]",              "dx",     this->getConversionFactorLength());
    parameters.row("Time step [s]",                 "dt",     this->getConversionFactorTime());
    parameters.row("Poisson ratio",                 "nu",     this->getPoissonRatio());
    parameters.row("Young's modulus [N/m^2]",       "E",      this->getYoungsModulus());
    parameters.row("Lame shear modulus [N/m^2]",    "mu",     this->getLameMu());
    parameters.row("Lame lambda [N/m^2]",           "lambda", this->getLameLambda());
    parameters.row("Bulk modulus [N/m^2]",          "K",      this->getBulkModulus());
    parameters.row("Density [kg/m^d]",              "rho",    this->getPhysDensity());
    parameters.row("Dilatational wave speed [m/s]", "cK",     this->getCK());
    parameters.row("Shear wave speed [m/s]",        "cMu",    this->getCMu());
    parameters.render(clout);

    clout << std::endl;
    clout << io::text::section("Lattice Units") << std::endl;
    io::text::Table lattice({{"Quantity"}, {"Symbol"}, {"Value", Align::Right}});
    lattice.row("Young's modulus",         "E",      this->getLatticeYoungsModulus());
    lattice.row("Lame shear modulus",      "mu",     this->getLatticeLameMu());
    lattice.row("Lame lambda",             "lambda", this->getLatticeLameLambda());
    lattice.row("Bulk modulus",            "K",      this->getLatticeBulkModulus());
    lattice.row("Dilatational wave speed", "cK",     this->getLatticeCK());
    lattice.row("Shear wave speed",        "cMu",    this->getLatticeCMu());
    lattice.render(clout);

    clout << std::endl;
    clout << io::text::section("Characteristic Units") << std::endl;
    io::text::Table characteristic({{"Quantity"}, {"Symbol"}, {"Value", Align::Right}});
    characteristic.row("Char. time [s]",    "charT", this->getCharPhysTime());
    characteristic.row("Char. length [m]",  "charL", this->getCharPhysLength());
    characteristic.row("Char. speed [m/s]", "charV", this->getCharPhysVelocity());
    characteristic.render(clout);

    clout << std::endl;
    clout << io::text::section("Conversion Factors") << std::endl;
    io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
    factors.row("Voxel length [m]",    "physDeltaX",       this->getConversionFactorLength());
    factors.row("Time step [s]",       "physDeltaT",       this->getConversionFactorTime());
    factors.row("Lattice speed [m/s]", "latticeSpeed",     this->getLatticeSpeed());
    factors.row("Velocity [m/s]",      "physVelocity",     this->getConversionFactorVelocity());
    factors.row("Displacement [m]",    "physDisplacement", this->getConversionFactorDisplacement());
    factors.row("Stress [N/m^2]",      "physStress",       this->getConversionFactorSolidStress());
    factors.row("Traction [N/m^2]",    "physTraction",     this->getConversionFactorTraction());
    factors.row("Force [N]",           "physForce",        this->getConversionFactorSolidForce());
    factors.render(clout);

    clout << std::endl;
    clout << io::text::section("Stability") << std::endl;
    io::text::Table stability({{"Criterion"}, {"Value", Align::Right}});
    stability.row("Simulation is stable?",      this->template elastodynamicsIsStable<DESCRIPTOR>());
    stability.row("2 * sqrt(cK^2 + cMu^2) / c", this->template elastodynamicsGetStability<DESCRIPTOR>());
    stability.row("Solid case",                 this->getSolidCase());
    stability.render(clout);
  }

  template <typename T, class DESCRIPTOR>
  void ElastodynamicsVectorLBMUnitConverter<T, DESCRIPTOR>::print() const {
    OstreamManager clout(std::cout, "ElastodynamicsVectorLBMUnitConverter");
    print(clout);
  }

  template <typename T, typename DESCRIPTOR>
  void ElastodynamicsVectorLBMUnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const {
    OstreamManager clout(std::cout, "ElastodynamicsVectorLBMUnitConverter");
    std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";
    if (singleton::mpi().isMainProcessor()) {
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        clout << "error write() function: can not open std::ofstream" << std::endl;
      }
      else {
        print( fout );
        fout.close();
      }
    }
  }

  template <typename T, typename DESCRIPTOR>
  void StationaryMRTUnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const {
    using Align = io::text::Table::Align;

    clout << io::text::section("Unit Converter") << std::endl;
    io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    parameters.row("Resolution",                 "N",      this->getResolution());
    parameters.row("Char. length [m]",           "charL",  this->getCharPhysLength());
    parameters.row("Char. time [s]",             "charT",  this->getCharPhysTime());
    parameters.row("Char. displacement [m]",     "charU",  this->getCharPhysDisplacement());
    parameters.row("Lame shear modulus [N/m^2]", "mu",     this->getLameMu());
    parameters.row("Lame lambda [N/m^2]",        "lambda", this->getLameLambda());
    parameters.row("Bulk modulus [N/m^2]",       "K",      this->getBulkModulus());
    parameters.row("Young's modulus [N/m^2]",    "E",      this->getYoungsModulus());
    parameters.row("Poisson ratio",              "nu",     this->getPoissonRatio());
    parameters.row("Damping factor",             "kappa",  this->getDampingFactor());
    parameters.render(clout);

    clout << std::endl;
    clout << io::text::section("Lattice Units") << std::endl;
    io::text::Table lattice({{"Quantity"}, {"Symbol"}, {"Value", Align::Right}});
    lattice.row("Lame shear modulus", "mu",     this->getLatticeLameMu());
    lattice.row("Lame lambda",        "lambda", this->getLatticeLameLambda());
    lattice.row("Bulk modulus",       "K",      this->getLatticeBulkModulus());
    lattice.row("Young's modulus",    "E",      this->getLatticeYoungsModulus());
    lattice.render(clout);

    clout << std::endl;
    clout << io::text::section("Conversion Factors") << std::endl;
    io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
    factors.row("Voxel length [m]", "physDeltaX",   this->getConversionFactorLength());
    factors.row("Time step [s]",    "physDeltaT",   this->getConversionFactorTime());
    factors.row("Stress [N/m^2]",   "physStress",   this->getConversionFactorSolidStress());
    factors.row("Force [N]",        "physForce",    this->getConversionFactorSolidForce());
    factors.row("Traction [N/m^2]", "physTraction", this->getConversionFactorTraction());
    factors.render(clout);
  }

  template <typename T, class DESCRIPTOR>
  void StationaryMRTUnitConverter<T, DESCRIPTOR>::print() const {
    OstreamManager clout(std::cout, "StationaryMRTUnitConverter");
    print(clout);
  }

  template <typename T, typename DESCRIPTOR>
  void StationaryMRTUnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const {
    OstreamManager clout(std::cout, "StationaryMRTUnitConverter");
    std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";
    if (singleton::mpi().isMainProcessor()) {
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        clout << "error write() function: can not open std::ofstream" << std::endl;
      }
      else {
        print( fout );
        fout.close();
      }
    }
  }

} // namespace olb

#endif
