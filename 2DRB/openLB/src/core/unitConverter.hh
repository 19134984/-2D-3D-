/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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


#ifndef UNITCONVERTER_HH
#define UNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "io/textTable.h"
#include "unitConverter.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const
{
  using Align = io::text::Table::Align;
  clout << io::text::section("Unit Converter") << std::endl;
  io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  parameters.row("Resolution",                  "N",            this->getResolution());
  parameters.row("Lattice velocity",            "latticeU",     this->getCharLatticeVelocity());
  parameters.row("Lattice relaxation frequency","omega",        this->getLatticeRelaxationFrequency());
  parameters.row("Lattice relaxation time",     "tau",          this->getLatticeRelaxationTime());
  parameters.row("Char. length [m]",            "charL",        this->getCharPhysLength());
  parameters.row("Char. speed [m/s]",           "charU",        this->getCharPhysVelocity());
  parameters.row("Kinematic viscosity [m^2/s]", "charNu",       this->getPhysViscosity());
  parameters.row("Density [kg/m^d]",            "charRho",      this->getPhysDensity());
  parameters.row("Char. pressure [N/m^2]",      "charPressure", this->getCharPhysPressure());
  parameters.row("Mach number",                 "Ma",           this->getMachNumber());
  parameters.row("Reynolds number",             "Re",           this->getReynoldsNumber());
  parameters.row("Knudsen number",              "Kn",           this->getKnudsenNumber());
  parameters.row("Char. CFL number",            "CFL",          this->getCharCFLnumber());
  parameters.render(clout);

  clout << std::endl;
  clout << io::text::section("Conversion Factors") << std::endl;
  io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
  factors.row("Voxel length [m]",     "physDeltaX",    this->getConversionFactorLength());
  factors.row("Time step [s]",        "physDeltaT",    this->getConversionFactorTime());
  factors.row("Velocity [m/s]",       "physVelocity",  this->getConversionFactorVelocity());
  factors.row("Density [kg/m^3]",     "physDensity",   this->getConversionFactorDensity());
  factors.row("Mass [kg]",            "physMass",      this->getConversionFactorMass());
  factors.row("Viscosity [m^2/s]",    "physViscosity", this->getConversionFactorViscosity());
  factors.row("Force [N]",            "physForce",     this->getConversionFactorForce());
  factors.row("Pressure [N/m^2]",     "physPressure",  this->getConversionFactorPressure());
  factors.render(clout);

  // Print thermal section if thermal parameters are initialized
  if (this->_latticeThermalRelaxationTime.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Thermal Parameters") << std::endl;
    io::text::Table thermalParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    thermalParameters.row("Thermal lattice relaxation time", "tau_AD",  this->getLatticeThermalRelaxationTime());
    thermalParameters.row("Thermal conductivity [W/m/K]",    "k",       this->getThermalConductivity());
    thermalParameters.row("Specific heat capacity [J/kg/K]", "cp",      this->getPhysSpecificHeatCapacity());
    thermalParameters.row("Thermal expansion coeff. [K^-1]", "beta",    this->getPhysThermalExpansionCoefficient());
    thermalParameters.row("Char. low temperature [K]",       "T_low",   this->getCharPhysLowTemperature());
    thermalParameters.row("Char. high temperature [K]",      "T_high",  this->getCharPhysHighTemperature());
    thermalParameters.row("Prandtl number",                  "Pr",      this->getPrandtlNumber());
    thermalParameters.row("Rayleigh number",                 "Ra",      this->getRayleighNumber());
    thermalParameters.render(clout);

    clout << std::endl;
    clout << io::text::section("Thermal Conversion Factors") << std::endl;
    io::text::Table thermalFactors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
    thermalFactors.row("Temperature [K]",             "deltaT",              this->getConversionFactorTemperature());
    thermalFactors.row("Thermal diffusivity [m^2/s]", "alpha",               this->getConversionFactorThermalDiffusivity());
    thermalFactors.row("Specific heat capacity [J/kg/K]","cp",               this->getConversionFactorSpecificHeatCapacity());
    thermalFactors.row("Thermal conductivity [W/m/K]", "k",                  this->getConversionFactorThermalConductivity());
    thermalFactors.row("Heat flux [W/m^2]",           "q",                   this->getConversionFactorHeatFlux());
    thermalFactors.render(clout);
  }

  // Print ADE section if ADE parameters are initialized
  if (this->_latticeAdeRelaxationTime.hasValue()) {
    clout << std::endl;
    clout << io::text::section("ADE Parameters") << std::endl;
    io::text::Table adeParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    adeParameters.row("ADE lattice relaxation time",     "tau_ADE", this->getLatticeAdeRelaxationTime());
    adeParameters.row("Diffusivity [m^2/s]",             "D",       this->getPhysDiffusivity());
    adeParameters.row("Peclet number",                   "Pe",      this->getPecletNumber());
    adeParameters.render(clout);
  }

  // Print multiphase section if multiphase parameters are initialized
  if (this->_physSurfaceTension.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Multiphase Parameters") << std::endl;
    io::text::Table multiphaseParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    multiphaseParameters.row("EoS parameter a",         "a",        this->getPhysEoSa());
    multiphaseParameters.row("EoS parameter b",         "b",        this->getPhysEoSb());
    multiphaseParameters.row("Molar mass [kg/mol]",     "M",        this->getPhysMolarMass());
    multiphaseParameters.row("Surface tension [N/m]",   "sigma",    this->getPhysSurfaceTension());
    multiphaseParameters.row("Characteristic temp [K]", "T",        this->getCharPhysTemperature());
    multiphaseParameters.row("Lattice surface tension", "sigma*",   this->getLatticeSurfaceTension());
    multiphaseParameters.render(clout);

    clout << std::endl;
    clout << io::text::section("Multiphase Conversion Factors") << std::endl;
    io::text::Table multiphaseFactors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
    multiphaseFactors.row("EoS a",             "convEoSa", this->getConversionFactorEoSa());
    multiphaseFactors.row("EoS b",             "convEoSb", this->getConversionFactorEoSb());
    multiphaseFactors.row("Molar mass",        "convM",    this->getConversionFactorMolarMass());
    multiphaseFactors.row("Surface tension",   "convSigma",this->getConversionFactorSurfaceTension());
    multiphaseFactors.row("Gas constant",      "R",        this->getConversionFactorGasConstant());
    multiphaseFactors.row("Chemical potential","convMu",   this->getConversionFactorChemicalPotential());
    multiphaseFactors.render(clout);
  }

  // Print power law section if power law parameters are initialized
  if (this->_powerLawIndex.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Power Law Parameters") << std::endl;
    io::text::Table powerLawParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    powerLawParameters.row("Consistency coefficient", "K",        this->getPhysConsistencyCoeff());
    powerLawParameters.row("Power law index",         "n",        this->getPowerLawIndex());
    powerLawParameters.row("Lattice consistency coeff","K*",      this->getLatticeConsistencyCoeff());
    powerLawParameters.render(clout);
  }

  // Print radiative section if radiative parameters are initialized
  if (this->_physAbsorption.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Radiative Parameters") << std::endl;
    io::text::Table radiativeParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    radiativeParameters.row("Phys. absorption [1/m]",      "mu_a",     this->getPhysAbsorption());
    radiativeParameters.row("Phys. scattering [1/m]",      "mu_s",     this->getPhysScattering());
    radiativeParameters.row("Extinction [1/m]",            "mu_t",     this->getExtinction());
    radiativeParameters.row("Effective attenuation [1/m]", "mu_eff",   this->getEffectiveAttenuation());
    radiativeParameters.row("Phys. diffusion [m]",         "D",        this->getPhysDiffusion());
    radiativeParameters.row("Single scattering albedo",    "albedo",   this->getScatteringAlbedo());
    radiativeParameters.row("Anisotropy factor",           "g",        this->getAnisotropyFactor());
    radiativeParameters.row("Refractive relative",         "n_rel",    this->getRefractiveRelative());
    radiativeParameters.render(clout);

    clout << std::endl;
    clout << io::text::section("Radiative Lattice Quantities") << std::endl;
    io::text::Table radiativeLattice({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    radiativeLattice.row("Lattice diffusion",  "D*",       this->getLatticeDiffusion());
    radiativeLattice.row("Lattice absorption", "abs*",     this->getLatticeAbsorption());
    radiativeLattice.row("Lattice scattering", "scat*",    this->getLatticeScattering());
    radiativeLattice.render(clout);
  }

  // Print structural mechanics section if structural mechanics parameters are initialized
  if (this->_youngsModulus.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Structural Mechanics Parameters") << std::endl;
    io::text::Table structuralmechanicsParameter({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    structuralmechanicsParameter.row("Char. length [m]",        "charL",    this->getCharPhysLength());
    structuralmechanicsParameter.row("Char. time [s]",          "charT",    this->getCharPhysTime());
    structuralmechanicsParameter.row("Char. displacement [m]",  "charU",    this->getCharPhysDisplacement());
    structuralmechanicsParameter.row("Young's modulus (Phys)",  "E",        this->getPhysYoungsModulus());
    structuralmechanicsParameter.row("Young's modulus (Latt)",  "E*",       this->getLatticeYoungsModulus());
    structuralmechanicsParameter.row("Bulk modulus (Phys)",     "K",        this->getBulkModulus());
    structuralmechanicsParameter.row("Bulk modulus (Latt)",     "K*",       this->getLatticeBulkModulus());
    structuralmechanicsParameter.row("Shear modulus (Phys)",    "mue",      this->getLameMu());
    structuralmechanicsParameter.row("Shear modulus (Latt)",    "mue*",     this->getLatticeLameMu());
    structuralmechanicsParameter.row("Lambda (Phys)",           "lambda",   this->getLameLambda());
    structuralmechanicsParameter.row("Lambda (Latt)",           "lambda*",  this->getLatticeLameLambda());
    structuralmechanicsParameter.row("Poisson ratio",           "nu",       this->getPoissonRatio());
    structuralmechanicsParameter.row("Solid Case",              "Case",     this->getSolidCase());
    structuralmechanicsParameter.render(clout);
  }

  // Print Stationary Linear Elasticity MRT section if Stationary Linear Elasticity MRT parameters are initialized
  if (this->_dampingFactor.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Stationary Linear Elasticity MRT Parameters") << std::endl;
    io::text::Table stationaryLinElaMRTParameter({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    stationaryLinElaMRTParameter.row("Damping factor",          "kappa",    this->getDampingFactor());
    stationaryLinElaMRTParameter.render(clout);
  }

  // Print Linear Elastodynamics Vector LBM section if Linear Elastodynamics Vector LBM parameters are initialized
  if (this->_latticeSpeed.hasValue()) {
    clout << std::endl;
    clout << io::text::section("Linear Elastodynamics VLBM Parameters") << std::endl;
    io::text::Table elastodynamicsVLBMParameter({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
    elastodynamicsVLBMParameter.row("cK",           "cK",   this->getCK());
    elastodynamicsVLBMParameter.row("cK (Lattice)", "cK*",  this->getLatticeCK());
    elastodynamicsVLBMParameter.row("cMu",          "cMu",  this->getCMu());
    elastodynamicsVLBMParameter.row("cMu (Lattice)","cMu*", this->getLatticeCMu());
    elastodynamicsVLBMParameter.row("Lattice speed","c",    this->getLatticeSpeed());
    elastodynamicsVLBMParameter.render(clout);
  }

  // Stability warning
  if (this->getLatticeRelaxationTime() < T(0.55) && this->getCharLatticeVelocity() > (T(8) * (this->getLatticeRelaxationTime() - T(0.5)) + T(1.e-8)) ) {
    T tauStable = this->getCharLatticeVelocity() / T(8) + T(0.5);
    T timeToCell = this->getConversionFactorTime() / this->getConversionFactorLength();
    if (this->getCharLatticeVelocity() >= T(0.3)) {
      tauStable = T(0.15) / T(8) + T(0.5);
      timeToCell = T(0.15) / this->getCharPhysVelocity();
    }
    T dxNew = timeToCell * this->getPhysViscosity() * descriptors::invCs2<T,DESCRIPTOR>() / (tauStable - T(0.5));
    T dtNew = dxNew * timeToCell;
    clout << std::endl;
    clout << "WARNING:" << std::endl;
    clout << "Potentially UNSTABLE combination of relaxation time (tau=" << this->getLatticeRelaxationTime() << ")" << std::endl;
    clout << "and characteristical CFL number (lattice velocity) charCFLnumber=" << this->getCharCFLnumber() << "!" << std::endl;
    clout << "Potentially maximum characteristical CFL number (maxCharCFLnumber=" << T(8)*(this->getLatticeRelaxationTime() - T(0.5)) << ")" << std::endl;
    clout << "Actual characteristical CFL number (charCFLnumber=" << this->getCharCFLnumber() << ") > " << T(8)*(this->getLatticeRelaxationTime() - T(0.5)) << std::endl;
    if (this->getCharLatticeVelocity() >= T(0.3)) {
      clout << "Please make the CFL number smaller than 0.3!" << std::endl;
    }
    clout << "Please reduce the the cell size or the time step size!" << std::endl;
    if (this->getCharLatticeVelocity() >= T(0.3)) {
      clout << "We recommend to use the cell size of " << dxNew << " m and the time step size of " << dtNew << " s (CFL = 0.15)." << std::endl;
    } else {
      clout << "We recommend to use the cell size of " << dxNew << " m and the time step size of " << dtNew << " s." << std::endl;
    }
    clout << "-------------------------------------------------------------" << std::endl;
  }

}

template <typename T, class DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::print() const
{
  OstreamManager clout(std::cout, "UnitConverter");
  print(clout);
}

template <typename T, typename DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const
{
  OstreamManager clout(std::cout, "UnitConverter");
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

template<typename T, typename DESCRIPTOR>
UnitConverter<T, DESCRIPTOR>* createUnitConverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createUnitConverter");

  params.setWarningsOn(false);

  T physDeltaX{};
  T physDeltaT{};

  T charPhysLength{};
  T charPhysVelocity{};
  T physViscosity{};
  T physDensity{};
  T charPhysPressure = 0;

  int resolution{};
  T latticeRelaxationTime{};
  T charLatticeVelocity{};

  int counter = 0; // counting the number of Discretization parameters and returning a warning if more than 2 are provided

  // params[parameter].read(value) sets the value or returns false if the parameter can not be found
  params["Application"]["PhysParameters"]["CharPhysLength"].read(charPhysLength);
  params["Application"]["PhysParameters"]["CharPhysVelocity"].read(charPhysVelocity);
  params["Application"]["PhysParameters"]["PhysViscosity"].read(physViscosity);
  params["Application"]["PhysParameters"]["PhysDensity"].read(physDensity);
  params["Application"]["PhysParameters"]["CharPhysPressure"].read(charPhysPressure);

  std::vector<std::string> discretizationParam = {"PhysDeltaX", "Resolution",
    "CharLatticeVelocity", "PhysDeltaT", "LatticeRelaxationTime"};

  for (std::size_t i = 0; i < discretizationParam.size(); i++) {
    std::string test;
    if (params["Application"]["Discretization"][discretizationParam[i]].read(test, false)) {
      counter++;
    }
  }

  if (counter > 2) {
    clout << "WARNING: More than 2 discretization parameters provided" << std::endl;
  }

  if (!params["Application"]["Discretization"]["PhysDeltaX"].read(physDeltaX,false)) {
    if (!params["Application"]["Discretization"]["Resolution"].read<int>(resolution,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaX, resolution or charLatticeVelocity
        throw std::runtime_error("Error: Have not found PhysDeltaX, Resolution or CharLatticeVelocity in XML file.");

      }
      else {
        // found charLatticeVelocity
        if (params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
          physDeltaX = charPhysVelocity / charLatticeVelocity * physDeltaT;
        }
        else if (params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          physDeltaX = physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5);
        }
        else {
          throw std::runtime_error("Error: Only found CharLatticeVelocity, missing PhysDeltaT or LatticeRelaxationTime");
        }
      }
    }
    else {
      // found resolution
      physDeltaX = charPhysLength / resolution;

      if (params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        latticeRelaxationTime = physViscosity * charLatticeVelocity * descriptors::invCs2<T,DESCRIPTOR>() * resolution + 0.5;
      }
      else {
        if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          throw std::runtime_error("Error: Have not found LatticeRelaxationTime and was not able to derive it using CharLatticeVelocity");
        }
      }
    }
  }
  // found physDeltaX
  if (!params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
    if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaT, latticeRelaxationTime and charLatticeVelocity
        throw std::runtime_error("Error: Have not found PhysDeltaT, LatticeRelaxationTime or CharLatticeVelocity in XML file.");
      }
      else {
        // found charLatticeVelocity
        physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
      }
    }
    else {
      // found latticeRelaxationTime
      physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
    }
  }

  return new UnitConverter<T, DESCRIPTOR>(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physViscosity, physDensity, charPhysPressure);
}

}  // namespace olb

#endif
