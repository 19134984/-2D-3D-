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

#ifndef THERMALUNITCONVERTER_HH
#define THERMALUNITCONVERTER_HH

#include "thermalUnitConverter.h"
#include "io/textTable.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR, typename ThermalLattice>
void ThermalUnitConverter<T, DESCRIPTOR, ThermalLattice>::print() const
{
  using Align = io::text::Table::Align;
  OstreamManager clout(std::cout, "ThermalUnitConverter");

  clout << io::text::section("Unit Converter") << std::endl;
  io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  parameters.row("Resolution",                      "N",              this->getResolution());
  parameters.row("Lattice velocity",                "latticeU",       this->getCharLatticeVelocity());
  parameters.row("Lattice relaxation frequency",    "omega",          this->getLatticeRelaxationFrequency());
  parameters.row("Lattice relaxation time",         "tau",            this->getLatticeRelaxationTime());
  parameters.row("Thermal lattice relaxation frequency", "omega_AD",  this->getLatticeThermalRelaxationFrequency());
  parameters.row("Thermal lattice relaxation time", "tau_AD",         this->getLatticeThermalRelaxationTime());
  parameters.row("Char. length [m]",                "charL",          this->getCharPhysLength());
  parameters.row("Char. speed [m/s]",               "charU",          this->getCharPhysVelocity());
  parameters.row("Kinematic viscosity [m^2/s]",     "charNu",         this->getPhysViscosity());
  parameters.row("Density [kg/m^d]",                "charRho",        this->getPhysDensity());
  parameters.row("Char. pressure [N/m^2]",          "charPressure",   this->getCharPhysPressure());
  parameters.row("Reynolds number",                 "reynoldsNumber", this->getReynoldsNumber());
  parameters.render(clout);

  clout << std::endl;
  clout << io::text::section("Thermal Parameters") << std::endl;
  io::text::Table thermalParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  thermalParameters.row("Phys. Delta X [m]",             "physDeltaX",                     this->getPhysDeltaX());
  thermalParameters.row("Phys. Delta T [s]",             "physDeltaT",                     this->getPhysDeltaT());
  thermalParameters.row("Char. pressure [N/m^2]",        "charPressure",                   this->getCharPhysPressure());
  thermalParameters.row("Thermal conductivity [W/m/K]",  "physThermalConductivity",        this->getThermalConductivity());
  thermalParameters.row("Specific heat capacity [J/kg/K]","physSpecificHeatCapacity",      this->getPhysSpecificHeatCapacity());
  thermalParameters.row("Thermal expansion coeff. [K^-1]","physThermalExpansionCoefficent",this->getPhysThermalExpansionCoefficient());
  thermalParameters.row("Char. low temperature [K]",     "charPhysLowTemperature",         this->getCharPhysLowTemperature());
  thermalParameters.row("Char. high temperature [K]",    "charPhysHighTemperature",        this->getCharPhysHighTemperature());
  thermalParameters.row("Prandtl number",                "prandtlNumber",                  this->getPrandtlNumber());
  thermalParameters.row("Rayleigh number",               "rayleighNumber",                 this->getRayleighNumber());
  thermalParameters.render(clout);

  clout << std::endl;
  clout << io::text::section("Conversion Factors") << std::endl;
  io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
  factors.row("Voxel length [m]",  "physDeltaX",    this->getConversionFactorLength());
  factors.row("Time step [s]",     "physDeltaT",    this->getConversionFactorTime());
  factors.row("Velocity [m/s]",    "physVelocity",  this->getConversionFactorVelocity());
  factors.row("Density [kg/m^3]",  "physDensity",   this->getConversionFactorDensity());
  factors.row("Mass [kg]",         "physMass",      this->getConversionFactorMass());
  factors.row("Viscosity [m^2/s]", "physViscosity", this->getConversionFactorViscosity());
  factors.row("Force [N]",         "physForce",     this->getConversionFactorForce());
  factors.row("Pressure [N/m^2]",  "physPressure",  this->getConversionFactorPressure());
  factors.render(clout);

  clout << std::endl;
  clout << io::text::section("Thermal Conversion Factors") << std::endl;
  io::text::Table thermalFactors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
  thermalFactors.row("Temperature [K]",             "temperature",             this->getConversionFactorTemperature());
  thermalFactors.row("Thermal diffusivity [m^2/s]", "physThermalDiffusity",    this->getConversionFactorThermalDiffusivity());
  thermalFactors.row("Specific heat capacity [J/kg]","physSpecificHeatCapacity",this->getConversionFactorSpecificHeatCapacity());
  thermalFactors.row("Thermal conductivity [W/m/K]","physThermalConductivity", this->getConversionFactorThermalConductivity());
  thermalFactors.row("Heat flux [W]",               "physHeatFlux",            this->getConversionFactorHeatFlux());
  thermalFactors.render(clout);
}

}

#endif
