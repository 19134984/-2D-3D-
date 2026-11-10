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


#ifndef MULTIPHASEUNITCONVERTER_HH
#define MULTIPHASEUNITCONVERTER_HH

#include "multiPhaseUnitConverter.h"
#include "io/textTable.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void MultiPhaseUnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const
{
  using Align = io::text::Table::Align;
  clout << io::text::section("Multi-Phase Unit Converter") << std::endl;
  io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  parameters.row("Resolution",                  "N",                  this->getResolution());
  parameters.row("Lattice velocity",            "latticeU",           this->getCharLatticeVelocity());
  parameters.row("Lattice relaxation frequency","omega",              this->getLatticeRelaxationFrequency());
  parameters.row("Lattice relaxation time",     "tau",                this->getLatticeRelaxationTime());
  parameters.row("Char. length [m]",            "charL",              this->getCharPhysLength());
  parameters.row("Char. speed [m/s]",           "charU",              this->getCharPhysVelocity());
  parameters.row("Kinematic viscosity [m^2/s]", "charNu",             this->getPhysViscosity());
  parameters.row("Char. pressure [N/m^2]",      "charPressure",       this->getCharPhysPressure());
  parameters.row("Reynolds number",             "reynoldsNumber",     this->getReynoldsNumber());
  parameters.row("Phys. Delta X [m]",           "physDeltaX",         this->getPhysDeltaX());
  parameters.row("Phys. Delta T [s]",           "physDeltaT",         this->getPhysDeltaT());
  parameters.row("Surface tension [J/m^2]",     "physSurfaceTension", this->getPhysSurfaceTension());
  parameters.row("Char. temperature [K]",       "charPhysTemperature",this->getCharPhysTemperature());
  parameters.render(clout);

  clout << std::endl;
  clout << io::text::section("Conversion Factors") << std::endl;
  io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
  factors.row("Voxel length [m]",              "physDeltaX",         this->getConversionFactorLength());
  factors.row("Time step [s]",                 "physDeltaT",         this->getConversionFactorTime());
  factors.row("Velocity [m/s]",                "physVelocity",       this->getConversionFactorVelocity());
  factors.row("Density [kg/m^3]",              "physDensity",        this->getConversionFactorDensity());
  factors.row("Mass [kg]",                     "physMass",           this->getConversionFactorMass());
  factors.row("Viscosity [m^2/s]",             "physViscosity",      this->getConversionFactorViscosity());
  factors.row("Force [N]",                     "physForce",          this->getConversionFactorForce());
  factors.row("Pressure [N/m^2]",              "physPressure",       this->getConversionFactorPressure());
  factors.row("Equation of state a [Jm^3/mol^2]","physEoSa",         this->getConversionFactorEoSa());
  factors.row("Equation of state b [m^3/mol]", "physEoSb",           this->getConversionFactorEoSb());
  factors.row("Molar mass [kg/mol]",           "physMolarMass",      this->getConversionFactorMolarMass());
  factors.row("Surface tension [J/m^2]",       "physSurfaceTension", this->getConversionFactorSurfaceTension());
  factors.row("Temperature [K]",               "physTemperature",    this->getConversionFactorTemperature());
  factors.render(clout);
}

template <typename T, typename DESCRIPTOR>
void MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR>::print(std::ostream& clout) const
{
  using Align = io::text::Table::Align;
  clout << io::text::section("Parameters in Lattice Units") << std::endl;
  io::text::Table latticeParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  latticeParameters.row("Resolution",                  "N",        this->getResolution());
  latticeParameters.row("Lattice velocity",            "latticeU", this->getCharLatticeVelocity());
  latticeParameters.row("Lattice relaxation frequency","omega",    this->getLatticeRelaxationFrequency());
  latticeParameters.row("Lattice relaxation time",     "tau",      this->getLatticeRelaxationTime());
  latticeParameters.render(clout);

  clout << std::endl;
  clout << io::text::section("Parameters in Physical Units") << std::endl;
  io::text::Table physParameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  physParameters.row("Char. length [m]",            "charL",          this->getCharPhysLength());
  physParameters.row("Char. speed [m/s]",           "charU",          this->getCharPhysVelocity());
  physParameters.row("Kinematic viscosity [m^2/s]", "charNu",         this->getPhysViscosity());
  physParameters.row("Char. pressure [N/m^2]",      "charPressure",   this->getCharPhysPressure());
  physParameters.row("Reynolds number",             "reynoldsNumber", this->getReynoldsNumber());
  physParameters.row("Phys. Delta X [m]",           "physDeltaX",     this->getPhysDeltaX());
  physParameters.row("Phys. Delta T [s]",           "physDeltaT",     this->getPhysDeltaT());
  physParameters.render(clout);

  clout << std::endl;
  clout << io::text::section("Conversion Factors") << std::endl;
  io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
  factors.row("Voxel length [m]",        "physDeltaX",         this->getConversionFactorLength());
  factors.row("Time step [s]",           "physDeltaT",         this->getConversionFactorTime());
  factors.row("Velocity [m/s]",          "physVelocity",       this->getConversionFactorVelocity());
  factors.row("Density [kg/m^3]",        "physDensity",        this->getConversionFactorDensity());
  factors.row("Mass [kg]",               "physMass",           this->getConversionFactorMass());
  factors.row("Viscosity [m^2/s]",       "physViscosity",      this->getConversionFactorViscosity());
  factors.row("Force [N]",               "physForce",          this->getConversionFactorForce());
  factors.row("Pressure [N/m^2]",        "physPressure",       this->getConversionFactorPressure());
  factors.row("Surface tension [J/m^2]", "physSurfaceTension", this->getConversionFactorSurfaceTension());
  factors.render(clout);
}

}

#endif
