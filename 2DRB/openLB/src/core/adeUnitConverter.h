/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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
 * Unit conversion handling for Advection-Diffusion Problmes -- header file.
 */
#ifndef ADE_UNITCONVERTER_H
#define ADE_UNITCONVERTER_H

#include <math.h>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "io/ostreamManager.h"
#include "io/fileName.h"
#include "io/textTable.h"

#include "core/util.h"
#include "core/unitConverter.h"
#include "core/singleton.h"


namespace olb {

template <typename T, typename DESCRIPTOR>
struct AdeUnitConverter : public UnitConverter<T,DESCRIPTOR> {
  AdeUnitConverter(
    T physDeltaX,
    T physDeltaT,
    T charPhysLength,
    T charPhysVelocity,
    T physDiffusivity,
    T physDensity)
    : UnitConverter<T,DESCRIPTOR>(physDeltaX, physDeltaT,
                                  charPhysLength, charPhysVelocity,
                                  physDiffusivity, physDensity)
  {
    this->_physDiffusivity = physDiffusivity;
    this->_conversionDiffusivity = physDeltaX * physDeltaX / physDeltaT;
    this->_latticeAdeRelaxationTime = (physDiffusivity / this->_conversionDiffusivity * descriptors::invCs2<T,DESCRIPTOR>()) + 0.5;
  }

  void print(std::ostream& clout) const override;
};

template <typename T, class DESCRIPTOR>
void AdeUnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const
{
  using Align = io::text::Table::Align;
  clout << io::text::section("ADE Unit Converter") << std::endl;
  io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  parameters.row("Resolution",                  "N",             this->getResolution());
  parameters.row("Lattice velocity",            "latticeU",      this->getCharLatticeVelocity());
  parameters.row("Lattice relaxation frequency","omega",         this->getLatticeRelaxationFrequency());
  parameters.row("Lattice relaxation time",     "tau",           this->getLatticeRelaxationTime());
  parameters.row("Char. length [m]",            "charL",         this->getCharPhysLength());
  parameters.row("Char. speed [m/s]",           "charU",         this->getCharPhysVelocity());
  parameters.row("Density [kg/m^d]",            "charRho",       this->getPhysDensity());
  parameters.row("Peclet number",               "Pe",            this->getPecletNumber());
  parameters.row("Diffusivity [m^2/s]",         "charNu",        this->getPhysDiffusivity());
  parameters.row("Mach number",                 "machNumber",    this->getMachNumber());
  parameters.row("Knudsen number",              "knudsenNumber", this->getKnudsenNumber());
  parameters.render(clout);

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
}

template <typename T, typename DESCRIPTOR>
struct AdeUnitConverterFromResolutionAndRelaxationTime : public AdeUnitConverter<T, DESCRIPTOR> {
  AdeUnitConverterFromResolutionAndRelaxationTime(
    std::size_t resolution,
    T latticeRelaxationTime,
    T charPhysLength,
    T charPhysVelocity,
    T physDiffusivity,
    T physDensity): AdeUnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow((charPhysLength/resolution),2) / physDiffusivity,
        charPhysLength,
        charPhysVelocity,
        physDiffusivity,
        physDensity)
  { }
};

template <typename T, typename DESCRIPTOR>
struct AdeUnitConverterFromResolutionAndLatticeVelocity : public AdeUnitConverter<T, DESCRIPTOR> {
  AdeUnitConverterFromResolutionAndLatticeVelocity(
    std::size_t resolution,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity) : AdeUnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (charLatticeVelocity / charPhysVelocity * charPhysLength / resolution),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity)
  { }
};

}
#endif
