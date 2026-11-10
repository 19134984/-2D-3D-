#ifndef OLB_APPS_FLORIAN_ADSORPTIONCONVERTER_H_
/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Florian Raichle
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
 * Unit conversion handling -- header file.
 */

#define OLB_APPS_FLORIAN_ADSORPTIONCONVERTER_H_

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

template<typename T, typename DESCRIPTOR>
struct AdsorptionConverter : public UnitConverter<T, DESCRIPTOR> {
  /** Converter for adsorption reactions primarily for ADE lattices based on concentration.
   * Renames viscosity to diffusivity and also takes fluid parameters.
   *
   * @param physDeltaX
   * @param physDeltaT
   * @param charPhysLength
   * @param charPhysVelocity
   * @param physDiffusivity
   * @param charConcentration
   * @param particleConcentration
   * @param fluidPhysViscosity
   */
  AdsorptionConverter(
      T physDeltaX,
      T physDeltaT,
      T charPhysLength,
      T charPhysVelocity,
      T physDiffusivity,
      T charConcentration,
      T particleConcentration,
      T fluidPhysViscosity) : UnitConverter<T, DESCRIPTOR>(physDeltaX,
                                                           physDeltaT,
                                                           charPhysLength,
                                                           charPhysVelocity,
                                                           physDiffusivity,
                                                           charConcentration)
  {
    this->_physViscosityAdsorption = fluidPhysViscosity;
    this->_conversionViscosityAdsorption = physDeltaX * physDeltaX / physDeltaT;
    this->_conversionViscosity = physDeltaX * physDeltaX / physDeltaT;
    this->_particleConcentrationAdsorption = particleConcentration;

    this->_physDiffusivity = this->getPhysViscosity();
    this->_conversionDiffusivity = this->getConversionFactorViscosity();
  }

  void print(std::ostream &clout) const override {
      using Align = io::text::Table::Align;
      clout << io::text::section("Adsorption Unit Converter") << std::endl;
      io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
      parameters.row("Resolution",                  "N",              this->getResolution());
      parameters.row("Lattice velocity",            "latticeU",       this->getCharLatticeVelocity());
      parameters.row("Lattice relaxation frequency","omega",          this->getLatticeRelaxationFrequency());
      parameters.row("Lattice relaxation time",     "tau",            this->getLatticeRelaxationTime());
      parameters.row("Char. length [m]",            "charL",          this->getCharPhysLength());
      parameters.row("Char. speed [m/s]",           "charU",          this->getCharPhysVelocity());
      parameters.row("Concentration [kg/m^d]",      "charC",          this->getPhysDensity());
      parameters.row("Diffusivity [m^2/s]",         "Diffusivity",    this->getPhysDiffusivity());
      parameters.row("Viscosity [m^2/s]",           "Viscosity",      this->getPhysViscosity());
      parameters.row("Schmidt number",              "Sc",             this->getSchmidtNumber());
      parameters.row("Reynolds number",             "reynoldsNumber", this->getReynoldsNumber());
      parameters.row("Mach number",                 "machNumber",     this->getMachNumber());
      parameters.row("Fourier number",              "fourierNumber",  this->getFourierNumber());
      parameters.render(clout);

      clout << std::endl;
      clout << io::text::section("Conversion Factors") << std::endl;
      io::text::Table factors({{"Factor"}, {"Symbol"}, {"Value", Align::Right}});
      factors.row("Voxel length [m]",       "physDeltaX",      this->getConversionFactorLength());
      factors.row("Time step [s]",          "physDeltaT",      this->getConversionFactorTime());
      factors.row("Velocity [m/s]",         "physVelocity",    this->getConversionFactorVelocity());
      factors.row("Density [kg/m^3]",       "physDensity",     this->getConversionFactorDensity());
      factors.row("Mass [kg]",              "physMass",        this->getConversionFactorMass());
      factors.row("Viscosity [m^2/s]",      "physViscosity",   this->getConversionFactorViscosity());
      factors.row("Diffusivity [m^2/s]",    "physDiffusivity", this->getConversionFactorDiffusivity());
      factors.row("Force [N]",              "physForce",       this->getConversionFactorForce());
      factors.row("Pressure [N/m^2]",       "physPressure",    this->getConversionFactorPressure());
      factors.row("Particle density factor","",                this->getConversionFactorParticleDensity());
      factors.render(clout);
  }

  using UnitConverter<T,DESCRIPTOR>::print;

};

template<typename T, typename DESCRIPTOR>
class AdsorptionConverterFromSchmidtNumberAndRelaxation : public AdsorptionConverter<T, DESCRIPTOR> {
 public:
  /** Calculate parameters from nondimensionalized values
   *
   * @param schmidtNumber
   * @param reynoldsNumber
   * @param relaxationTimeADV
   * @param charPhysLength
   * @param particleLength
   * @param resolution
   * @param charPhysVelocity
   * @param physDensity
   */
  constexpr AdsorptionConverterFromSchmidtNumberAndRelaxation(
      T schmidtNumber,
      T reynoldsNumber,
      T relaxationTimeADV,
      T charPhysLength,
      T particleLength,
      T resolution,
      T charPhysVelocity,
      T physDensity,
      T particleConcentration) : AdsorptionConverter<T, DESCRIPTOR>(charPhysLength / resolution,
                                                          (relaxationTimeADV - 0.5) / descriptors::invCs2<T,DESCRIPTOR>()
                                                              * util::pow((charPhysLength / resolution), 2)
                                                              / (charPhysVelocity * particleLength
                                                              / (schmidtNumber * reynoldsNumber)),
                                                          particleLength,
                                                          charPhysVelocity,
                                                          charPhysVelocity * particleLength
                                                              / (schmidtNumber * reynoldsNumber),
                                                          physDensity,
                                                          particleConcentration,
                                                          charPhysVelocity*particleLength/(reynoldsNumber)) {};

  using UnitConverter<T,DESCRIPTOR>::print;
};

}
#endif // OLB_APPS_FLORIAN_ADSORPTIONCONVERTER_H_
