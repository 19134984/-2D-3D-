/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017-2018 Max Gaedtke, Albert Mink, Davide Dapelo
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


#ifndef PL_UNITCONVERTER_HH
#define PL_UNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "io/textTable.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void PowerLawUnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const
{
  using Align = io::text::Table::Align;
  clout << io::text::section("Power-Law Unit Converter") << std::endl;
  io::text::Table parameters({{"Parameter"}, {"Symbol"}, {"Value", Align::Right}});
  parameters.row("Resolution",                        "N",              this->getResolution());
  parameters.row("Lattice velocity",                  "latticeU",       this->getCharLatticeVelocity());
  parameters.row("Lattice relaxation frequency",      "omega",          this->getLatticeRelaxationFrequency());
  parameters.row("Lattice relaxation time",           "tau",            this->getLatticeRelaxationTime());
  parameters.row("Char. length [m]",                  "charL",          this->getCharPhysLength());
  parameters.row("Char. speed [m/s]",                 "charU",          this->getCharPhysVelocity());
  parameters.row("Kinematic viscosity [m^2/s]",       "charNu",         this->getPhysViscosity());
  parameters.row("Consistency coeff. [m^2 s^(n-2)]",  "charM",          this->getPhysConsistencyCoeff());
  parameters.row("Power-law index",                   "n",              this->getPowerLawIndex());
  parameters.row("Density [kg/m^d]",                  "charRho",        this->getPhysDensity());
  parameters.row("Char. pressure [N/m^2]",            "charPressure",   this->getCharPhysPressure());
  parameters.row("Reynolds number",                   "reynoldsNumber", this->getReynoldsNumber());
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

template<typename T, typename DESCRIPTOR>
PowerLawUnitConverter<T, DESCRIPTOR>* createPowerLawUnitConverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createUnitConverter");
  params.setWarningsOn(false);

  T physDeltaX;
  T physDeltaT;

  T charPhysLength;
  T charPhysVelocity;
  T physViscosity;
  T physConsistencyCoeff;
  T powerLawIndex;
  T physDensity;
  T charPhysPressure = 0;

  int resolution;
  T latticeRelaxationTime;
  T charLatticeVelocity;

  int counter = 0; // counting the number of Discretization parameters and returning a warning if more than 2 are provided


  // params[parameter].read(value) sets the value or returns false if the parameter can not be found
  params["Application"]["PhysParameters"]["CharPhysLength"].read(charPhysLength);
  params["Application"]["PhysParameters"]["CharPhysVelocity"].read(charPhysVelocity);
  params["Application"]["PhysParameters"]["PhysConsistencyCoeff"].read(physConsistencyCoeff);
  params["Application"]["PhysParameters"]["powerLawIndex"].read(powerLawIndex);
  params["Application"]["PhysParameters"]["PhysDensity"].read(physDensity);
  params["Application"]["PhysParameters"]["CharPhysPressure"].read(charPhysPressure);

  physViscosity = physConsistencyCoeff * util::pow(charPhysVelocity / (2*charPhysLength), powerLawIndex-1);

  std::vector<std::string> discretizationParam = {"PhysDeltaX", "Resolution",
    "CharLatticeVelocity", "PhysDeltaT", "LatticeRelaxationTime"};

  for(int i = 0; i<discretizationParam.size(); i++){
    std::string test;
    if(params["Application"]["Discretization"][discretizationParam[i]].read(test,false)){
      counter++;
    }
  }
  if(counter>2){
    clout << "WARNING: More than 2 discretization parameters provided" << std::endl;
  }

  if (!params["Application"]["Discretization"]["PhysDeltaX"].read(physDeltaX,false)) {
    if (!params["Application"]["Discretization"]["Resolution"].read<int>(resolution,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaX, resolution or charLatticeVelocity
        clout << "Error: Have not found PhysDeltaX, Resolution or CharLatticeVelocity in XML file."
              << std::endl;
        exit (1);
      }
      else {
        // found charLatticeVelocity
        if (params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
          physDeltaX = charPhysVelocity / charLatticeVelocity * physDeltaT;
        }
        else if (params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          physDeltaX = physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5);
        }
      }
    }
    else {
      // found resolution
      physDeltaX = charPhysLength / resolution;
    }
  }
  // found physDeltaX
  if (!params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
    if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaT, latticeRelaxationTime and charLatticeVelocity
        clout << "Error: Have not found PhysDeltaT, LatticeRelaxationTime or CharLatticeVelocity in XML file."
              << std::endl;
        exit (1);
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

  return new PowerLawUnitConverter<T, DESCRIPTOR>(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physConsistencyCoeff, powerLawIndex, physDensity, charPhysPressure);
}

}  // namespace olb

#endif
