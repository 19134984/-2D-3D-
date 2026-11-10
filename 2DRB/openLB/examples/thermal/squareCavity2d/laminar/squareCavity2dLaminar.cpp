/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Florian Kaiser
 *                2008 Orestis Malaspinas
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

/* squareCavity2dLaminar.cpp
 * The reference is the paper in "Gaedtke, M., Wachter, S., Raedle, M., Nirschl, H., & Krause, M. J. (2018).
 * Application of a lattice Boltzmann method combined with a Smagorinsky turbulence model to spatially resolved heat flux inside a refrigerated vehicle.
 * Computers & Mathematics with Applications, 76(10), 2315-2329."
 */

// natural convection of air in a square cavity in 2D

#include <olb.h>
#include "case.h"
#include "referenceData/referenceData.h"

void compareToStudy(MyCase& myCase)
{
  OstreamManager clout(std::cout,"compareToStudy");

  using T                         = MyCase::value_t_of<NavierStokes>;

  auto& parameters                = myCase.getParameters();
  auto& converter                 = myCase.getLattice(NavierStokes{}).getUnitConverter();

  const std::size_t Ra            = parameters.get<parameters::RAYLEIGH>();
  const T physThermalDiffusivity  = converter.getPhysThermalDiffusivity();
  const T charL                   = converter.getCharPhysLength();

  ReferenceData<T> ref{};

  auto simValues  = parameters.get<parameters::OUTPUT_VALUE_WRAPPER>();

  if (ref.hasRayleigh(Ra) && singleton::mpi().isMainProcessor()) {
    T xVelComp      = ref.getUx_max(Ra);
    T yVelComp      = ref.getUy_max(Ra);
    T xCoordComp    = ref.getX_max(Ra);
    T yCoordComp    = ref.getY_max(Ra);
    T nusseltComp   = ref.getNusselt(Ra);

    clout << "Comparison against De Vahl Davis (1983):" << std::endl;
    clout << "xVelocity in yDir="   <<  simValues[0] / physThermalDiffusivity * charL   << "; error(rel)=" << (T) util::fabs((xVelComp - simValues[0] / physThermalDiffusivity * charL) / xVelComp)                 << std::endl;
    clout << "yVelocity in xDir="   <<  simValues[1] / physThermalDiffusivity * charL   << "; error(rel)=" << (T) util::fabs((yVelComp - simValues[1] / physThermalDiffusivity * charL) / yVelComp)                 << std::endl;
    clout << "yMaxVel / xMaxVel="   <<  simValues[1] / simValues[0]                     << "; error(rel)=" << (T) util::fabs((simValues[1] / simValues[0]) - (yVelComp / xVelComp)) / (yVelComp / xVelComp)         << std::endl;
    clout << "yCoord of xMaxVel="   <<  simValues[2] / charL                            << "; error(rel)=" << (T) util::fabs((yCoordComp - simValues[2] / charL) / yCoordComp)                                      << std::endl;
    clout << "xCoord of yMaxVel="   <<  simValues[3] / charL                            << "; error(rel)=" << (T) util::fabs((xCoordComp - simValues[3] / charL) / xCoordComp)                                      << std::endl;
    clout << "Nusselt="             <<  simValues[4]                                    << "; error(rel)=" << (T) util::fabs((nusseltComp - simValues[4]) / nusseltComp)                                            << std::endl;
    std::fstream fs;
    fs.open("output.txt",
      std::fstream::in | std::fstream::out | std::fstream::app);
    fs    << "xVelocity in yDir="   <<  simValues[0] / physThermalDiffusivity * charL   << "; error(rel)=" << (T) util::fabs((xVelComp - simValues[0] / physThermalDiffusivity * charL) / xVelComp)                 << std::endl;
    fs    << "yVelocity in xDir="   <<  simValues[1] / physThermalDiffusivity * charL   << "; error(rel)=" << (T) util::fabs((yVelComp - simValues[1] / physThermalDiffusivity * charL) / yVelComp)                 << std::endl;
    fs    << "yMaxVel / xMaxVel="   <<  simValues[1] / simValues[0]                     << "; error(rel)=" << (T) util::fabs((simValues[1] / simValues[0]) - (yVelComp / xVelComp)) / (yVelComp / xVelComp)         << std::endl;
    fs    << "yCoord of xMaxVel="   <<  simValues[2] / charL                            << "; error(rel)=" << (T) util::fabs((yCoordComp - simValues[2] / charL) / yCoordComp)                                      << std::endl;
    fs    << "xCoord of yMaxVel="   <<  simValues[3] / charL                            << "; error(rel)=" << (T) util::fabs((xCoordComp - simValues[3] / charL) / xCoordComp)                                      << std::endl;
    fs    << "Nusselt="             <<  simValues[4]                                    << "; error(rel)=" << (T) util::fabs((nusseltComp - simValues[4]) / nusseltComp)                                            << std::endl;
    fs.close();
  } else {
    clout << "Comparison against De Vahl Davis (1983) not possible. Received Ra = "     << Ra << std::endl;
    clout << "xVelocity in yDir="   <<  simValues[0] / physThermalDiffusivity * charL         << std::endl;
    clout << "yVelocity in xDir="   <<  simValues[1] / physThermalDiffusivity * charL         << std::endl;
    clout << "yMaxVel / xMaxVel="   <<  simValues[1] / simValues[0]                           << std::endl;
    clout << "yCoord of xMaxVel="   <<  simValues[2] / charL                                  << std::endl;
    clout << "xCoord of yMaxVel="   <<  simValues[3] / charL                                  << std::endl;
    clout << "Nusselt="             <<  simValues[4]                                          << std::endl;
  }
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  using T = MyCase::value_t_of<NavierStokes>;

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_MAX_T                 >(1e4);
    myCaseParameters.set<RESOLUTION                 >(64);
    myCaseParameters.set<GRAVITATIONAL_ACC          >(9.81);
    myCaseParameters.set<PHYS_CHAR_VELOCITY         >(1);
    myCaseParameters.set<DOMAIN_ORIGIN              >({0., 0.});
    myCaseParameters.set<LATTICE_RELAXATION_TIME    >(0.65);
    myCaseParameters.set<CONVERGENCE_PRECISION      >(1e-4);
    myCaseParameters.set<PHYS_THERMAL_EXPANSION     >(3.41e-3);
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY  >(25.684e-3);
    myCaseParameters.set<PHYS_HEAT_CAPACITY         >(1.01309e3);
    myCaseParameters.set<PHYS_CHAR_DENSITY          >(1.19);
    myCaseParameters.set<PHYS_KINEMATIC_VISCOSITY   >(1.5126e-5);

    myCaseParameters.set<PRANDTL                    >(0.71);
    myCaseParameters.set<RAYLEIGH                   >(1e3);

    myCaseParameters.set<T_HOT >(285.15);
    myCaseParameters.set<T_COLD>(275.15);
    myCaseParameters.set<T_MEAN>(
        (myCaseParameters.get<T_HOT>() + myCaseParameters.get<T_COLD>()) / 2.0
    );

    myCaseParameters.set<CONV_ITER                  >(1000);
    myCaseParameters.set<PHYS_VTK_ITER_T            >(5);
    myCaseParameters.set<PHYS_STAT_ITER_T           >(1);
    myCaseParameters.set<CONVERGENCE_PRECISION      >(1e-5);
  }

  const int N              = myCaseParameters.get<parameters::RESOLUTION               >();
  const std::size_t Ra     = myCaseParameters.get<parameters::RAYLEIGH                 >();
  const T physViscosity    = myCaseParameters.get<parameters::PHYS_KINEMATIC_VISCOSITY >();
  const T Pr               = myCaseParameters.get<parameters::PRANDTL                  >();
  const T g                = myCaseParameters.get<parameters::GRAVITATIONAL_ACC        >();
  const T Thot             = myCaseParameters.get<parameters::T_HOT                    >();
  const T Tcold            = myCaseParameters.get<parameters::T_COLD                   >();
  const T thermalExpansion = myCaseParameters.get<parameters::PHYS_THERMAL_EXPANSION   >();

  // lx from Rayleigh number
  const T lx               = util::pow(Ra * physViscosity * physViscosity / (Pr * g * (Thot - Tcold) * thermalExpansion),
                                      (T) 1 / 3);  // length of the square
  myCaseParameters.set<parameters::PHYS_CHAR_LENGTH>(lx);
  myCaseParameters.set<parameters::PHYS_DELTA_X>(lx / N);
  myCaseParameters.set<parameters::DOMAIN_EXTENT>([](MyCase::ParametersD& p) -> Vector<T, MyCase::d> {
    return {
      p.get<parameters::PHYS_CHAR_LENGTH>(), p.get<parameters::PHYS_CHAR_LENGTH>()
    };
  });

  const T charU = physViscosity / lx * sqrt((T) Ra / Pr);
  myCaseParameters.set<parameters::PHYS_CHAR_VELOCITY>(charU);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

  compareToStudy(myCase);
}
