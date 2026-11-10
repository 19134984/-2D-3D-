/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Stephan Simonis, Davide Dapelo,
 *                2024 Marc Heinzelmann
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

/*  Transient Longitudinal Mixing:
 *  The displacement of an initially homogeneous solute
 *  by the introduction of a solvent.
 *
 *  The numerical setup and analytical solution are taken from
 *  Long Ju, Chunhua Zhang, and Zhaoli Guo. “Local reactive boundary scheme for
 *  irregular geometries in lattice Boltzmann method”. In: International Journal of Heat
 *  and Mass Transfer 150 (2020), p. 119314. doi: 10.1016/j.ijheatmasstransfer.2020.119314.
 *
 *  Numerical solution can be calculated using two approaches.
 *  Error norms are calculated for three norms.
 */

#include "case.h"

int main(int argc, char *argv[]) {
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(32);
    myCaseParameters.set<ERROR_TIME>(1.);
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.04);
    myCaseParameters.set<PECLET>(1.);
    myCaseParameters.set<PHYS_DIFFUSIVITY>(0.01);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.53);
    myCaseParameters.set<PHYS_MAX_T>(1.01);
    myCaseParameters.set<C_EQ>(50.);
    myCaseParameters.set<PHYS_DENSITY>(1.);
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);

  clout <<"Errors with N="<<(myCaseParameters.get<parameters::RESOLUTION>())
        <<" after "<< myCaseParameters.get<parameters::ERROR_TIME>() << " s: "
        <<" log L1: " << myCaseParameters.get<parameters::L1_ERROR>()
        <<" log L2: "<< myCaseParameters.get<parameters::L2_ERROR>()
        <<" log Linf: "<< myCaseParameters.get<parameters::LINF_ERROR>()
        << std::endl;
}
