
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

/*  adePeriodic2d:
 *  The solution to a linear, scalar, two-dimensional advection-diffusion
 *  equation is approximated.
 *  The numerical setup and the analytical solution are taken from
 *  Simonis, S., Frank, M., and Krause M. J. Applied Mathematics Letters
 *  (2023) 135:108484, DOI: https://doi.org/10.1016/j.aml.2022.108484.
 *  Error norms are calculated for three subsequent resolutions and stored
 *  in the respective /tmp folders. A python script is provided to calculate
 *  the experimental order of convergence towards the analytical solution.
 */

#include "../case.h"

using namespace olb;
using namespace olb::graphics;
using namespace olb::names;


int main(int argc, char *argv[])
{
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<NUMBER_OF_RUNS>(3);
    myCaseParameters.set<RESOLUTION>(50);
    myCaseParameters.set<LATTICE_STAT_ITER_T>(50);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(2.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.);
    myCaseParameters.set<PHYS_DIFFUSIVITY>(0.05);
    myCaseParameters.set<PECLET>(100.);
    myCaseParameters.set<PULSE_DIFF_BOUND>(1e-1);
  }

  const MyCase::value_t N0 = myCaseParameters.get<parameters::RESOLUTION>();
  const MyCase::value_t statIter0 = myCaseParameters.get<parameters::LATTICE_STAT_ITER_T>();
  myCaseParameters.set<parameters::RESOLUTION_MAX>(
    util::pow(2, myCaseParameters.get<parameters::NUMBER_OF_RUNS>()-1) * myCaseParameters.get<parameters::RESOLUTION>()
  );

  CSV<MyCase::value_t> csvWriter("averageSimL2RelErr");

  for (std::size_t i = 0; i < myCaseParameters.get<parameters::NUMBER_OF_RUNS>(); ++i) {

    //Adjust Resolution and output interval
    myCaseParameters.set<parameters::RESOLUTION>(util::pow(2,i) * N0);
    myCaseParameters.set<parameters::LATTICE_STAT_ITER_T>(util::pow(4,i) * statIter0);


    Mesh mesh = createMesh(myCaseParameters);

    MyCase myCase(myCaseParameters, mesh);

    prepareGeometry(myCase);

    prepareLattice(myCase);

    setInitialValues(myCase);

    simulate(myCase);

    // this outputs into ./tmp/gnuplotData/data/averageSimL2RelErr
    singleton::directories().setOutputDir("./tmp/");
    csvWriter.writeDataFile(myCaseParameters.get<parameters::RESOLUTION>(),
                            myCaseParameters.get<parameters::L2_ERROR>(),
                            "averageSimL2RelErr");
  }
}
