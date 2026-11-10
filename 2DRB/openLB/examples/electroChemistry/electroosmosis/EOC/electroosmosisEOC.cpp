/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev
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

#include "../electroosmosis.h"
using T = MyCase::value_t;

static Gnuplot<T> gplot(
  "ElectroosmosisEOC",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<T>::LOGLOG,
  Gnuplot<T>::LINREG);

namespace olb::parameters {


}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION_START>(20);
    myCaseParameters.set<RESOLUTION_END>(80);
    myCaseParameters.set<E_FIELD>(250);
    myCaseParameters.set<POISSON_RELAXATION_TIME>(0.9);
    myCaseParameters.set<IONS_RELAXATION_TIME>(0.7);
    myCaseParameters.set<NSE_RELAXATION_TIME>(0.7);
    myCaseParameters.set<NSE_PHYS_CHAR_VELOCITY>(0.5);
    myCaseParameters.set<NPE_PHYS_CHAR_VELOCITY>(0.5);
    myCaseParameters.set<POISSON_PHYS_CHAR_VELOCITY>(0.001);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1.e-6);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1000.);
    myCaseParameters.set<PHYS_MAX_T>(10);
    myCaseParameters.set<RESIDUUM>(1.e-9);
    myCaseParameters.set<DIFFUSION>(1.e-8);
    myCaseParameters.set<VALENCE>(1);
    myCaseParameters.set<parameters::TEMPERATURE>(293.15);
    myCaseParameters.set<DIELECTRIC_CONSTANT>(6.95e-10);
    myCaseParameters.set<C_0>(0.01);
    myCaseParameters.set<PSI_BC>(-0.02);
    myCaseParameters.set<DEBYE>([] (MyCase::ParametersD& p) {
      return util::debyeLength<T>(p.get<DIELECTRIC_CONSTANT>(),
                                  p.get<VALENCE>(),
                                  p.get<parameters::TEMPERATURE>(),
                                  p.get<C_0>());
    });
    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) {
      return Vector<T,3>(13./4.*p.get<DEBYE>(), 13.*p.get<DEBYE>(),0.);
    });
    myCaseParameters.set<CB_CATION>([] (MyCase::ParametersD& p) {
      return p.get<C_0>()*util::exp(
             -physConstants::elementaryCharge<T>()*p.get<VALENCE>()
             *p.get<PSI_BC>()/physConstants::boltzmannConstant<T>()
             /p.get<parameters::TEMPERATURE>());
    });
    myCaseParameters.set<CB_ANION>([] (MyCase::ParametersD& p) {
      return p.get<C_0>()*util::exp(
             physConstants::elementaryCharge<T>()*p.get<VALENCE>()
             *p.get<PSI_BC>()/physConstants::boltzmannConstant<T>()
             /p.get<parameters::TEMPERATURE>());
    });
    myCaseParameters.set<ERROR_NORMS_PSI>({0.,0.,0.});
    myCaseParameters.set<ERROR_NORMS_CATION>({0.,0.,0.});
    myCaseParameters.set<ERROR_NORMS_ANION>({0.,0.,0.});
    myCaseParameters.set<ERROR_NORMS_NSE>({0.,0.,0.});
    myCaseParameters.set<CONVERGED>(false);
  }

  gplot.setLabel("Resolution test", "average Error");

  for(int simuN =  myCaseParameters.get<parameters::RESOLUTION_START>();
      simuN <=  myCaseParameters.get<parameters::RESOLUTION_END>();
      simuN *= 2) {
    myCaseParameters.set<parameters::RESOLUTION>(simuN);

    /// === Step 3: Create Mesh ===
    Mesh mesh = createMesh(myCaseParameters);

    /// === Step 4: Create Case ===
    MyCase myCase(myCaseParameters, mesh);

    /// === Step 5: Prepare Geometry ===
    prepareGeometry(myCase);

    /// === Step 6: Prepare Lattice ===
    prepareLatticePoisson(myCase);
    prepareLatticeNernstPlanck<0>(myCase);
    prepareLatticeNernstPlanck<1>(myCase);
    prepareLatticeNSE(myCase);

    prepareLatticeCoupling(myCase);

    /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
    setInitialValuesPoisson(myCase);
    setInitialValuesCation(myCase);
    setInitialValuesAnion(myCase);
    setInitialValuesNSE(myCase);

    /// === Step 8: Simulate ===
    simulate(myCase);

    auto errorsPsi = myCaseParameters.get<parameters::ERROR_NORMS_PSI>();
    auto errorsCation = myCaseParameters.get<parameters::ERROR_NORMS_CATION>();
    auto errorsAnion = myCaseParameters.get<parameters::ERROR_NORMS_ANION>();
    auto errorsNSE = myCaseParameters.get<parameters::ERROR_NORMS_NSE>();
    gplot.setData (
          MyCase::value_t(simuN),
          { errorsPsi[0], errorsPsi[1], errorsPsi[2],
            errorsCation[0], errorsCation[1], errorsCation[2],
            errorsAnion[0], errorsAnion[1], errorsAnion[2],
            errorsNSE[0], errorsNSE[1], errorsNSE[2]},
          { "psi L1 Rel Error","psi L2 Rel Error",
            "psi Linf Rel error","conc L1 Rel Error","conc L2 Rel Error",
            "conc Linf Rel error","conc2 L1 Rel Error","conc2 L2 Rel Error",
            "conc2 Linf Rel error","vel L1 Rel Error","vel  L2 Rel Error",
            "vel  Linf Rel error"},
          "top right",
          { 'p','p','p','p','p','p','p','p','p','p','p','p' } );
  }
  gplot.writePNG();
}
