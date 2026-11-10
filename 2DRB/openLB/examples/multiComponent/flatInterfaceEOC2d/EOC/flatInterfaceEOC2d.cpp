/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
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

/** flatInterface2d.cpp
 * In this example a flat interface test is performed for several multi-
 * phase models such as incompressible Allen-Cahn or well-balanced Cahn-
 * Hilliard. A rectangular domain of fluid 2 is immersed in fluid 1.
 * A diffusive interface forms with the profile of a hyperbolic tangent
 * whose accuracy is measured for multiple resolutions in order to test
 * the models experimental order of convergence. The equilibrium pressure
 * is also investigated in a similar way. The grid refinement can be
 * performed with either a constant or a decreasing Cahn number. Multiple
 * different error norms can then be written into an output EOC.dat file.
 *
 * This example shows the simplest application of the hybrid (and local)
 * phase field Allen-Cahn model with periodic boundaries, based on:
 *
 * Liu, Xi, Zhenhua Chai, and Baochang Shi. "Improved hybrid Allen-Cahn
 * phase-field-based lattice Boltzmann method for incompressible two-phase
 * flows." Physical Review E 107.3 (2023): 035308.
 *
 * It also shows the same application of the well-balanced Cahn-Hilliard
 * model that has no spurious currents (machine accuracy), from:
 *
 * Ju, Long, et al. "A well-balanced lattice Boltzmann model for binary
 * fluids based on the incompressible phase-field theory." arXiv preprint
 * arXiv:2311.10827 (2023).
 */

// Choose model from loading three different cases
#include "../flatInterface2d/flatInterfaceConservativePhaseField2d.h"
// #include "../flatInterface2d/flatInterfaceHybridAllenCahn2d.h"
// #include "../flatInterface2d/flatInterfaceCahnHilliard2d.h"

const bool Cahn_const = true;  // scaling with const Cahn number?
const int Ny = 100;            // initial resolution
const double w = 5.;           // initial interface thickness
const double sigma = 0.01;     // initial surface tension

using T = MyCase::value_t;

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );

  std::ofstream outfile;
  outfile.open ("EOC.dat");

  // Arrays to store data for plotting
  std::vector<T> resolutions;
  std::vector<T> orderParamL1Errors, orderParamL2Errors, orderParamLinfErrors;
  std::vector<T> pressureL1Errors, pressureL2Errors, pressureLinfErrors;
  std::vector<T> velocityLinfErrors;

  //Update the geometry, lattice, and other dependent variables for each Ny
  for (size_t i=0; i<4; ++i) {
    singleton::directories().setOutputDir("./tmp/sub" + std::to_string(i) + "/");
    /// === Step 2: Set Parameters ===
    MyCase::ParametersD myCaseParameters;
    {
      using namespace olb::parameters;
      myCaseParameters.set<DOMAIN_EXTENT          >({2e-6, 100e-6});    // domain size [physical units]
      myCaseParameters.set<RESOLUTION             >(Ny*util::pow(2,i)); // resolution y [lattice units]
      myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8);       // tau liquid [lattice units]
      myCaseParameters.set<RELAXATION_TIME_PF     >(0.8);            // tau mobility [lattice units]
      myCaseParameters.set<REYNOLDS               >(0.);             // Reynolds number
      myCaseParameters.set<RHO_VAPOR              >(1.);             // physDensity gas/vapor [physical units]
      myCaseParameters.set<RHO_LIQUID             >(1000.);          // physDensity liquid [physical units]
      myCaseParameters.set<NU_VAPOR               >(9e-7);           // physViscosity gas/vapor [physical units]
      myCaseParameters.set<NU_LIQUID              >(9e-6);           // physViscosity liquid [physical units]
      myCaseParameters.set<PHYS_CHAR_PRESSURE     >(1e5);            // physPressure [physical units]
      myCaseParameters.set<SURFACE_TENSION        >(
        sigma/util::pow(2,i));                                  // surface tension [lattice units]
      myCaseParameters.set<parameters::INTERFACE_WIDTH>(
        w*util::pow(2,i*(Cahn_const)));                         // interface width [lattice units]
      myCaseParameters.set<C_RHO             >(1000.);          // conversion factor density [physical units]
    }

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
    simulate(myCase, Cahn_const);

    // Store data for plotting
    T resolution = Ny*util::pow(2,i);
    resolutions.push_back(resolution);

    if (Cahn_const) {
      orderParamL1Errors.push_back(orderParameterL1RelError);
      orderParamL2Errors.push_back(orderParameterL2RelError);
      orderParamLinfErrors.push_back(orderParameterLinfRelError);
      outfile << resolution << "," << orderParameterL2RelError << "," << pressureL2RelError << "," << velocityLinfAbsError << "\n";
    }
    else {
      orderParamL1Errors.push_back(orderParameterL1RelError);
      orderParamL2Errors.push_back(orderParameterL2RelError);
      orderParamLinfErrors.push_back(orderParameterLinfRelError);
      outfile << resolution << "," << orderParameterL1RelError << "," << pressureL2RelError << "," << velocityLinfAbsError << "\n";
    }

    pressureL1Errors.push_back(pressureL1RelError);
    pressureL2Errors.push_back(pressureL2RelError);
    pressureLinfErrors.push_back(pressureLinfRelError);
    velocityLinfErrors.push_back(velocityLinfAbsError);
  }
  outfile.close();

  singleton::directories().setOutputDir("./tmp/");

  // --- helpers -------------------------------------------------------------
  // Least-squares slope of log10(error) vs log10(resolution), sign-flipped so
  // that a positive number means "error falls with refinement". This is the
  // overall rate across all resolutions, not a pairwise ratio between
  // neighbouring grids, so a single outlier can't dominate it.
  auto convergenceRate = [](const std::vector<T>& x, const std::vector<T>& y) -> T {
    T sx = 0., sy = 0., sxx = 0., sxy = 0.;
    int n = 0;
    for (size_t i = 0; i < x.size() && i < y.size(); ++i) {
      if (x[i] <= T(0) || y[i] <= T(0) || !std::isfinite(y[i])) {
        continue;   // skip zeros / round-off floors: log10 would blow up
      }
      const T lx = util::log10(x[i]);
      const T ly = util::log10(y[i]);
      sx += lx; sy += ly; sxx += lx*lx; sxy += lx*ly; ++n;
    }
    if (n < 2) { return T(0); }
    const T denom = T(n)*sxx - sx*sx;
    if (util::nearZero(denom)) { return T(0); }
    return -(T(n)*sxy - sx*sy) / denom;
  };

  // "Pressure L1" -> "Pressure L1, Convergence = 1.98"
  auto legendName = [](const std::string& tag, T rate) {
    std::ostringstream oss;
    oss << tag << ", Convergence = " << std::fixed << std::setprecision(2)
        << BaseType<T>(rate);
    return oss.str();
  };

  // --- one plot holding every dataset --------------------------------------
  Gnuplot<T> gplotEOC("eoc_error", false,
    "set terminal png size 720, 720 font 'Helvetica,10'",
    Gnuplot<T>::LOGLOG, Gnuplot<T>::LINREG);

  const std::vector<std::string> names = {
    legendName("Order Parameter L1",   convergenceRate(resolutions, orderParamL1Errors)),
    legendName("Order Parameter L2",   convergenceRate(resolutions, orderParamL2Errors)),
    legendName("Order Parameter Linf", convergenceRate(resolutions, orderParamLinfErrors)),
    legendName("Pressure L1",          convergenceRate(resolutions, pressureL1Errors)),
    legendName("Pressure L2",          convergenceRate(resolutions, pressureL2Errors)),
    legendName("Pressure Linf",        convergenceRate(resolutions, pressureLinfErrors))
  };

  gplotEOC.setLabel("Resolution", "Error (order parameter / pressure relative, velocity absolute)");

  const std::vector<char> markers(names.size(), 'p');

  for (size_t i = 0; i < resolutions.size(); ++i) {
    gplotEOC.setData(
      T(resolutions[i]),
      { orderParamL1Errors[i], orderParamL2Errors[i], orderParamLinfErrors[i],
        pressureL1Errors[i],   pressureL2Errors[i],   pressureLinfErrors[i] },
      (i == 0) ? names : std::vector<std::string>(names.size(), ""),
      "bottom left",
      markers
    );
  }
  gplotEOC.writePNG();

  // Echo the rates to the terminal as well
  {
    OstreamManager clout(std::cout, "EOC");
    for (const auto& n : names) { clout << n << std::endl; }
  }
}
