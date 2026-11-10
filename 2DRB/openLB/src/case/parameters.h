/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#ifndef CASE_PARAMETERS_H
#define CASE_PARAMETERS_H

#include "io/cliReader.h"
#include "core/meta.h"
#include "core/parameters.h"

namespace olb {

namespace parameters {

// Characteristic numbers
struct REYNOLDS     : public descriptors::FIELD_BASE<1> { };
struct RAYLEIGH     : public descriptors::FIELD_BASE<1> { };
struct PRANDTL      : public descriptors::FIELD_BASE<1> { };
struct TAYLOR       : public descriptors::FIELD_BASE<1> { };
struct SMAGORINSKY  : public descriptors::FIELD_BASE<1> { };
struct SCHMIDT      : public descriptors::FIELD_BASE<1> { };
struct PECLET       : public descriptors::FIELD_BASE<1> { };
struct KNUDSEN      : public descriptors::FIELD_BASE<1> { };
struct STOKES       : public descriptors::FIELD_BASE<1> { };
struct NUSSELT      : public descriptors::FIELD_BASE<1> { };
struct PRANDTL_TURB : public descriptors::FIELD_BASE<1> { };
struct ARCHIMEDES   : public descriptors::FIELD_BASE<1> { };
struct STEFAN       : public descriptors::FIELD_BASE<1> { };
struct MACH         : public descriptors::FIELD_BASE<1> { };

// Converter-related parameters
struct RESOLUTION : public descriptors::TYPED_FIELD_BASE<int,1> {
  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct TIME_RESOLUTION : public descriptors::TYPED_FIELD_BASE<int,1> { };

struct PHYS_DELTA_X : public descriptors::FIELD_BASE<1> { };
struct PHYS_DELTA_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_LENGTH : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_TIME : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_VISCOSITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_DENSITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_PRESSURE : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_TEMPERATURE : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_HIGH_TEMPERATURE : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_LOW_TEMPERATURE : public descriptors::FIELD_BASE<1> { };
struct PHYS_THERMAL_CONDUCTIVITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_HEAT_CAPACITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_THERMAL_DIFFUSIVITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_THERMAL_EXPANSION : public descriptors::FIELD_BASE<1> { };

struct VOLUME : public descriptors::FIELD_BASE<1> { };

struct LATTICE_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct LATTICE_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
struct LATTICE_RELAXATION_FREQUENCY : public descriptors::FIELD_BASE<1> { };

struct CFL : public descriptors::FIELD_BASE<1> {};

struct TIME_STEPS: public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };

//Geometry
struct DOMAIN_EXTENT : public descriptors::FIELD_BASE<0,1> { };
struct DOMAIN_ORIGIN : public descriptors::FIELD_BASE<0,1> { };
struct ORIGIN : public descriptors::FIELD_BASE<0,1> { };
struct DOMAIN_L : public descriptors::FIELD_BASE<1> { };
struct RADIUS : public descriptors::FIELD_BASE<1> { };
struct CENTER : public descriptors::FIELD_BASE<0,1> { };
struct ANGLE : public descriptors::FIELD_BASE<1> { };

//BC/IC/Constraint Flow
struct BULK_VELOCITY     : public descriptors::FIELD_BASE<1> { };
struct PHYS_MAX_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct PRESSURE          : public descriptors::FIELD_BASE<1> { };
struct AVERAGE_PRESSURE  : public descriptors::FIELD_BASE<1> { };
struct PRESSURE_DROP     : public descriptors::FIELD_BASE<1> { };
struct SHEAR_RATE        : public descriptors::FIELD_BASE<1> { };
struct STRAIN_RATE       : public descriptors::FIELD_BASE<1> { };

struct COEFFICIENT_OF_DRAG : public descriptors::FIELD_BASE<1> { };

//BC/IC/Constraints Heat
struct TEMPERATURE : public descriptors::FIELD_BASE<1> { };
struct T_HOT  : public descriptors::FIELD_BASE<1> { };
struct T_COLD : public descriptors::FIELD_BASE<1> { };
struct T_MEAN   : public descriptors::FIELD_BASE<1> { };
struct T_PERTURBATION : public descriptors::FIELD_BASE<1> { };

struct RHO_1  : public descriptors::FIELD_BASE<1> { };
struct RHO_2  : public descriptors::FIELD_BASE<1> { };

struct PHYS_CP : public descriptors::FIELD_BASE<1> { };
struct PHYS_LAMBDA : public descriptors::FIELD_BASE<1> { };

struct PHYS_DENSITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_KINEMATIC_VISCOSITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_DYNAMIC_VISCOSITY : public descriptors::FIELD_BASE<1> { };
struct PHYS_DIFFUSIVITY : public descriptors::FIELD_BASE<1> { };

// Solids (NCE)
struct YOUNGS_MODULUS : public descriptors::FIELD_BASE<1> { };
struct POISSON_RATIO : public descriptors::FIELD_BASE<1> { };
struct PHYS_CHAR_DISPLACEMENT : public descriptors::FIELD_BASE<1> { };
struct OMEGA_SOLID : public descriptors::FIELD_BASE<0,0,1> { };
struct SOLID_CASE             : public descriptors::TYPED_FIELD_BASE<std::string, 1> {};
struct DAMPING_FACTOR         : public descriptors::FIELD_BASE<1> {};

// Multicomponent
struct NU : public descriptors::FIELD_BASE<1> { };
struct NU_VAPOR        : public descriptors::FIELD_BASE<1> { };
struct NU_LIQUID       : public descriptors::FIELD_BASE<1> { };
struct RHO_VAPOR       : public descriptors::FIELD_BASE<1> { };
struct RHO_VAPOR_ANALYTICAL : public descriptors::FIELD_BASE<1> { };
struct RHO_LIQUID      : public descriptors::FIELD_BASE<1> { };
struct RHO_LIQUID_ANALYTICAL : public descriptors::FIELD_BASE<1> { };
struct THETA           : public descriptors::FIELD_BASE<1> { };
struct SURFACE_TENSION : public descriptors::FIELD_BASE<1> { };
struct INTERFACE_WIDTH : public descriptors::FIELD_BASE<1> { };
struct INTERFACE_VELOCITY_ANALYTICAL : public descriptors::FIELD_BASE<1> { };
struct LIQUID_PHASE_LENGTH : public descriptors::FIELD_BASE<1> { };
struct THICKNESS : public descriptors::FIELD_BASE<1> { };
struct BOUNDARY_CHEMICAL_POTENTIAL : public descriptors::FIELD_BASE<1> { };

struct COUPLING_G       : public descriptors::FIELD_BASE<1> { };
struct NOISE            : public descriptors::FIELD_BASE<1> { };
struct FORCE            : public descriptors::FIELD_BASE<1> { };
struct ZERO             : public descriptors::FIELD_BASE<1> { };

// Multiphase parameters
struct PHYS_EOS_A : public descriptors::FIELD_BASE<1> {};
struct PHYS_EOS_B : public descriptors::FIELD_BASE<1> {};
struct PHYS_MOLAR_MASS : public descriptors::FIELD_BASE<1> {};
struct LATTICE_EOS_A : public descriptors::FIELD_BASE<1> {};

// Power law parameters
struct PHYS_CONSISTENCY_COEFF : public descriptors::FIELD_BASE<1> {};
struct POWER_LAW_INDEX : public descriptors::FIELD_BASE<1> {};

// Radiative parameters
struct PHYS_ABSORPTION : public descriptors::FIELD_BASE<1> {};
struct PHYS_SCATTERING : public descriptors::FIELD_BASE<1> {};
struct REFRACTIVE_RELATIVE : public descriptors::FIELD_BASE<1> {};

struct VTK_ENABLED     : public descriptors::TYPED_FIELD_BASE<bool,1> { };
struct GNUPLOT_ENABLED : public descriptors::TYPED_FIELD_BASE<bool,1> { };

// Meta/IO/Time
struct PHYS_BOUNDARY_VALUE_UPDATE_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_START_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_UPDATE_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_INTERVAL_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_MAX_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_SAVE_ITER : public descriptors::FIELD_BASE<1> { };
struct PHYS_VTK_ITER_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_STAT_ITER_T       : public descriptors::FIELD_BASE<1> { };
struct PHYS_IMAGE_ITER_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_GNUPLOT_ITER_T : public descriptors::FIELD_BASE<1> { };
struct PHYS_START_PERIOD : public descriptors::FIELD_BASE<1> { };
struct PHYS_TIMER_ITER_T : public descriptors::FIELD_BASE<1>{};
struct PHYS_ITER_T : public descriptors::FIELD_BASE<1>{};

struct LATTICE_INTERVAL_T : public descriptors::FIELD_BASE<1> { };
struct LATTICE_MAX_T : public descriptors::TYPED_FIELD_BASE<int, 1> { };
struct LATTICE_STAT_ITER_T    : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct LATTICE_VTK_ITER_T     : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct LATTICE_LOG_ITER_T     : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct LATTICE_ITER_T : public descriptors::TYPED_FIELD_BASE<int,1>{};

struct RAMP_UP_T : public descriptors::FIELD_BASE<1> { };

struct STL_SIZE : public descriptors::FIELD_BASE<1> { };
struct OUTPUT_DIRECTORY : public descriptors::TYPED_FIELD_BASE<std::string, 1> { };

// Meta/EOC
struct CONVERGED : public descriptors::TYPED_FIELD_BASE<bool,1> { };
struct CONVERGENCE_PRECISION : public descriptors::FIELD_BASE<1> { };
struct EPSILON : public descriptors::FIELD_BASE<1> { };
struct RESIDUUM : public descriptors::FIELD_BASE<1> { };
struct CONV_ITER : public descriptors::FIELD_BASE<1> { };
struct PHYS_CONVERGENCE_CHECK_T : public descriptors::FIELD_BASE<1> { };
struct NUMBER_OF_RUNS : public descriptors::TYPED_FIELD_BASE<int, 1> { };
struct RESOLUTION_START : public descriptors::TYPED_FIELD_BASE<int, 1> { };
struct RESOLUTION_END : public descriptors::TYPED_FIELD_BASE<int, 1> { };
struct RESOLUTION_STEP : public descriptors::TYPED_FIELD_BASE<int, 1> { };
struct RESOLUTION_MAX : public descriptors::TYPED_FIELD_BASE<int, 1> { };

// Statistics
struct L1_ERROR : public descriptors::FIELD_BASE<1> { };
struct L2_ERROR : public descriptors::FIELD_BASE<1> { };
struct LINF_ERROR : public descriptors::FIELD_BASE<1> { };
struct ERROR : public descriptors::FIELD_BASE<1> { };

// Gravity
struct GRAVITY : public descriptors::FIELD_BASE<0,1> { };
struct GRAVITATIONAL_ACC  : public descriptors::FIELD_BASE<1> { };
struct INITIAL_FALLING_SPEED : public descriptors::FIELD_BASE<0,1> { };

// Turbulence Parameters
struct TURBULENCE_INTENSITY    : public descriptors::FIELD_BASE<1> { };
struct TURBULENCE_N_SEEDS      : public descriptors::FIELD_BASE<1> { };
struct TURBULENCE_SIGMA        : public descriptors::FIELD_BASE<1> { };

//Radiative Transport
struct ABSORPTION : public descriptors::FIELD_BASE<1> {};
struct SCATTERING : public descriptors::FIELD_BASE<1> {};
struct MU_EFF : public descriptors::FIELD_BASE<1> {};
struct ANISOTROPY_FACTOR : public descriptors::FIELD_BASE<1> {};
struct MCVALUE : public descriptors::FIELD_BASE<1> {};
struct TOTAL_ENERGY : public descriptors::FIELD_BASE<1> {};
struct INTENSITY : public descriptors::FIELD_BASE<1> {};

struct PART_RADIUS : public descriptors::FIELD_BASE<1> {};
struct PART_RHO : public descriptors::FIELD_BASE<1> {};

// Concentration/component vector parameters
struct PHYS_DENSITY_COMPONENTS : public descriptors::FIELD_BASE<3> { };
struct DIFFUSION_COMPONENTS : public descriptors::FIELD_BASE<3> { };

// Electricity
struct DIFFUSION : public descriptors::FIELD_BASE<1> { };
struct VALENCE : public descriptors::FIELD_BASE<1> { };
struct DIELECTRIC_CONSTANT : public descriptors::FIELD_BASE<1> { };
struct DEBYE : public descriptors::FIELD_BASE<1> { };

// Species-specific time parameters (laminarReactiveTmixer)
struct PHYS_MAX_T_AD : public descriptors::FIELD_BASE<1> { };

// Optimization
struct INITIAL_CONTROL_SCALAR : public descriptors::FIELD_BASE<1> { };
struct REGULARIZATION_FACTOR : public descriptors::FIELD_BASE<1> { };
struct TARGET_VOLUME : public descriptors::FIELD_BASE<1> { };
struct FDQ_STENCIL_WIDTH : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<FDQ_STENCIL_WIDTH>()>(T(1e-8));
  }
};
struct CDQ_STENCIL_WIDTH : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<FDQ_STENCIL_WIDTH>()>(T(5e-6));
  }
};

// Primitives
struct FACTOR : public descriptors::FIELD_BASE<1> { };
struct OMEGA : public descriptors::FIELD_BASE<1> { };
struct DT : public descriptors::FIELD_BASE<1> { };
struct DX : public descriptors::FIELD_BASE<1> { };
struct LATTICE_VISCOSITY : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_PRESSURE : public descriptors::FIELD_BASE<1> { };
struct REG_ALPHA : public descriptors::FIELD_BASE<1> { };
struct NORMALIZE : public descriptors::FIELD_BASE<1> { };
struct BOUNDARY_FACES : public descriptors::FIELD_BASE<0,  1, 0> { };
struct CONVERSION_FORCE : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_DRAG : public descriptors::FIELD_BASE<1> { };
struct LATTICE_TIME : public descriptors::TYPED_FIELD_BASE<std::size_t,1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct LATTICE_REFERENCE_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct LATTICE_REFERENCE_PRESSURE : public descriptors::FIELD_BASE<1> { };
struct SMAGORINSKY_CONSTANT          : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_DISPLACEMENT : public descriptors::FIELD_BASE<1> { };

}

}

#endif
