/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Shota Ito, Johannes Grafen
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

#ifndef PHYS_F_H
#define PHYS_F_H

namespace olb {

namespace functors {

/// Computes the physical velocity
struct VelocityF {

  using parameters_t = meta::list<parameters::CONVERSION_VELOCITY>;
  using result_t = descriptors::VELOCITY;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    V conversion = parameters.template get<parameters::CONVERSION_VELOCITY>();
    V u[DESCRIPTOR::d]{V{0}};
    cell.computeU(u);

    FieldD<V,DESCRIPTOR,VelocityF::result_t> physU{};
    for (auto iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
      physU[iDim] = u[iDim] * conversion;
    }
    return physU;
  }
};

struct ConcentrationF {

  using parameters_t = meta::list<>;
  using result_t = descriptors::RHO;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    FieldD<V,DESCRIPTOR,ConcentrationF::result_t> concentration{cell.computeRho()};
    return concentration;
  }
};

struct DissipationF {

  using parameters_t = meta::list<parameters::OMEGA,
                                  parameters::DT,
                                  parameters::PHYS_CHAR_VISCOSITY>;

  using result_t = descriptors::DISSIPATION;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V omega = parameters.template get<parameters::OMEGA>();
    const V dt = parameters.template get<parameters::DT>();
    const V physViscosity = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();

    V rho = V{0}; V u[DESCRIPTOR::d]{V{0}}; V pi[util::TensorVal<DESCRIPTOR>::n]{V{0}};
    cell.computeAllMomenta(rho, u, pi);

    V PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
    if (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4] + pi[5] * pi[5];
    }

    FieldD<V,DESCRIPTOR,descriptors::DISSIPATION> dissipation{};
    dissipation[0] = PiNeqNormSqr *
                     util::pow(omega * descriptors::invCs2<V,DESCRIPTOR>() / rho / dt, 2) /
                     V(2) * physViscosity;
    return dissipation;
  }
};

struct PorousDissipationF {
  using parameters_t = meta::list<parameters::OMEGA,
                                  parameters::DX,
                                  parameters::LATTICE_VISCOSITY,
                                  parameters::CONVERSION_VELOCITY,
                                  parameters::PHYS_CHAR_VISCOSITY>;
  using result_t = descriptors::POROUS_DISSIPATION;
  using fields_t = meta::list<descriptors::POROSITY>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const V omega = parameters.template get<parameters::OMEGA>();
    const V dx = parameters.template get<parameters::DX>();
    const V viscosity = parameters.template get<parameters::LATTICE_VISCOSITY>();
    const V physViscosity = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
    const V conversionVelocity = parameters.template get<parameters::CONVERSION_VELOCITY>();

    V u[DESCRIPTOR::d]{V{0}};
    cell.computeU(u);
    const V gridTerm = dx * dx * viscosity / omega;
    const V invPermeability = (V(1) - porosity) / gridTerm;
    const V uNormSq = util::euklidN2(u, DESCRIPTOR::d) * conversionVelocity * conversionVelocity;

    FieldD<V,DESCRIPTOR,descriptors::POROUS_DISSIPATION> dissipation{};
    dissipation[0] = physViscosity * invPermeability * uNormSq;
    return dissipation;
  }
};

struct BoundaryForceF {
  using parameters_t = meta::list<parameters::CONVERSION_FORCE>;

  using result_t = descriptors::BOUNDARY_FORCE;
  using fields_t = meta::list<descriptors::MATERIAL>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V conversionForce = parameters.template get<parameters::CONVERSION_FORCE>();

    FieldD<V,DESCRIPTOR,BoundaryForceF::result_t> force{};
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      const Vector<int,DESCRIPTOR::d> c = descriptors::c<DESCRIPTOR>(iPop);
      if (cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::MATERIAL>(0) == 1 &&
          cell.template getFieldComponent<descriptors::MATERIAL>(0) != 1) {
        V f = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop))[iPop];
        f += cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        // Update force
        for (int i = 0; i < DESCRIPTOR::d; ++i) {
          force[i] -= c[i] * f;
        }
      }
    }

    force *= conversionForce;
    return force;
  }
};


// Drag force comuptation based on the momentum exchange algorithm.
struct DragForceF {
  using parameters_t = typename BoundaryForceF::parameters_t::template include<parameters::CONVERSION_DRAG,
                                                                             parameters::DX,
                                                                             parameters::BOUNDARY_FACES>;

  using result_t = descriptors::DRAG_FORCE;
  using fields_t = typename BoundaryForceF::fields_t;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V conversionDrag = parameters.template get<parameters::CONVERSION_DRAG>();
    const V dx = parameters.template get<parameters::DX>();
    const auto faces = parameters.template get<parameters::BOUNDARY_FACES>();

    FieldD<V,DESCRIPTOR,BoundaryForceF::result_t> boundaryForce = BoundaryForceF{}.compute(cell, parameters);
    FieldD<V,DESCRIPTOR,DragForceF::result_t> force{};
    force = conversionDrag * boundaryForce / (util::pow(dx, DESCRIPTOR::d));
    for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] /= faces[iD];
    }
    return force;
  }
};

// Far field drag approximation introduced by Onorato et al. 1986, DOI: 10.1016/0167-6105(86)90094-2
// This functor computes the far field drag force component in FLOW_DIRECTION.
// The reference state is the state of an undisturbed free flow, see paper for more detailed explanation.
// Cartesian direction is passed via DIRECTION with ()
template <unsigned int DIRECTION = 0>
struct FarFieldDragF {
  using parameters_t = meta::list<parameters::CONVERSION_FORCE,
                                  parameters::LATTICE_REFERENCE_VELOCITY,
                                  parameters::LATTICE_REFERENCE_PRESSURE>;

  using result_t = descriptors::DRAG;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    static_assert(DIRECTION <= DESCRIPTOR::d, "functors::FarFieldDragF: Invalid direction.");

    const V conversionForce = parameters.template get<parameters::CONVERSION_FORCE>();
    const V u_0 = parameters.template get<parameters::LATTICE_REFERENCE_VELOCITY>();
    const V p_0 = parameters.template get<parameters::LATTICE_REFERENCE_PRESSURE>();

    V rho{V{0}};
    V u[DESCRIPTOR::d]{V{0}};
    V u_x{V{0}};
    cell.computeRhoU(rho, u);
    rho -= V{1.0}; // Undo the population shift
    V p = rho / descriptors::invCs2<V,DESCRIPTOR>(); // Compute pressure

    // Extract the velocity component in flow direction
    if constexpr (DIRECTION == 0) {
      u_x = u[0];
    } else if (DIRECTION == 1) {
      u_x = u[1];
    } else {
      u_x = u[2];
    }

    FieldD<V,DESCRIPTOR,FarFieldDragF::result_t> drag{};
    drag = (rho * (u_x*u_0 - u_x*u_x) + (p_0 - p)) * conversionForce;
    return drag;
  }
};

template <unsigned int DIRECTION = 0>
struct FarFieldDragForceF {
  using parameters_t = typename FarFieldDragF<DIRECTION>::parameters_t::template include<parameters::CONVERSION_DRAG,
                                                                                         parameters::DX,
                                                                                         parameters::BOUNDARY_FACES>;

  using result_t = descriptors::DRAG_FORCE_COMPONENT;
  using fields_t = typename FarFieldDragF<DIRECTION>::fields_t;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V conversionDrag = parameters.template get<parameters::CONVERSION_DRAG>();
    const V dx = parameters.template get<parameters::DX>();
    const auto faces = parameters.template get<parameters::BOUNDARY_FACES>();

    FieldD<V,DESCRIPTOR,typename FarFieldDragF<DIRECTION>::result_t> boundaryForce = FarFieldDragF<DIRECTION>{}.compute(cell, parameters);
    FieldD<V,DESCRIPTOR,typename FarFieldDragForceF<DIRECTION>::result_t> force{};
    force = conversionDrag * boundaryForce / (util::pow(dx, DESCRIPTOR::d));
    force /= faces[DIRECTION];
    return force;
  }
};

template <typename DYNAMICS>
struct SmagorinskyTurbulentTauF {

  using parameters_t = meta::list<parameters::OMEGA, parameters::SMAGORINSKY_CONSTANT>;
  using result_t = descriptors::AVERAGE_TURBULENT_TAU;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    using MomentaF = typename DYNAMICS::MomentaF;
    V piNeqNormSqr { };
    MomentaF().computePiNeqNormSqr(cell, piNeqNormSqr);
    const V rho = MomentaF().computeRho(cell);
    const V omega = parameters.template get<parameters::OMEGA>();
    const V smagorinsky = parameters.template get<parameters::SMAGORINSKY_CONSTANT>();
    V piNeqNorm = util::sqrt(piNeqNormSqr);
    V preFactor = smagorinsky*smagorinsky
                * descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()
                * 2 * util::sqrt(2);
    /// Molecular realaxation time
    V tauMol = V{1} / omega;
    /// Turbulent realaxation time
    V tauTurb = V{0.5} * (util::sqrt(tauMol*tauMol + preFactor / rho * piNeqNorm) - tauMol);

    return FieldD<V,DESCRIPTOR,SmagorinskyTurbulentTauF::result_t>{tauTurb};
  }
};

template <unsigned MATERIAL1, unsigned MATERIAL2>
struct PressureDropF {
  using parameters_t = meta::list<parameters::CONVERSION_PRESSURE>;
  using result_t = descriptors::PRESSURE_LEVEL;
  using fields_t = meta::list<descriptors::MATERIAL>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    V rho = cell.computeRho();
    V conversion = parameters.template get<parameters::CONVERSION_PRESSURE>();
    V pressure = V(1) / descriptors::invCs2<V,DESCRIPTOR>() * (rho - 1) * conversion;
    int material = cell.template getField<descriptors::MATERIAL>();
    FieldD<V,DESCRIPTOR,descriptors::PRESSURE_LEVEL> pressuredrop{};
    pressuredrop[0] = pressure * (V(material == MATERIAL1) - V(material == MATERIAL2));
    return pressuredrop;
  }
};

/// Computes the physical displacement of the quasistationary linear elasticity MRT scheme
/**
 * Mirrors the moment space collision of structuralmechanics::StationaryLinearElasticityMRT
 * (see src/structuralmechanics/quasistationary/boolakee_quasistationary.h): the displacement
 * is the first order bared moment, recomputed here from the cell's populations so that it
 * stays an explicit, differentiable function of the state.
 */
struct DisplacementF {
  using parameters_t = meta::list<parameters::CONVERSION_DISPLACEMENT, parameters::OMEGA_SOLID>;
  using result_t = descriptors::DISPLACEMENT;
  using fields_t = meta::list<descriptors::FORCE, opti::CONTROL_VECTOR>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V          = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    const V                        conversion = parameters.template get<parameters::CONVERSION_DISPLACEMENT>();
    const Vector<V, DESCRIPTOR::q> Omega      = parameters.template get<parameters::OMEGA_SOLID>();
    const V                        theta      = V(1) / descriptors::invCs2<V, DESCRIPTOR>();

    const Vector<V, DESCRIPTOR::d> f1 = cell.template getField<descriptors::FORCE>();
    // Optional component-wise scaling of the force, one unless set (see opti::CONTROL_VECTOR)
    const Vector<V, DESCRIPTOR::d> control = cell.template getField<opti::CONTROL_VECTOR>();
    Vector<V, DESCRIPTOR::q> force = {
      control[0] * f1[0],
      control[1] * f1[1],
      V(0.), V(0.), V(0.), V(0.), V(0.), V(0.)
    };

    Vector<Vector<V, DESCRIPTOR::q>, DESCRIPTOR::q> M;
    M[0] = Vector<V, DESCRIPTOR::q>(1,  0, -1,  0, 1, -1, -1,  1);
    M[1] = Vector<V, DESCRIPTOR::q>(0,  1,  0, -1, 1,  1, -1, -1);
    M[2] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1, -1,  1, -1);
    M[3] = Vector<V, DESCRIPTOR::q>(1,  1,  1,  1, 2,  2,  2,  2);
    M[4] = Vector<V, DESCRIPTOR::q>(1, -1,  1, -1, 0,  0,  0,  0);
    M[5] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1, -1, -1,  1);
    M[6] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1,  1, -1, -1);
    M[7] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1,  1,  1,  1);

    Vector<V, DESCRIPTOR::q> cell_vec([&cell](unsigned i) { return cell[i]; });
    Vector<V, DESCRIPTOR::q> m = M * cell_vec + V(0.5) * force;

    Vector<Vector<V, DESCRIPTOR::d>, DESCRIPTOR::q> momentaFactors;
    momentaFactors[0] = Vector<V, DESCRIPTOR::d>(    1,     0);
    momentaFactors[1] = Vector<V, DESCRIPTOR::d>(    0,     1);
    momentaFactors[2] = Vector<V, DESCRIPTOR::d>(    0,     0);
    momentaFactors[3] = Vector<V, DESCRIPTOR::d>(    0,     0);
    momentaFactors[4] = Vector<V, DESCRIPTOR::d>(    0,     0);
    momentaFactors[5] = Vector<V, DESCRIPTOR::d>(theta,     0);
    momentaFactors[6] = Vector<V, DESCRIPTOR::d>(    0, theta);
    momentaFactors[7] = Vector<V, DESCRIPTOR::d>(    0,     0);

    Vector<V, 2> m_first_two = {
      m[0],
      m[1]
    };

    Vector<V, DESCRIPTOR::q> meq = momentaFactors * m_first_two;

    Vector<V, DESCRIPTOR::q> bracket{};
    for (unsigned iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      bracket[iPop] = V(0.5) * force[iPop] + m[iPop] + Omega[iPop] * (meq[iPop] - m[iPop]);
    }

    Vector<V, DESCRIPTOR::q> m_bared = V(0.5) * (m + bracket);

    FieldD<V, DESCRIPTOR, DisplacementF::result_t> physDisp {};
    for (unsigned iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
      physDisp[iDim] = m_bared[iDim] * conversion;
    }
    return physDisp;
  }
};

} // namespace functors
}

#endif
