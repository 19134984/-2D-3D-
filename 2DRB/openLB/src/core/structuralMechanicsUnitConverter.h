/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Florian Kaiser
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

#ifndef STRUCTURAL_MECHANICS_UNIT_CONVERTER_H
#define STRUCTURAL_MECHANICS_UNIT_CONVERTER_H

#include <math.h>
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"
#include "core/unitConverter.h"

namespace olb {

/** Conversion between physical and lattice units, as well as discretization specialized for applications in structural mechanics
  * Be aware of the nomenclature:
  * We distingish between physical (dimensioned) and lattice (dimensionless) values.
  * A specific conversion factor maps the two different scopes,
 * e.g. __physLength = conversionLength * latticeLength__
**/

/**
* @brief Base unit converter providing core physical and material parameters for structural mechanics solvers.
* @tparam T Floating-point data type.
* @tparam DESCRIPTOR Lattice descriptor defining velocity set and spatial dimension.
*/
  SolidCase solidCaseFromString(const std::string& name)
  {
    if (name == "planeStrain") { return SolidCase::planeStrain; }
    if (name == "planeStress") { return SolidCase::planeStress; }
    if (name == "continua")    { return SolidCase::continua;    }
    if (name == "standard")    { return SolidCase::standard;    }
    throw std::runtime_error("Unknown solid case: " + name);
  }

  template <typename T, typename DESCRIPTOR>
  class StructuralMechanicsUnitConverter : public UnitConverter<T, DESCRIPTOR> {
    public:
      constexpr StructuralMechanicsUnitConverter(
        T physDeltaX,
        T physDeltaT,
        T charPhysLength,
        T charPhysTime,
        T charPhysFirstMoment,
        T youngsModulus,
        T poissonRatio,
        T physDensity,
        const std::string& solidCase
      ) : UnitConverter<T, DESCRIPTOR>(
        physDeltaX,
        physDeltaT,
        charPhysLength,
        charPhysFirstMoment,
        poissonRatio,
        physDensity,
        youngsModulus
      )
      {
        this->_charPhysTime  = charPhysTime;
        this->_poissonRatio  = poissonRatio;
        this->_youngsModulus = youngsModulus;
        this->_solidCase     = solidCaseFromString(solidCase);

        switch (this->_solidCase) {
          case SolidCase::planeStress:
            this->_lameLambda = youngsModulus * poissonRatio / (T{1.0} - poissonRatio * poissonRatio);
            break;
          case SolidCase::planeStrain:
          case SolidCase::standard:
            this->_lameLambda = youngsModulus * poissonRatio
                              / ((T{1.0} + poissonRatio) * (T{1.0} - T{2.0} * poissonRatio));
            break;
          case SolidCase::continua:
            this->_lameLambda = youngsModulus * poissonRatio
                              / ((T{1.0} + poissonRatio) * (T{1.0} - poissonRatio));
            break;
        }
        this->_lameMu      = youngsModulus / (T{2.0} * (T{1.0} + poissonRatio));
        this->_bulkModulus = this->_lameLambda + T{2.0} / T(DESCRIPTOR::d) * this->_lameMu;
      };

      const std::string getSolidCase() const {
        switch (this->_solidCase) {
          case SolidCase::planeStrain:
            return "planeStrain";
          case SolidCase::planeStress:
            return "planeStress";
          case SolidCase::continua:
            return "continua";
          case SolidCase::standard:
            return "standard";
        }
        return std::string{};
      }

      constexpr T getCharPhysLength(  ) const {
        return this->_charPhysLength;
      };

      constexpr T getCharPhysTime(  ) const {
        return this->_charPhysTime;
      };

      constexpr T getPoissonRatio(  ) const {
        return this->_poissonRatio;
      };

      constexpr T getYoungsModulus(  ) const {
        return this->_youngsModulus;
      };

      constexpr T getPhysDensity(  ) const {
        return this->_physDensity;
      };

      constexpr T getBulkModulus(  ) const {
        return UnitConverterBase<T>::getBulkModulus();
      }

      constexpr T getLameMu(  ) const {
        return UnitConverterBase<T>::getLameMu();
      }

      constexpr T getLameLambda(  ) const {
        return UnitConverterBase<T>::getLameLambda();
      }

      void print() const override;
      void print(std::ostream& fout) const override;
      void write(const std::string& fileName) const override;

  };

  /**
 * @brief Unit converter for the vector-based elastodynamics LBM scheme (Boolakee et al., 2025).
 *
 * Inherits from StructuralMechanicsUnitConverter to compute and manage conversion factors
 * between physical and lattice units specifically for the vector formulation.
 */
  template <typename T, typename DESCRIPTOR>
  class ElastodynamicsVectorLBMUnitConverter : public StructuralMechanicsUnitConverter<T, DESCRIPTOR> {
    public:

      constexpr ElastodynamicsVectorLBMUnitConverter(
        T physDeltaX,
        T physDeltaT,
        T charPhysLength,
        T charPhysTime,
        T charPhysVelocity,
        T youngsModulus,
        T poissonRatio,
        T physDensity,
        const std::string& solidCase
      ) : StructuralMechanicsUnitConverter<T, DESCRIPTOR>(
          physDeltaX,
          physDeltaT,
          charPhysLength,
          charPhysTime,
          charPhysVelocity,
          youngsModulus,
          poissonRatio,
          physDensity,
          solidCase
        )
        {
          this->_latticeSpeed           =  physDeltaX / physDeltaT                                          ;
          this->_charPhysDisplacement   =  charPhysVelocity * charPhysTime                                  ;
          this->_conversionDisplacement =  physDeltaX                                                       ;
          this->_conversionSolidStress  = -physDeltaX / physDeltaT * physDensity                            ;
          this->_conversionSolidForce   =  physDeltaX / physDeltaT * physDensity / physDeltaT               ;
          this->_conversionTraction     =  physDeltaX * physDeltaX * physDensity / (physDeltaT * physDeltaT);
          this->_conversionLame         =  physDeltaX * physDeltaX * physDensity / (physDeltaT * physDeltaT);
          this->_cMu                    =  util::sqrt(T(this->_lameMu)      / physDensity);
          this->_cK                     =  util::sqrt(T(this->_bulkModulus) / physDensity);
        };

      void print() const override;
      void print(std::ostream& fout) const override;
      void write(const std::string& fileName) const override;
  };

  /**
 * @brief Unit converter for the quasistationary MRT LBM scheme (Boolakee et al., 2023).
 *
 * Inherits from StructuralMechanicsUnitConverter to compute and manage conversion factors
 * between physical and lattice units specifically for the quasistationary MRT formulation.
 */
  template <typename T, typename DESCRIPTOR>
  class StationaryMRTUnitConverter : public StructuralMechanicsUnitConverter<T, DESCRIPTOR> {
    public:

      constexpr StationaryMRTUnitConverter(
        T physDeltaX,
        T physDeltaT,
        T charPhysLength,
        T charPhysTime,
        T charPhysDisplacement,
        T youngsModulus,
        T poissonRatio,
        T physDensity,
        T dampingFactor,
        const std::string& solidCase)
        : StructuralMechanicsUnitConverter<T, DESCRIPTOR>(
            physDeltaX,
            physDeltaT,
            charPhysLength,
            charPhysTime,
            charPhysDisplacement,
            youngsModulus,
            poissonRatio,
            physDensity,
            solidCase
      )
      {
        this->_charPhysDisplacement   = charPhysDisplacement;
        this->_conversionDisplacement = physDeltaX;

        this->_conversionSolidStress  = dampingFactor * physDeltaX * physDeltaX / physDeltaT;
        this->_conversionSolidForce   = dampingFactor * physDeltaX / physDeltaT;
        this->_conversionTraction     = dampingFactor * physDeltaX * physDeltaX / physDeltaT;

        this->_conversionLame         = dampingFactor * physDeltaX * physDeltaX / physDeltaT;
        this->_dampingFactor          = dampingFactor;
        this->_epsilon = physDeltaX / charPhysLength;
      };

      T computeRelaxationTime(const std::string& id) const {
        if (id == "11") {
          return getLatticeLameMu() / (T(1) / descriptors::invCs2<T, DESCRIPTOR>());
        }
        else if (id == "s") {
          return getLatticeBulkModulus() * T(2) / (T(1) + T(1) / descriptors::invCs2<T, DESCRIPTOR>());
        }
        else if (id == "d") {
          return getLatticeBulkModulus() * T(2) / (T(1) - T(1) / descriptors::invCs2<T, DESCRIPTOR>());
        }
        else if (id == "12" || id == "21" || id == "22" || id == "f") {
          return T(0.5);
        }
        else {
          return T(0);
        }
      }

      constexpr T getGamma(  ) const {
        T tau_s = computeRelaxationTime("s");
        T tau_f = computeRelaxationTime("f");
        T theta = T(1) / descriptors::invCs2<T, DESCRIPTOR>();
        return theta * tau_f / ((T(1) + theta) * (tau_s - tau_f));
      }

      constexpr T getCharPhysDisplacement(  ) const {
        return this->_charPhysDisplacement;
      };

      template < typename DISPLACEMENT >
      constexpr DISPLACEMENT getLatticeDisplacement( DISPLACEMENT physDisplacement ) const {
        return this->getLatticeDisplacement( physDisplacement );
      }

      template < typename DISPLACEMENT >
      constexpr DISPLACEMENT getPhysDisplacement( DISPLACEMENT latticeDisplacement ) const {
        return this->getPhysDisplacement( latticeDisplacement );
      }

      constexpr T getConversionFactorDisplacement( ) const {
        return this->_conversionDisplacement;
      }

      constexpr T getLatticeSolidStress( T physStress ) const {
        return this->getLatticeSolidStress( physStress );
      }

      constexpr T getPhysSolidStress( T latticeStress ) const {
        return this->getPhysSolidStress( latticeStress );
      }

      constexpr T getConversionFactorSolidStress( ) const {
        return this->_conversionSolidStress;
      }

      template < typename FORCE >
      constexpr FORCE getLatticeForce( FORCE physForce ) const {
        return this->getSolidLatticeForce( physForce );
      }

      template < typename FORCE >
      constexpr FORCE getPhysForce( FORCE latticeForce ) const {
        return this->getPhysSolidForce( latticeForce );
      }

      constexpr T getConversionFactorSolidForce() const {
        return this->_conversionSolidForce;
      }

      template < typename TRACTION >
      constexpr TRACTION getLatticeTraction( TRACTION physTraction ) const {
        return this->getLatticeTraction( physTraction );
      }

      template < typename TRACTION >
      constexpr TRACTION getPhysTraction( TRACTION latticeTraction ) const {
        return this->getPhysTraction( latticeTraction );
      }

      constexpr T getConversionFactorTraction() const {
        return this->_conversionTraction;
      }

      constexpr T getLatticeLameMu( ) const {
        return UnitConverterBase<T>::getLatticeLameMu();
      }

      constexpr T getLatticeLameLambda(  ) const {
        return UnitConverterBase<T>::getLatticeLameLambda();
      }

      constexpr T getLatticeBulkModulus(  ) const {
        return UnitConverterBase<T>::getLatticeBulkModulus();
      }

      constexpr T getConversionFactorLame(  ) const {
        return this->_conversionLame;
      }

      constexpr T getLatticeYoungsModulus(  ) const {
        return UnitConverterBase<T>::getLatticeYoungsModulus();
      }

      constexpr T getDampingFactor(  ) const {
        return this->_dampingFactor;
      }

      const std::string isStable() const {
        return UnitConverterBase<T>::quasistationaryIsStable();
      }

      constexpr T getMaxDt() const {
        return UnitConverterBase<T>::quasistationaryGetMaxDt();
      }

      void print() const override;
      void print(std::ostream& fout) const override;
      void write(const std::string& fileName) const override;
  };



}; // namespace olb

#endif
