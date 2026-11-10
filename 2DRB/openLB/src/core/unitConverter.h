/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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

#ifndef CORE_UNITCONVERTER_H
#define CORE_UNITCONVERTER_H

#include "utilities/omath.h"
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"
#include "descriptor/fields.h"
#include "utilities/optionalValue.h"
#include "core/operator.h"
#include "core/concepts.h"
#include "case/parametersD.h"

// known design issues
//    1. How can we prevent abuse of constructur by mixing up parameters?
//    2. physical problems may have different names for viscosity, e.g. diffusity,  temperature conductivity
//    4. Feedback about stability or comment the chosen discretization
//    5. Explain why Desctiptor as template
//    6. Is it worth to introduce invConversionDensity to avoid division


namespace olb {

enum class SolidCase {
  planeStress,
  planeStrain,
  continua,
  standard
};

namespace fields {

namespace converter {

struct PHYS_VELOCITY : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct PHYS_FORCE : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct PHYS_LENGTH : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
 };
struct PHYS_DELTA_X : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct PHYS_DELTA_T : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct PHYS_PRESSURE : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct LATTICE_TIME : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};

struct LATTICE_VISCOSITY : public descriptors::FIELD_BASE<1> { };

}

}

template <typename T>
class UnitConverterBase {
public:
  virtual ~UnitConverterBase() = default;

  virtual void print() const = 0;
  virtual void print(std::ostream& fout) const = 0;
  virtual void write(std::string const& fileName = "unitConverter") const = 0;

  /// return resolution
  int getResolution(  ) const
  {
    return _resolution;
  }
  /// return relaxation time in lattice units
  T getLatticeRelaxationTime(  ) const
  {
    return _latticeRelaxationTime;
  }
  /// return relaxation frequency in lattice units
  T getLatticeRelaxationFrequency(  ) const
  {
    return 1./_latticeRelaxationTime;
  }
  /// return relaxation frequency in lattice units computed from given physical diffusivity in __m^2 / s__
  template <typename DESCRIPTOR_>
  T getLatticeRelaxationFrequencyFromDiffusivity( const T physDiffusivity ) const
  {
    return 1.0 / ( physDiffusivity / _conversionViscosity * descriptors::invCs2<T,DESCRIPTOR_>() + 0.5 );
  }
  /// return characteristic length in physical units
  T getCharPhysLength(  ) const
  {
    return _charPhysLength;
  }
  /// return characteristic velocity in physical units
  T getCharPhysVelocity(  ) const
  {
    return _charPhysVelocity;
  }
  /// return characteristic velocity in lattice units
  T getCharLatticeVelocity(  ) const
  {
    return _charLatticeVelocity;
  }
  /// return characteristic CFL number
  T getCharCFLnumber(  ) const
  {
    return _charLatticeVelocity;
  }
  /// return viscosity in physical units
  T getPhysViscosity(  ) const
  {
    return _physViscosity;
  }
  /// return density in physical units
  T getPhysDensity(  ) const
  {
    return _physDensity;
  }
  /// return characteristic pressure in physical units
  T getCharPhysPressure(  ) const
  {
    return _charPhysPressure;
  }
  /// return Reynolds number
  T getReynoldsNumber(  ) const
  {
    return _charPhysVelocity * _charPhysLength / _physViscosity;

    // Power law
    // Calculation according to Metzner and Reed (1955): 10.1002/aic.690010409
    //return util::pow(this->_charPhysVelocity, T{2} - _powerLawIndex) * util::pow(this->_charPhysLength, _powerLawIndex)
    // / _physConsistencyCoeff;
  }
  /// return Mach number
  T getMachNumber() const
  {
    return getCharLatticeVelocity() * util::sqrt(_invCs2.value());
  }
  /// return Knudsen number
  virtual T getKnudsenNumber(  ) const
  {
    // This calculates the lattice Knudsen number.
    // See e.g. (7.22) in "The Lattice Boltzmann Method: Principles and Practice" [kruger2017lattice].
    return getMachNumber() / getReynoldsNumber();

    // ADE unit converter
    //return this->getMachNumber()/getPecletNumber();
  }
  /// conversion from lattice to  physical length
  T getPhysLength( int latticeLength ) const
  {
    return _conversionLength * latticeLength;
  }
  /// conversion from physical to lattice length, returns number of voxels for given physical length
  int getLatticeLength( T physLength ) const
  {
    return int( physLength / _conversionLength + 0.5 );
  }
  /// access (read-only) to private member variable
  T getConversionFactorLength() const
  {
    return _conversionLength;
  }
  /// returns grid spacing (voxel length) in __m__
  T getPhysDeltaX() const
  {
    return _conversionLength;
  }

  /// conversion from lattice to  physical time
  T getPhysTime( size_t latticeTime ) const
  {
    return _conversionTime * latticeTime;
  }
  /// conversion from physical to lattice time
  size_t getLatticeTime( T physTime ) const
  {
    return size_t(physTime / _conversionTime + 0.5);
  }
  /// access (read-only) to private member variable
  T getConversionFactorTime() const
  {
    return _conversionTime;
  }
  /// returns time spacing (timestep length) in __s__
  T getPhysDeltaT() const
  {
    return _conversionTime;
  }

  /// conversion from lattice to  physical velocity
  T getPhysVelocity( T latticeVelocity ) const
  {
    return _conversionVelocity * latticeVelocity;
  }
  /// conversion from physical to lattice velocity
  T getLatticeVelocity( T physVelocity ) const
  {
    return physVelocity / _conversionVelocity;
  }
  /// conversion from physical to lattice velocity
  template <unsigned D>
  Vector<T,D> getLatticeVelocity(Vector<T,D> physU) const
  {
    return Vector<T,D>([&](std::size_t iD) -> T {
      return this->getLatticeVelocity(physU[iD]);
    });
  }
  /// access (read-only) to private member variable
  T getConversionFactorVelocity() const
  {
    return _conversionVelocity;
  }

  /// conversion from lattice to  physical density
  T getPhysDensity( T latticeDensity ) const
  {
    return _conversionDensity * latticeDensity;
  }
  /// conversion from physical to lattice density
  T getLatticeDensity( T physDensity ) const
  {
    return physDensity / _conversionDensity;
  }
  T getLatticeDensityFromPhysPressure( T physPressure ) const
  {
    return getLatticePressure(physPressure) * _invCs2 + 1.0;
  }
  /// access (read-only) to private member variable
  T getConversionFactorDensity() const
  {
    return _conversionDensity;
  }

  /// conversion from lattice to  physical mass
  T getPhysMass( T latticeMass ) const
  {
    return _conversionMass * latticeMass;
  }
  /// conversion from physical to lattice mass
  T getLatticeMass( T physMass ) const
  {
    return physMass / _conversionMass;
  }
  /// access (read-only) to private member variable
  T getConversionFactorMass() const
  {
    return _conversionMass;
  }

  /// conversion from lattice to  physical viscosity
  T getPhysViscosity( T latticeViscosity ) const
  {
    return _conversionViscosity * latticeViscosity;
  }
  /// conversion from physical to lattice viscosity
  T getLatticeViscosity(  ) const
  {
    return _physViscosity / _conversionViscosity;
  }
  /// access (read-only) to private member variable
  T getConversionFactorViscosity() const
  {
    return _conversionViscosity;
  }

  /// conversion from lattice to  physical force
  T getPhysForce( T latticeForce ) const
  {
    return _conversionForce * latticeForce;
  }
  /// conversion from lattice to  physical force vector
  template <unsigned D>
  Vector<T,D> getPhysForce(Vector<T,D> latticeForce) const
  {
    return Vector<T,D>([&](std::size_t iD) -> T {
      return this->getPhysForce(latticeForce[iD]);
    });
  }
  /// conversion from physical to lattice force
  T getLatticeForce( T physForce ) const
  {
    return physForce / _conversionForce;
  }
  /// access (read-only) to private member variable
  T getConversionFactorForce() const
  {
    return _conversionForce;
  }

  /// conversion from lattice to  physical torque
  T getPhysTorque ( T latticeTorque ) const
  {
    return _conversionTorque * latticeTorque;
  }
  /// conversion from lattice to  physical force vector
  template <unsigned D>
  Vector<T,D> getPhysTorque(Vector<T,D> latticeTorque) const
  {
    return Vector<T,D>([&](std::size_t iD) -> T {
      return this->getPhysTorque(latticeTorque[iD]);
    });
  }
  /// conversion from physical to lattice torque
  T getLatticeTorque( T physTorque ) const
  {
    return physTorque / _conversionTorque;
  }
  /// access (read-only) to private member variable
  T getConversionFactorTorque() const
  {
    return _conversionTorque;
  }

  /// conversion from lattice to  physical pressure
  T getPhysPressure( T latticePressure ) const
  {
    return _conversionPressure * latticePressure + _charPhysPressure;
  }
  /// conversion from physical to lattice pressure
  T getLatticePressure( T physPressure ) const
  {
    return ( physPressure - _charPhysPressure ) / _conversionPressure;
  }
  /// access (read-only) to private member variable
  T getConversionFactorPressure() const
  {
    return _conversionPressure;
  }

  // from thermalUnitConverter
  /// return thermal relaxation time in lattice units
  T getLatticeThermalRelaxationTime(  ) const
  {
    return _latticeThermalRelaxationTime;
  };
  /// return thermal relaxation frequency in lattice units
  T getLatticeThermalRelaxationFrequency(  ) const
  {
    return 1.0 / _latticeThermalRelaxationTime;
  };

  /// return characteristic low temperature in physical units
  T getCharPhysLowTemperature(  ) const
  {
    return _charPhysLowTemperature;
  };
  /// return characteristic high temperature in physical units
  T getCharPhysHighTemperature(  ) const
  {
    return _charPhysHighTemperature;
  };
  /// return characteristic temperature difference in physical units
  T getCharPhysTemperatureDifference(  ) const
  {
    return _charPhysTemperatureDifference;
  };
  /// return thermal expansion coefficient in physical units
  T getPhysThermalExpansionCoefficient(  ) const
  {
    return _physThermalExpansionCoefficient;
  };
  /// return thermal diffusivity in physical units
  T getPhysThermalDiffusivity(  ) const
  {
    return _physThermalDiffusivity;
  };
  /// return specific heat capacity in physical units
  T getPhysSpecificHeatCapacity(  ) const
  {
    return _physSpecificHeatCapacity;
  };
  /// return thermal conductivity in physical units
  T getThermalConductivity(  ) const
  {
    return _physThermalConductivity;
  };

  /// conversion from lattice to physical temperature
  T getPhysTemperature( T latticeTemperature ) const
  {
    return _conversionTemperature * (latticeTemperature - 0.5) + _charPhysLowTemperature;
  };
  /// conversion from physical to lattice temperature
  T getLatticeTemperature( T physTemperature ) const
  {
    return (physTemperature - _charPhysLowTemperature) / _conversionTemperature + 0.5;
  };

  /// conversion from lattice to physical thermal diffusivity
  T getPhysThermalDiffusivity( T latticeThermalDiffusivity ) const
  {
    return _conversionThermalDiffusivity * latticeThermalDiffusivity;
  };
  /// conversion from physical to lattice thermal diffusivity
  T getLatticeThermalDiffusivity( T physThermalDiffusivity ) const
  {
    return physThermalDiffusivity / _conversionThermalDiffusivity;
  };
  /// access (read-only) to private member variable
  T getConversionFactorThermalDiffusivity() const
  {
    return _conversionThermalDiffusivity;
  };


  /// conversion from lattice to physical specific heat capacity
  T getPhysSpecificHeatCapacity( T latticeSpecificHeatCapacity ) const
  {
    return _conversionSpecificHeatCapacity * latticeSpecificHeatCapacity;
  };
  /// conversion from physical to lattice specific heat capacity
  T getLatticeSpecificHeatCapacity( T physSpecificHeatCapacity ) const
  {
    return physSpecificHeatCapacity / _conversionSpecificHeatCapacity;
  };
  /// access (read-only) to private member variable
  T getConversionFactorSpecificHeatCapacity() const
  {
    return _conversionSpecificHeatCapacity;
  };

  /// conversion from lattice to physical thermal  conductivity
  T getPhysThermalConductivity( T latticeThermalConductivity ) const
  {
    return _conversionThermalConductivity * latticeThermalConductivity;
  };
  /// conversion from physical to lattice thermal  conductivity
  T getLatticeThermalConductivity( T physThermalConductivity ) const
  {
    return physThermalConductivity / _conversionThermalConductivity;
  };
  /// access (read-only) to private member variable
  T getConversionFactorThermalConductivity() const
  {
    return _conversionThermalConductivity;
  };

  /// conversion from lattice to physical heat flux
  T getPhysHeatFlux( T latticeHeatFlux ) const
  {
    return _conversionHeatFlux * latticeHeatFlux;
  };
  /// conversion from physical to lattice heat flux
  T getLatticeHeatFlux( T physHeatFlux ) const
  {
    return physHeatFlux / _conversionHeatFlux;
  };
  /// access (read-only) to private member variable
  T getConversionFactorHeatFlux() const
  {
    return _conversionHeatFlux;
  };
  T getPrandtlNumber() const
  {
    return this->_physViscosity/_physThermalDiffusivity;
  };
  T getRayleighNumber() const
  {
    return 9.81 * _physThermalExpansionCoefficient/this->_physViscosity/_physThermalDiffusivity * (_charPhysHighTemperature - _charPhysLowTemperature) * util::pow(this->_charPhysLength.value(),3);
  };

  // from multiPhaseUnitConverter
  /// return characteristic temperature in physical units
  T getCharPhysTemperature(  ) const {
    return _charPhysTemperature;
  };
  /// return equation of state parameter a in physical units
  T getPhysEoSa(  ) const {
    return _physEoSa;
  };
  /// return equation of state parameter b in physical units
  T getPhysEoSb(  ) const {
    return _physEoSb;
  };
  /// return molar mass in physical units
  T getPhysMolarMass(  ) const {
    return _physMolarMass;
  };
  /// return surface tension in physical units
  T getPhysSurfaceTension(  ) const {
    return _physSurfaceTension;
  };
  /// return characteristic temperature in physical units
  T getPhysTemperature(  ) const {
    return _charPhysTemperature;
  };
  /// access (read-only) to private member variable
  T getConversionFactorEoSa() const {
    return _conversionEoSa;
  };
  /// access (read-only) to private member variable
  T getConversionFactorEoSb() const {
    return _conversionEoSb;
  };
  /// access (read-only) to private member variable
  T getConversionFactorMolarMass() const {
    return _conversionMolarMass;
  };
  /// access (read-only) to private member variable
  T getConversionFactorGasConstant() const {
    return _conversionGasConstant;
  };
  /// access (read-only) to private member variable
  T getConversionFactorSurfaceTension() const {
    return _conversionSurfaceTension;
  };
  /// access (read-only) to private member variable
  T getConversionFactorTemperature(  ) const {
    return _conversionTemperature;
  };
  /// return lattice surface tension for parameter fitting
  T getLatticeSurfaceTension(  ) const {
    return _latticeSurfaceTension;
  };
  /// access (read-only) to private member variable
  T getConversionFactorChemicalPotential() const {
    return _conversionChemicalPotential;
  };
  /// compute relaxation time from physical viscosity
  T computeRelaxationTimefromPhysViscosity(T userViscosity) const {
    return 0.5 + _invCs2 *
                 userViscosity / this->_conversionViscosity;
  };
  /// compute lattice surface tension from physical one
  T computeLatticeSurfaceTension(T userSurfaceTension) const {
    return userSurfaceTension / this->_conversionSurfaceTension;
  };
  /// compute Reynolds from kinematic viscosity
  T computeReynolds(T userVelocity, T userLength, T userViscosity) const {
    return userVelocity * userLength / userViscosity;
  };
  /// compute Weber
  T computeWeber(T userVelocity, T userLength, T userSurfaceTension) const {
    return userLength * userVelocity * userVelocity / userSurfaceTension;
  };

  //from powerLawUnitConverter
  /// return consistency coefficient in physical units
  T getPhysConsistencyCoeff( ) const {
    return _physConsistencyCoeff;
  }
  /// conversion from lattice to  physical consistency coefficient
  T getPhysConsistencyCoeff( T latticeConsistencyCoeff ) const {
    return _conversionConsistencyCoeff * latticeConsistencyCoeff;
  }
  /// conversion from physical to lattice consistency coefficient
  T getLatticeConsistencyCoeff(  ) const {
    return _physConsistencyCoeff / _conversionConsistencyCoeff;
  }
  /// access (read-only) to private member variable
  T getConversionFactorConsistencyCoeff() const {
    return _conversionConsistencyCoeff;
  }
  /// access (read-only) to private member variable
  T getPowerLawIndex() const {
    return _powerLawIndex;
  }

  // from adeUnitConverter
  /// return thermal relaxation time in lattice units
  T getLatticeAdeRelaxationTime() const {
    return _latticeAdeRelaxationTime;
  };
  /// return thermal relaxation frequency in lattice units
  T getLatticeAdeRelaxationFrequency() const {
    return 1.0 / _latticeAdeRelaxationTime;
  };
  T getPhysDiffusivity() const {
    return _physDiffusivity;
  }
  T getLatticeDiffusivity() const {
    return _physDiffusivity / _conversionDiffusivity ;
  }
  T getConversionFactorDiffusivity() const {
    return _conversionDiffusivity;
  }
  T getPecletNumber() const {
    return getCharPhysVelocity()  * getCharPhysLength() / getPhysDiffusivity();
  };

  // from adsorptionUnitConverter
  /// conversion factor to convert particle density from lattice units to kg/m^3
  T getConversionFactorParticleDensity() const {
    return _particleConcentrationAdsorption;
  }
  T getPhysParticleConcentration(T c) const {
    return c * _particleConcentrationAdsorption;
  }
  T getPhysConcentration(T c) const {
    return c * this->getConversionFactorDensity();
  }
  T getPhysLoading(T Cq) const {
    return Cq * _particleConcentrationAdsorption;
  }
  T getSchmidtNumber() const {
    return getPhysViscosity() / getPhysDiffusivity();
  };
  T getFourierNumber() const {
    return getPhysDiffusivity() * getPhysDeltaT() / getPhysDeltaX() / getPhysDeltaX();
  };

  // from radiativeUnitConverter
  T getPhysAbsorption() const {
    return _physAbsorption;
  };
  T getPhysScattering() const {
    return _physScattering;
  };
  T getAnisotropyFactor() const {
    return _anisotropyFactor;
  };
  T getExtinction() const {
    return _extinction;
  };
  T getScatteringAlbedo() const {
    return _scatteringAlbedo;
  };
  T getPhysDiffusion() const {
    return _physDiffusion;
  };
  T getEffectiveAttenuation() const {
    return _effectiveAttenuation;
  };
  T getLatticeAbsorption() const {
    return _latticeAbsorption;
  };
  T getLatticeScattering() const {
    return _latticeScattering;
  };
  T getLatticeDiffusion() const {
    return _latticeDiffusion;
  };
  T getRefractiveRelative() const {
    return _refractiveRelative;
  };

  // From Structural Mechanics

  // --- Base Solid Unit Converter ---- //
  T getCharPhysDisplacement() const {
    return _charPhysDisplacement;
  };

  T getConversionFactorDisplacement() const {
    return _conversionDisplacement;
  };

  T getConversionFactorTraction() const {
    return _conversionTraction;
  };

  T getCharPhysTime() const {
    return _charPhysTime;
  };

  T getPhysYoungsModulus() const {
    return _youngsModulus;
  };

  T getQuasistationarySolidDampingFactor() const {
    return _dampingFactor;
  };

  T getPoissonRatio() const {
    return _poissonRatio;
  };

  T getLatticeDisplacement( T physDisplacement ) const {
    return physDisplacement / _conversionDisplacement;
  };

  template < unsigned D >
  Vector<T, D> getLatticeDisplacement( Vector<T, D> physDisplacement ) const {
    return physDisplacement / _conversionDisplacement;
  };

  T getPhysDisplacement( T latticeDisplacement ) const {
    return latticeDisplacement * _conversionDisplacement;
  };

  template < unsigned D >
  Vector<T, D> getPhysDisplacement( Vector<T, D> latticeDisplacement ) const {
    return latticeDisplacement * _conversionDisplacement;
  };

  T getLatticeSolidStress( T physStress ) const {
    return physStress / _conversionSolidStress;
  };

  T getPhysSolidStress( T latticeStress ) const {
    return latticeStress * _conversionSolidStress;
  };

  T getConversionFactorSolidStress( ) const {
    return _conversionSolidStress;
  };

  T getLatticeSolidForce( T physForce ) const {
    return physForce / _conversionSolidForce;
  };

  template < unsigned D >
  Vector<T, D> getLatticeSolidForce( Vector<T, D> physForce ) const {
    return physForce / _conversionSolidForce;
  };

  T getPhysSolidForce( T latticeForce ) const {
    return latticeForce * _conversionSolidForce;
  };

  template < unsigned D >
  Vector<T, D> getPhysSolidForce( Vector<T, D> latticeForce ) const {
    return latticeForce * _conversionSolidForce;
  };

  T getConversionFactorSolidForce() const {
    return _conversionSolidForce;
  };

  T getLatticeTraction( T physTraction ) const {
    return physTraction / _conversionTraction;
  };

  template < unsigned D >
  Vector<T, D> getLatticeTraction( Vector<T, D> physTraction ) const {
    return physTraction / _conversionTraction;
  };

  T getPhysTraction( T latticeTraction ) const {
    return latticeTraction * _conversionTraction;
  };

  template < unsigned D >
  Vector<T, D> getPhysTraction( Vector<T, D> latticeTraction ) const {
    return latticeTraction * _conversionTraction;
  };

  T getBulkModulus(  ) const {
    return _bulkModulus;
  };

  T getLameMu(  ) const {
    return _lameMu;
  };

  T getLameLambda(  ) const {
    return _lameLambda;
  };

  T getLatticeLameMu( ) const {
    return _lameMu / _conversionLame;
  };

  T getLatticeLameLambda(  ) const {
    return _lameLambda / _conversionLame;
  };

  T getLatticeBulkModulus(  ) const {
    return _bulkModulus / _conversionLame;
  };

  T getConversionFactorLame(  ) const {
    return _conversionLame;
  };

  T getLatticeYoungsModulus(  ) const {
    return _youngsModulus / _conversionLame;
  };

  std::string getSolidCase(  ) const {
    if (! _solidCase.hasValue()) {
        throw std::runtime_error("Solid case was not set in the unit converter!");
    }
    switch (_solidCase) {
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
  };
  // --- End of: Base Solid Unit Converter ---- //


  // ---- Quasistationary Boolakee ---- //
  T getDampingFactor(  ) const {
    return _dampingFactor;
  };

  std::string quasistationaryIsStable() const {
    const T E_lattice = getLatticeYoungsModulus();
    return ((E_lattice > T(0.3) && E_lattice < T(1)) ? "true" : "false");
  }
  // ---- End of: Quasistationary Boolakee ---- //


  // ---- Elastodynamics Boolakee ---- //
  T getLatticeSpeed(  ) const {
    return _latticeSpeed;
  };

  T getCMu(  ) const {
    return _cMu;
  }

  T getCK(  ) const {
    return _cK;
  }

  T getLatticeCMu(  ) const {
    return _cMu / _conversionVelocity;
  }

  T getLatticeCK(  ) const {
    return _cK / _conversionVelocity;
  }

  template< typename DESCRIPTOR >
  std::string elastodynamicsIsStable() const {
    const T latticeCK   = getLatticeCK();
    const T latticeCMu  = getLatticeCMu();
    return (descriptors::invCs2<T, DESCRIPTOR>() * util::sqrt(latticeCK * latticeCK + latticeCMu * latticeCMu) < this->_latticeSpeed ? "true" : "false");
  }

  template< typename DESCRIPTOR >
  T elastodynamicsGetStability() const {
    const T latticeCK = getLatticeCK();
    const T latticeCMu = getLatticeCMu();
    return descriptors::invCs2<T,DESCRIPTOR>() * util::sqrt(latticeCK * latticeCK + latticeCMu * latticeCMu) / this->_latticeSpeed;
  }
  // ---- End of: Elastodynamics Boolakee ---- //


protected:

  OptionalValue<T> _invCs2;

  // conversion factors
  OptionalValue<T> _conversionLength;      // m
  OptionalValue<T> _conversionTime;        // s
  OptionalValue<T> _conversionVelocity;    // m / s
  OptionalValue<T> _conversionDensity;     // kg / m^3
  OptionalValue<T> _conversionMass;        // kg
  OptionalValue<T> _conversionViscosity;   // m^2 / s
  OptionalValue<T> _conversionForce;       // kg m / s^2
  OptionalValue<T> _conversionTorque;      // kg m^2 / s^2
  OptionalValue<T> _conversionPressure;    // kg / m s^2

  // physical units, e.g characteristic or reference values
  OptionalValue<T> _charPhysLength;        // m
  OptionalValue<T> _charPhysVelocity;      // m / s
  OptionalValue<T> _physViscosity;         // m^2 / s
  OptionalValue<T> _physDensity;           // kg / m^3
  OptionalValue<T> _charPhysPressure;      // kg / m s^2

  // lattice units, discretization parameters
  OptionalValue<std::size_t> _resolution;
  OptionalValue<T> _latticeRelaxationTime;
  OptionalValue<T> _charLatticeVelocity;   //

  // conversion factors
  OptionalValue<T> _conversionTemperature; // K
  OptionalValue<T> _conversionThermalDiffusivity; // m^2 / s
  OptionalValue<T> _conversionSpecificHeatCapacity; // J / kg K = m^2 / s^2 K
  OptionalValue<T> _conversionThermalConductivity; // W / m K = kg m / s^3 K
  OptionalValue<T> _conversionHeatFlux; // W / m^2 = kg / s^3

  // physical units, e.g characteristic or reference values
  OptionalValue<T> _charPhysLowTemperature; // K
  OptionalValue<T> _charPhysHighTemperature; // K
  OptionalValue<T> _charPhysTemperatureDifference; // K
  OptionalValue<T> _physThermalExpansionCoefficient; // 1 / K
  OptionalValue<T> _physThermalDiffusivity; // m^2 / s
  OptionalValue<T> _physSpecificHeatCapacity; // J / kg K = m^2 / s^2 K
  OptionalValue<T> _physThermalConductivity; // W / m K = kg m / s^3 K

  // lattice units, discretization parameters
  OptionalValue<T> _latticeThermalRelaxationTime; // -

  // ADE lattice units
  OptionalValue<T> _physDiffusivity;
  OptionalValue<T> _conversionDiffusivity;
  OptionalValue<T> _latticeAdeRelaxationTime;

  // MultiPhase conversion factors
  OptionalValue<T> _conversionEoSa;                      // kg m^5 / s^2 mol^2
  OptionalValue<T> _conversionEoSb;                      // m^3 / mol
  OptionalValue<T> _conversionMolarMass;                 // kg / mol
  OptionalValue<T> _conversionGasConstant = T{8.314462618}; // J / mol K = kg m^2 / s^2 mol K
  OptionalValue<T> _conversionSurfaceTension;            // J / m^2 = kg / s^2
  OptionalValue<T> _physEoSa;                            // kg m^5 / s^2 mol^2
  OptionalValue<T> _physEoSb;                            // m^3 / mol
  OptionalValue<T> _physEoSav;
  OptionalValue<T> _physEoSbv;
  OptionalValue<T> _physMolarMass;                       // kg / mol
  OptionalValue<T> _physSurfaceTension;                  // J / m^2 = kg / s^2
  OptionalValue<T> _charPhysTemperature;                 // K
  OptionalValue<T> _latticeSurfaceTension;               // -
  OptionalValue<T> _conversionChemicalPotential;         // J / kg = m^2 / s^2

  // power law conversion factors
  OptionalValue<T> _conversionConsistencyCoeff;          // m^2 s^(n-2)
  OptionalValue<T> _powerLawIndex;
  OptionalValue<T> _physConsistencyCoeff;                // m^2 s^(n-2)

  // adsorption factors
  OptionalValue<T> _physViscosityAdsorption;
  OptionalValue<T> _conversionViscosityAdsorption;
  OptionalValue<T> _particleConcentrationAdsorption;

  // radiative factors
  OptionalValue<T> _physAbsorption;
  OptionalValue<T> _physScattering;
  OptionalValue<T> _anisotropyFactor;
  OptionalValue<T> _extinction;
  OptionalValue<T> _scatteringAlbedo;
  OptionalValue<T> _physDiffusion;
  OptionalValue<T> _effectiveAttenuation;
  OptionalValue<T> _refractiveRelative;
  OptionalValue<T> _latticeAbsorption;
  OptionalValue<T> _latticeScattering;
  OptionalValue<T> _latticeDiffusion;

  // Structural mechanics factors
  OptionalValue<T> _youngsModulus;
  OptionalValue<T> _poissonRatio;
  OptionalValue<T> _lameMu;
  OptionalValue<T> _lameLambda;
  OptionalValue<T> _bulkModulus;

  OptionalValue<T> _cMu;
  OptionalValue<T> _cK;

  OptionalValue<T> _charPhysTime;
  OptionalValue<T> _charPhysDisplacement;
  OptionalValue<T> _epsilon;
  OptionalValue<T> _conversionLame;
  OptionalValue<T> _dampingFactor;
  OptionalValue<T> _conversionDisplacement;
  OptionalValue<T> _conversionTraction;
  OptionalValue<T> _conversionSolidForce;
  OptionalValue<T> _conversionSolidStress;

  OptionalValue<T> _latticeSpeed;
  OptionalValue<SolidCase> _solidCase;
};

template <typename T, typename DESCRIPTOR>
struct UnitConverter : public UnitConverterBase<T> {
private:
// Checks to see if parameters include input variables needed for conversion of thermal values
template <typename _T, typename _DESCRIPTOR>
bool thermalConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
 return params.template provides<parameters::PHYS_THERMAL_CONDUCTIVITY>() &&
        params.template provides<parameters::PHYS_HEAT_CAPACITY>() &&
        params.template provides<parameters::PHYS_THERMAL_EXPANSION>() &&
        params.template provides<parameters::PHYS_CHAR_HIGH_TEMPERATURE>() &&
        params.template provides<parameters::PHYS_CHAR_LOW_TEMPERATURE>();
}

template <typename _T, typename _DESCRIPTOR>
void setThermalConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
 const T physThermalConductivity = params.template get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
 const T physSpecificHeatCapacity = params.template get<parameters::PHYS_HEAT_CAPACITY>();
 const T physThermalExpansionCoefficient = params.template get<parameters::PHYS_THERMAL_EXPANSION>();
 const T charPhysHighTemperature = params.template get<parameters::PHYS_CHAR_HIGH_TEMPERATURE>();
 const T charPhysLowTemperature = params.template get<parameters::PHYS_CHAR_LOW_TEMPERATURE>();

 this->_conversionTemperature = charPhysHighTemperature - charPhysLowTemperature;
 this->_conversionThermalDiffusivity = this->_conversionViscosity;
 this->_conversionSpecificHeatCapacity = this->_conversionVelocity * this->_conversionVelocity / this->_conversionTemperature;
 this->_conversionThermalConductivity = this->_conversionForce / this->_conversionTime / this->_conversionTemperature;
 this->_conversionHeatFlux = this->_conversionMass / util::pow(this->_conversionTime.value(), 3);
 this->_charPhysLowTemperature = charPhysLowTemperature;
 this->_charPhysHighTemperature = charPhysHighTemperature;
 this->_charPhysTemperatureDifference = charPhysHighTemperature - charPhysLowTemperature;
 this->_physThermalExpansionCoefficient = physThermalExpansionCoefficient;
 this->_physThermalDiffusivity = physThermalConductivity / (this->_physDensity * physSpecificHeatCapacity);
 this->_physSpecificHeatCapacity = physSpecificHeatCapacity;
 this->_physThermalConductivity = physThermalConductivity;
 this->_latticeThermalRelaxationTime = (this->_physThermalDiffusivity / this->_conversionThermalDiffusivity * descriptors::invCs2<T,DESCRIPTOR>()) + 0.5;
}

// Checks to see if parameters include input variables needed for conversion of ADE values
template <typename _T, typename _DESCRIPTOR>
bool adeConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
 return params.template provides<parameters::PHYS_DIFFUSIVITY>();
}

template <typename _T, typename _DESCRIPTOR>
void setAdeConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
 const T physDiffusivity = params.template get<parameters::PHYS_DIFFUSIVITY>();

 this->_physDiffusivity = physDiffusivity;
 this->_conversionDiffusivity = this->_conversionLength * this->_conversionLength / this->_conversionTime;
 this->_latticeAdeRelaxationTime = (physDiffusivity / this->_conversionDiffusivity * descriptors::invCs2<T,DESCRIPTOR>()) + 0.5;
}

// Checks to see if parameters include input variables needed for conversion of multiphase values
template <typename _T, typename _DESCRIPTOR>
bool multiphaseConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
 return params.template provides<parameters::PHYS_EOS_A>() &&
        params.template provides<parameters::PHYS_EOS_B>() &&
        params.template provides<parameters::PHYS_MOLAR_MASS>() &&
        params.template provides<parameters::SURFACE_TENSION>() &&
        params.template provides<parameters::PHYS_CHAR_TEMPERATURE>();
}

template <typename _T, typename _DESCRIPTOR>
void setMultiphaseConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
 const T physEoSa = params.template get<parameters::PHYS_EOS_A>();
 const T physEoSb = params.template get<parameters::PHYS_EOS_B>();
 const T physMolarMass = params.template get<parameters::PHYS_MOLAR_MASS>();
 const T physSurfaceTension = params.template get<parameters::SURFACE_TENSION>();
 const T charPhysTemperature = params.template get<parameters::PHYS_CHAR_TEMPERATURE>();
 const T latticeEoSa = params.template getOrFallback<parameters::LATTICE_EOS_A>(T{1.0});

 this->_conversionEoSa = physEoSa / latticeEoSa;
 this->_conversionEoSb = physEoSb * T{10.5};
 this->_conversionMolarMass = physMolarMass;
 this->_conversionSurfaceTension = this->_conversionLength * this->_conversionEoSa * util::pow(this->_conversionEoSb.value(), -2);
 this->_conversionTemperature = this->_conversionEoSa / this->_conversionEoSb / this->_conversionGasConstant;
 this->_conversionChemicalPotential = this->_conversionVelocity * this->_conversionVelocity;
 this->_physEoSav = physEoSa;
 this->_physEoSbv = physEoSb;
 this->_physMolarMass = physMolarMass;
 this->_physSurfaceTension = physSurfaceTension;
 this->_charPhysTemperature = charPhysTemperature;
 this->_latticeSurfaceTension = physSurfaceTension / this->_conversionSurfaceTension;
}

// Checks to see if parameters include input variables needed for conversion of power law values
template <typename _T, typename _DESCRIPTOR>
bool powerLawConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
 return params.template provides<parameters::PHYS_CONSISTENCY_COEFF>() &&
        params.template provides<parameters::POWER_LAW_INDEX>();
}

template <typename _T, typename _DESCRIPTOR>
void setPowerLawConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
 const T physConsistencyCoeff = params.template get<parameters::PHYS_CONSISTENCY_COEFF>();
 const T powerLawIndex = params.template get<parameters::POWER_LAW_INDEX>();

 this->_conversionConsistencyCoeff = util::pow(this->_conversionTime.value(), powerLawIndex - 2) * util::pow(this->_conversionLength.value(), 2);
 this->_powerLawIndex = powerLawIndex;
 this->_physConsistencyCoeff = physConsistencyCoeff;
}

// Checks to see if parameters include input variables needed for conversion of radiative values
template <typename _T, typename _DESCRIPTOR>
bool radiativeConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
 return params.template provides<parameters::PHYS_ABSORPTION>() &&
        params.template provides<parameters::PHYS_SCATTERING>();
}

template <typename _T, typename _DESCRIPTOR>
void setRadiativeConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
 const T physAbsorption = params.template get<parameters::PHYS_ABSORPTION>();
 const T physScattering = params.template get<parameters::PHYS_SCATTERING>();
 const T anisotropyFactor = params.template getOrFallback<parameters::ANISOTROPY_FACTOR>(T{0});
 const T refractiveRelative = params.template getOrFallback<parameters::REFRACTIVE_RELATIVE>(T{1});

 this->_physAbsorption = physAbsorption;
 this->_physScattering = physScattering;
 this->_anisotropyFactor = anisotropyFactor;
 this->_extinction = physAbsorption + physScattering;
 this->_scatteringAlbedo = physScattering / (physAbsorption + physScattering);
 this->_physDiffusion = T{1.0} / (T{3.0} * (physAbsorption + physScattering));
 this->_effectiveAttenuation = util::sqrt(T{3.0} * physAbsorption * (physAbsorption + physScattering));
 this->_refractiveRelative = refractiveRelative;
 this->_latticeAbsorption = physAbsorption * this->_conversionLength;
 this->_latticeScattering = physScattering * this->_conversionLength;
 this->_latticeDiffusion = this->_physDiffusion / this->_conversionLength;
}

// ---------------------------------------------------------------------------
// Common linear-elastic material parameters (StructuralMechanicsUnitConverter)
// ---------------------------------------------------------------------------

template <typename _T, typename _DESCRIPTOR>
bool linearElasticConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
  return params.template provides<parameters::YOUNGS_MODULUS>()         &&
         params.template provides<parameters::POISSON_RATIO>()          &&
         params.template provides<parameters::PHYS_CHAR_TIME>()         &&
         params.template provides<parameters::PHYS_CHAR_DISPLACEMENT>() &&
         params.template provides<parameters::SOLID_CASE>();
}

SolidCase solidCaseFromString(const std::string& name)
{
  if (name == "planeStrain") { return SolidCase::planeStrain; }
  if (name == "planeStress") { return SolidCase::planeStress; }
  if (name == "continua")    { return SolidCase::continua;    }
  if (name == "standard")    { return SolidCase::standard;    }
  throw std::runtime_error("Unknown solid case: " + name);
}

template <typename _T, typename _DESCRIPTOR>
void setLinearElasticConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
  const T         youngsModulus       = params.template get<parameters::YOUNGS_MODULUS>();
  const T         poissonRatio        = params.template get<parameters::POISSON_RATIO>();
  const T         charPhysDisplacement= params.template get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const SolidCase   solidCase         = solidCaseFromString(params.template get<parameters::SOLID_CASE>());

  // physical material state
  this->_youngsModulus        = youngsModulus;
  this->_poissonRatio         = poissonRatio;
  this->_solidCase            = solidCase;
  this->_charPhysDisplacement = charPhysDisplacement;

  // Lamé parameters per solid case (same convention as the gradient converter)
  switch (solidCase) {
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
  this->_lameMu           = youngsModulus / (T{2.0} * (T{1.0} + poissonRatio));
  this->_bulkModulus  = this->_lameLambda + T{2.0} / T(_DESCRIPTOR::d) * this->_lameMu;

}

// ---------------------------------------------------------------------------
// Vector-based elastodynamics LBM (Boolakee et al., 2025)
// ---------------------------------------------------------------------------

template <typename _T, typename _DESCRIPTOR>
bool elastodynamicsVectorLBMConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
  return linearElasticConversionInputsProvided(params) &&
         params.template provides<parameters::PHYS_CHAR_VELOCITY>() &&
         params.template provides<parameters::PHYS_CHAR_TIME>();
}

template <typename _T, typename _DESCRIPTOR>
void setElastodynamicsVectorLBMConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
  setLinearElasticConversionFactors(params);

  const T dx  = this->_conversionLength;
  const T dt  = this->_conversionTime;
  const T rho = this->_physDensity;

  // in the vector formulation the characteristic displacement slot carries the velocity
  this->_charPhysVelocity       = params.template get<parameters::PHYS_CHAR_VELOCITY>();
  this->_charPhysTime           = params.template get<parameters::PHYS_CHAR_TIME>();
  this->_charPhysDisplacement   = this->_charPhysVelocity;

  this->_latticeSpeed           =  dx / dt;
  this->_conversionVelocity     =  dx / dt;
  this->_conversionDisplacement =  dx / dt;
  this->_conversionSolidStress  = -dx / dt * rho;
  this->_conversionSolidForce   =  dx / dt * rho / dt;
  this->_conversionTraction     =  dx * dx * rho / (dt * dt);
  this->_conversionLame         =  dx * dx * rho / (dt * dt);

  this->_cMu                    = util::sqrt(this->_lameMu / rho);
  this->_cK                     = util::sqrt(this->_bulkModulus  / rho);
}

// ---------------------------------------------------------------------------
// Quasistationary MRT LBM (Boolakee et al., 2023)
// ---------------------------------------------------------------------------

template <typename _T, typename _DESCRIPTOR>
bool stationaryMRTConversionInputsProvided(const ParametersD<_T, _DESCRIPTOR>& params) const {
  return linearElasticConversionInputsProvided(params) &&
         params.template provides<parameters::DAMPING_FACTOR>();
}

template <typename _T, typename _DESCRIPTOR>
void setStationaryMRTConversionFactors(const ParametersD<_T, _DESCRIPTOR>& params) {
  setLinearElasticConversionFactors(params);

  const T dx            = this->_conversionLength;
  const T dt            = this->_conversionTime;
  const T dampingFactor = params.template get<parameters::DAMPING_FACTOR>();

  this->_conversionDisplacement = dx;
  this->_conversionSolidStress  = dampingFactor * dx * dx / dt;
  this->_conversionSolidForce   = dampingFactor * dx / dt;
  this->_conversionTraction     = dampingFactor * dx * dx / dt;
  this->_conversionLame         = dampingFactor * dx * dx / dt;
  this->_dampingFactor          = dampingFactor;
  this->_epsilon = this->_conversionLength / this->_charPhysLength;
}

public:

  /** Documentation of constructor:
    *  \param physDeltaX              spacing between two lattice cells in __m__
    *  \param physDeltaT              time step in __s__
    *  \param charPhysLength          reference/characteristic length of simulation geometry in __m__
    *  \param charPhysVelocity        maximal or highest expected velocity during simulation in __m / s__
    *  \param physViscosity           physical kinematic viscosity in __m^2 / s__
    *  \param physDensity             physical density in __kg / m^3__
    *  \param charPhysPressure        reference/characteristic physical pressure in Pa = kg / m s^2
    */
  UnitConverter( T physDeltaX, T physDeltaT, T charPhysLength, T charPhysVelocity,
                           T physViscosity, T physDensity, T charPhysPressure = 0 )
  {
    this->_invCs2 = descriptors::invCs2<T,DESCRIPTOR>();
    this->_conversionLength = physDeltaX;
    this->_conversionTime = physDeltaT,
    this->_conversionVelocity = this->_conversionLength / this->_conversionTime;
    this->_conversionDensity = physDensity;
    this->_conversionMass = this->_conversionDensity * util::pow(this->_conversionLength.value(), 3);
    this->_conversionViscosity = this->_conversionLength * this->_conversionLength / this->_conversionTime;
    this->_conversionForce = this->_conversionMass * this->_conversionLength / (this->_conversionTime * this->_conversionTime);
    this->_conversionTorque = this->_conversionMass * this->_conversionLength * this->_conversionLength / (this->_conversionTime * this->_conversionTime);
    this->_conversionPressure = this->_conversionForce / util::pow(this->_conversionLength.value(), 2);
    this->_charPhysLength = charPhysLength;
    this->_charPhysVelocity = charPhysVelocity;
    this->_physViscosity = physViscosity;
    this->_physDensity = physDensity;
    this->_charPhysPressure = charPhysPressure;
    this->_resolution = (std::size_t)(this->_charPhysLength / this->_conversionLength + 0.5);
    this->_latticeRelaxationTime = (this->_physViscosity / this->_conversionViscosity * descriptors::invCs2<T,DESCRIPTOR>()) + 0.5;
    this->_charLatticeVelocity = this->_charPhysVelocity / this->_conversionVelocity;
  }

  template <concepts::LatticeDescriptor _DESCRIPTOR>
  UnitConverter(const UnitConverter<T,_DESCRIPTOR>& rhs) {
    static_cast<UnitConverterBase<T>&>(*this) = static_cast<const UnitConverterBase<T>&>(rhs);
  }

  /** Unified constructor accepting named parameters via ParametersD.
    * Automatically determines initialization strategy based on provided parameters.
    *
    * Required physical characteristic parameters:
    *   - PHYS_CHAR_LENGTH
    *   - PHYS_CHAR_VELOCITY
    *   - PHYS_VISCOSITY
    *   - PHYS_DENSITY
    *
    * Discretization parameters (at least 2 required):
    *   - PHYS_DELTA_X
    *   - PHYS_DELTA_T
    *   - RESOLUTION
    *   - LATTICE_RELAXATION_TIME
    *   - LATTICE_CHAR_VELOCITY
    *
    * Optional:
    *   - PHYS_CHAR_PRESSURE (default: 0)
    *   - Thermal, Multiphase, ADE, PowerLaw, Radiative  unit converter parameters
    *
    * Computed values are written back to the ParametersD object.
    */
  template <typename _T, typename _DESCRIPTOR>
  UnitConverter(const ParametersD<_T, _DESCRIPTOR>& params) {
    // Validate required physical characteristic parameters are explicitly provided
    if (!params.template provides<parameters::PHYS_CHAR_LENGTH>()) {
      throw std::invalid_argument("UnitConverter: Required parameter PHYS_CHAR_LENGTH not provided");
    }
    if (!params.template provides<parameters::PHYS_CHAR_VELOCITY>()) {
      throw std::invalid_argument("UnitConverter: Required parameter PHYS_CHAR_VELOCITY not provided");
    }
    if (!params.template provides<parameters::PHYS_CHAR_VISCOSITY>()) {
      throw std::invalid_argument("UnitConverter: Required parameter PHYS_CHAR_VISCOSITY not provided");
    }
    if (!params.template provides<parameters::PHYS_CHAR_DENSITY>()) {
      throw std::invalid_argument("UnitConverter: Required parameter PHYS_CHAR_DENSITY not provided");
    }

    // Extract required physical characteristic parameters
    T charPhysLength = params.template get<parameters::PHYS_CHAR_LENGTH>();
    T charPhysVelocity = params.template get<parameters::PHYS_CHAR_VELOCITY>();
    T physViscosity = params.template get<parameters::PHYS_CHAR_VISCOSITY>();
    T physDensity = params.template get<parameters::PHYS_CHAR_DENSITY>();
    // PHYS_CHAR_PRESSURE defaults to 0 if not provided
    T charPhysPressure = params.template getOrFallback<parameters::PHYS_CHAR_PRESSURE>(T{0});

    // Initialize invCs2 first
    this->_invCs2 = descriptors::invCs2<T,DESCRIPTOR>();

    // Store physical parameters
    this->_charPhysLength = charPhysLength;
    this->_charPhysVelocity = charPhysVelocity;
    this->_physViscosity = physViscosity;
    this->_physDensity = physDensity;
    this->_charPhysPressure = charPhysPressure;

    // Determine which discretization parameters are explicitly provided
    bool hasPhysDeltaX = params.template provides<parameters::PHYS_DELTA_X>();
    bool hasPhysDeltaT = params.template provides<parameters::PHYS_DELTA_T>();
    bool hasResolution = params.template provides<parameters::RESOLUTION>();
    bool hasTau = params.template provides<parameters::LATTICE_RELAXATION_TIME>();
    bool hasLatticeU = params.template provides<parameters::LATTICE_CHAR_VELOCITY>();

    // Count provided discretization parameters
    int discretizationCount = (hasPhysDeltaX ? 1 : 0) + (hasPhysDeltaT ? 1 : 0) +
                              (hasResolution ? 1 : 0) + (hasTau ? 1 : 0) + (hasLatticeU ? 1 : 0);

    if (discretizationCount < 2) {
      throw std::invalid_argument("UnitConverter: At least 2 discretization parameters required, got " + std::to_string(discretizationCount));
    }

    // Compute discretization parameters based on provided combination
    T physDeltaX{}, physDeltaT{}, latticeRelaxationTime{}, charLatticeVelocity{};

    // Case A: PHYS_DELTA_X + PHYS_DELTA_T (direct initialization)
    if (hasPhysDeltaX && hasPhysDeltaT) {
      physDeltaX = params.template getOrFallback<parameters::PHYS_DELTA_X>(T{});
      physDeltaT = params.template getOrFallback<parameters::PHYS_DELTA_T>(T{});

      if (hasTau) {
        latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});
        // Validate consistency
        T expectedDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
        if (std::abs(expectedDeltaT - physDeltaT) > std::max(T{1e-10}, T{1e-6} * physDeltaT)) {
          throw std::runtime_error("UnitConverter: Conflicting PHYS_DELTA_T and LATTICE_RELAXATION_TIME");
        }
      } else {
        // Compute tau from physDeltaX and physDeltaT
        latticeRelaxationTime = physViscosity / (physDeltaX * physDeltaX) * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT + 0.5;
      }

      if (hasLatticeU) {
        charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});
        // Validate consistency
        T expectedDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
        if (std::abs(expectedDeltaT - physDeltaT) > std::max(T{1e-10}, T{1e-6} * physDeltaT)) {
          throw std::runtime_error("UnitConverter: Conflicting PHYS_DELTA_T and LATTICE_CHAR_VELOCITY");
        }
      } else {
        // Compute lattice velocity from physDeltaX and physDeltaT
        charLatticeVelocity = charPhysVelocity / physDeltaX * physDeltaT;
      }
    }
    // Case B: PHYS_DELTA_X + RESOLUTION
    else if (hasPhysDeltaX && hasResolution) {
      physDeltaX = params.template getOrFallback<parameters::PHYS_DELTA_X>(T{});
      int resolution = params.template getOrFallback<parameters::RESOLUTION>(int{});

      // Validate consistency
      T expectedDeltaX = charPhysLength / resolution;
      if (std::abs(expectedDeltaX - physDeltaX) > std::max(T{1e-10}, T{1e-6} * physDeltaX)) {
        throw std::runtime_error("UnitConverter: Conflicting PHYS_DELTA_X and RESOLUTION");
      }

      if (hasTau) {
        latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});
        physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
      } else if (hasLatticeU) {
        charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});
        physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
        // Also compute tau
        latticeRelaxationTime = physViscosity / (physDeltaX * physDeltaX) * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT + 0.5;
      } else {
        throw std::invalid_argument("UnitConverter: PHYS_DELTA_X + RESOLUTION requires LATTICE_RELAXATION_TIME or LATTICE_CHAR_VELOCITY");
      }
    }
    // Case C: PHYS_DELTA_X + LATTICE_RELAXATION_TIME
    else if (hasPhysDeltaX && hasTau) {
      physDeltaX = params.template getOrFallback<parameters::PHYS_DELTA_X>(T{});
      latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});

      physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;

      if (hasLatticeU) {
        charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});
        // Validate consistency
        T expectedDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
        if (std::abs(expectedDeltaT - physDeltaT) > std::max(T{1e-10}, T{1e-6} * physDeltaT)) {
          throw std::runtime_error("UnitConverter: Conflicting LATTICE_RELAXATION_TIME and LATTICE_CHAR_VELOCITY");
        }
      } else {
        charLatticeVelocity = charPhysVelocity / physDeltaX * physDeltaT;
      }
    }
    // Case D: PHYS_DELTA_X + LATTICE_CHAR_VELOCITY
    else if (hasPhysDeltaX && hasLatticeU) {
      physDeltaX = params.template getOrFallback<parameters::PHYS_DELTA_X>(T{});
      charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});

      physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;

      latticeRelaxationTime = physViscosity / (physDeltaX * physDeltaX) * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT + 0.5;
    }
    // Case E: RESOLUTION + LATTICE_RELAXATION_TIME
    else if (hasResolution && hasTau) {
      int resolution = params.template getOrFallback<parameters::RESOLUTION>(int{});
      latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});

      physDeltaX = charPhysLength / resolution;

      physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;

      charLatticeVelocity = charPhysVelocity / physDeltaX * physDeltaT;
    }
    // Case F: RESOLUTION + LATTICE_CHAR_VELOCITY
    else if (hasResolution && hasLatticeU) {
      int resolution = params.template getOrFallback<parameters::RESOLUTION>(int{});
      charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});

      physDeltaX = charPhysLength / resolution;

      physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;

      latticeRelaxationTime = physViscosity / (physDeltaX * physDeltaX) * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT + 0.5;
    }
    // Case G: LATTICE_RELAXATION_TIME + LATTICE_CHAR_VELOCITY
    else if (hasTau && hasLatticeU) {
      latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});
      charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});

      physDeltaX = physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5);

      physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
    }
    // Case H: PHYS_DELTA_T + RESOLUTION (less common but valid)
    else if (hasPhysDeltaT && hasResolution) {
      int resolution = params.template getOrFallback<parameters::RESOLUTION>(int{});
      physDeltaT = params.template getOrFallback<parameters::PHYS_DELTA_T>(T{});

      physDeltaX = charPhysLength / resolution;

      if (hasTau) {
        latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});
        // Validate
        T expectedDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
        if (std::abs(expectedDeltaT - physDeltaT) > std::max(T{1e-10}, T{1e-6} * physDeltaT)) {
          throw std::runtime_error("UnitConverter: Conflicting PHYS_DELTA_T and LATTICE_RELAXATION_TIME");
        }
      } else {
        latticeRelaxationTime = physViscosity / (physDeltaX * physDeltaX) * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT + 0.5;
      }

      charLatticeVelocity = charPhysVelocity / physDeltaX * physDeltaT;
    }
    // Case I: PHYS_DELTA_T + LATTICE_RELAXATION_TIME
    else if (hasPhysDeltaT && hasTau) {
      physDeltaT = params.template getOrFallback<parameters::PHYS_DELTA_T>(T{});
      latticeRelaxationTime = params.template getOrFallback<parameters::LATTICE_RELAXATION_TIME>(T{});

      // Solve for physDeltaX from tau formula: tau = nu / dx^2 * invCs2 * dt + 0.5
      // => dx^2 = nu * invCs2 * dt / (tau - 0.5)
      physDeltaX = std::sqrt(physViscosity * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT / (latticeRelaxationTime - 0.5));

      if (hasLatticeU) {
        charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});
        // Validate
        T expectedDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
        if (std::abs(expectedDeltaT - physDeltaT) > std::max(T{1e-10}, T{1e-6} * physDeltaT)) {
          throw std::runtime_error("UnitConverter: Conflicting PHYS_DELTA_T and LATTICE_CHAR_VELOCITY");
        }
      } else {
        charLatticeVelocity = charPhysVelocity / physDeltaX * physDeltaT;
      }
    }
    // Case J: PHYS_DELTA_T + LATTICE_CHAR_VELOCITY
    else if (hasPhysDeltaT && hasLatticeU) {
      physDeltaT = params.template getOrFallback<parameters::PHYS_DELTA_T>(T{});
      charLatticeVelocity = params.template getOrFallback<parameters::LATTICE_CHAR_VELOCITY>(T{});

      physDeltaX = charPhysVelocity / charLatticeVelocity * physDeltaT;

      latticeRelaxationTime = physViscosity / (physDeltaX * physDeltaX) * descriptors::invCs2<T,DESCRIPTOR>() * physDeltaT + 0.5;
    }
    else {
      throw std::invalid_argument("UnitConverter: Unsupported combination of discretization parameters");
    }

    // Initialize all conversion factors using computed values
    this->_conversionLength = physDeltaX;
    this->_conversionTime = physDeltaT;
    this->_conversionVelocity = this->_conversionLength / this->_conversionTime;
    this->_conversionDensity = physDensity;
    this->_conversionMass = this->_conversionDensity * util::pow(this->_conversionLength.value(), 3);
    this->_conversionViscosity = this->_conversionLength * this->_conversionLength / this->_conversionTime;
    this->_conversionForce = this->_conversionMass * this->_conversionLength / (this->_conversionTime * this->_conversionTime);
    this->_conversionTorque = this->_conversionMass * this->_conversionLength * this->_conversionLength / (this->_conversionTime * this->_conversionTime);
    this->_conversionPressure = this->_conversionForce / util::pow(this->_conversionLength.value(), 2);
    this->_resolution = (std::size_t)(this->_charPhysLength / this->_conversionLength + 0.5);
    this->_latticeRelaxationTime = latticeRelaxationTime;
    this->_conversionSurfaceTension = this->_conversionDensity * this->_conversionViscosity * this->_conversionViscosity / this->_conversionLength;
    this->_conversionChemicalPotential = this->_conversionVelocity * this->_conversionVelocity;
    this->_charLatticeVelocity = charLatticeVelocity;

    // Initialize optional parameters required by specific converters
    if (thermalConversionInputsProvided(params)) {
      setThermalConversionFactors(params);
    }
    if (adeConversionInputsProvided(params)) {
      setAdeConversionFactors(params);
    }
    if (multiphaseConversionInputsProvided(params)) {
      setMultiphaseConversionFactors(params);
    }
    if (powerLawConversionInputsProvided(params)) {
      setPowerLawConversionFactors(params);
    }
    if (radiativeConversionInputsProvided(params)) {
      setRadiativeConversionFactors(params);
    }
    if (linearElasticConversionInputsProvided(params)) {
      const bool isElastodynamicsVectorLBM = elastodynamicsVectorLBMConversionInputsProvided(params);
      const bool isStationaryMRT           = stationaryMRTConversionInputsProvided(params);

      if (isElastodynamicsVectorLBM && isStationaryMRT) {
        throw std::runtime_error(
          "Ambiguous structural mechanics parameters: both the elastodynamics vector LBM "
          "(PHYS_CHAR_VELOCITY) and the quasistationary MRT LBM (DAMPING_FACTOR) inputs were "
          "provided. These schemes are mutually exclusive - add only the parameters belonging "
          "to the scheme you intend to use.");
      }

      if (isElastodynamicsVectorLBM) {
        setElastodynamicsVectorLBMConversionFactors(params);
      }
      else if (isStationaryMRT) {
        setStationaryMRTConversionFactors(params);
      }
      else {
        throw std::runtime_error(
          "Could not figure out which Structural Mechanics scheme to use based on the provided parameters. "
          "Therefore, no Unit Converter could be set. "
          "Please provide parameters tailored for the intended LBM scheme.");
      }
    }
  }

  virtual void print(std::ostream& fout) const;
  virtual void write(std::string const& fileName = "unitConverter") const;

  void print() const;

};

template <typename T, typename DESCRIPTOR>
UnitConverter<T,DESCRIPTOR> convectivelyRefineUnitConverter(
  const UnitConverter<T,DESCRIPTOR>& converter,
  unsigned scale = 2)
{
  const T refinementFactor = T{1} / scale;
  return UnitConverter<T,DESCRIPTOR>(
    converter.getPhysDeltaX() * refinementFactor,
    converter.getPhysDeltaT() * refinementFactor,
    converter.getCharPhysLength(),
    converter.getCharPhysVelocity(),
    converter.getPhysViscosity(),
    converter.getPhysDensity(),
    converter.getCharPhysPressure()
  );
}

/// creator function with data given by a XML file
template <typename T, typename DESCRIPTOR>
UnitConverter<T, DESCRIPTOR>* createUnitConverter(XMLreader const& params);

template <typename T, typename DESCRIPTOR>
struct UnitConverterFromResolutionAndRelaxationTime : public UnitConverter<T, DESCRIPTOR> {
  UnitConverterFromResolutionAndRelaxationTime(
    size_t resolution,
    T latticeRelaxationTime,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow((charPhysLength/resolution),2) / physViscosity,
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  { }
};

template <typename T, typename DESCRIPTOR>
struct UnitConverterFromResolutionAndLatticeVelocity : public UnitConverter<T, DESCRIPTOR> {
  UnitConverterFromResolutionAndLatticeVelocity(
    size_t resolution,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (charLatticeVelocity / charPhysVelocity * charPhysLength / resolution),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  { }
};

template <typename T, typename DESCRIPTOR>
struct UnitConverterFromRelaxationTimeAndLatticeVelocity : public UnitConverter<T, DESCRIPTOR> {
  UnitConverterFromRelaxationTimeAndLatticeVelocity(
    T latticeRelaxationTime,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5)),
        (physViscosity * charLatticeVelocity * charLatticeVelocity / charPhysVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5)),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  {
  }
};

template<typename T, typename DESCRIPTOR, typename MEMBER, typename... ARGS>
std::unique_ptr<UnitConverter<T,DESCRIPTOR>> createUnitConverter(ARGS&&... args) {
  return std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
}

}

#include "thermalUnitConverter.h"
#include "multiPhaseUnitConverter.h"
#include "structuralMechanicsUnitConverter.h"
#endif
