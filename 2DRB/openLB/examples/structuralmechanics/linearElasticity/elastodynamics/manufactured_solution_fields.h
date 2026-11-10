#ifndef LINEARELASTICITY_ELASTODYNAMICS_MANUFACTURED_SOLUTION_FIELDS_H
#define LINEARELASTICITY_ELASTODYNAMICS_MANUFACTURED_SOLUTION_FIELDS_H

namespace olb {

  template <typename T>
  class AnalyticalForce : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _rho0;
      T _lambda;
      T _mu;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalForce(T ts, T rho0, T lambda, T mu, T convFactor = T(1), T A = T(1), T k = T(2 * M_PI))
          : AnalyticalF2D<T, T>(2),
          _ts(ts),
          _rho0(rho0),
          _lambda(lambda),
          _mu(mu),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override {
        T x = input[0];
        T y = input[1];

        output[0] = (T(8)*_A*(_k)*(_k)*_mu*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y)) - _lambda*(-T(2)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(_k*(-T(0.10000)*_ts + y)) - T(4)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))) - T(2)*_mu*(-T(1.00000000000000)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(_k*(-T(0.10000)*_ts + y)) - T(0.50000)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))) + _rho0*(T(3.20000)*_A*(_k)*(_k)*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::cos(T(2)*_k*(_ts + T(-0.10000))) - T(0.96000)*_A*(_k)*(_k)*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) - T(5.00000000000000)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y)) - T(2.40000)*_A*(_k)*(_k)*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x))*util::cos(T(2)*_k*(_ts + T(-0.10000))))) / _convFactor;
        output[1] = (T(2)*_A*(_k)*(_k)*_mu*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x)) - _lambda*(-T(2)*_A*(_k)*(_k)*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) - _A*(_k)*(_k)*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))) - T(2)*_mu*(-T(1.00000000000000)*_A*(_k)*(_k)*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) - T(2.00000000000000)*_A*(_k)*(_k)*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))) + _rho0*(T(5.60000)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::cos(T(2)*_k*(_ts + T(0.40000))) - T(0.28000)*_A*(_k)*(_k)*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(_k*(-T(0.10000)*_ts + y)) - T(5.97000)*_A*(_k)*(_k)*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x)) - T(0.40000)*_A*(_k)*(_k)*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y))*util::cos(T(2)*_k*(_ts + T(0.40000))))) / _convFactor;

        return true;
      };
  };

  template <typename T>
  class AnalyticalRho : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _rho0;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalRho(T ts, T rho0, T convFactor = T(1), T A = T(1), T k = T(2* M_PI)) : AnalyticalF2D<T, T>(1),
          _ts(ts),
          _rho0(rho0),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override
      {
        T x = input[0];
        T y = input[1];

        output[0] = (_rho0*(-T(2)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) - _A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)) + T(1))) / _convFactor;

        return true;
      };
  };

template <typename T>
class AnalyticalRhoDeviation : public AnalyticalF2D<T,T> {
  AnalyticalRho<T> _anaRho;
  T _rho0;
public:
  AnalyticalRhoDeviation(T ts, T rho0)
    : AnalyticalF2D<T,T>(1),
    _anaRho(ts, rho0, T(1)),
    _rho0(rho0) {}

  bool operator()(T output[], const T input[]) override {
    _anaRho(output, input);

    output[0] -= _rho0;   // -rho0_lat * div u

    return true;
  }
};


  template <typename T>
  class AnalyticalVelocity : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalVelocity(T ts, T convFactor = T(1), T A = T(1), T k = T(2* M_PI)) : AnalyticalF2D<T, T>(2),
          _ts(ts),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override
      {
        T x = input[0];
        T y = input[1];

        output[0] = (T(0.80000)*_A*_k*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000))) + T(2)*_A*_k*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(_ts + T(-0.10000))) - T(0.60000)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x))) / _convFactor;
        output[1] = (T(1.40000)*_A*_k*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000))) + T(2)*_A*_k*util::sin(_k*(-T(0.10000)*_ts + y))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(T(2)*_k*(_ts + T(0.40000))) - T(0.10000)*_A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y))) / _convFactor;

        return true;
      };
  };

  template <typename T>
  class AnalyticalMomentum : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _rho0;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalMomentum(T ts, T rho0, T convFactor = T(1), T A = T(1), T k = T(2 * M_PI))
        : AnalyticalF2D<T, T>(2),
          _ts(ts),
          _rho0(rho0),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override
      {
        T x = input[0];
        T y = input[1];

        output[0] = (_rho0*(T(0.80000)*_A*_k*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000))) + T(2)*_A*_k*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(_ts + T(-0.10000))) - T(0.60000)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)))) / _convFactor;
        output[1] = (_rho0*(T(1.40000)*_A*_k*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000))) + T(2)*_A*_k*util::sin(_k*(-T(0.10000)*_ts + y))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(T(2)*_k*(_ts + T(0.40000))) - T(0.10000)*_A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)))) / _convFactor;

        return true;
      };
  };


  template <typename T>
  class AnalyticalDisplacement : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalDisplacement(T ts, T convFactor = T(1), T A = T(1), T k = T(2 * M_PI))
        : AnalyticalF2D<T, T>(2),
          _ts(ts),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override
      {
        T x = input[0];
        T y = input[1];

        output[0] = (_A*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))) / _convFactor;
        output[1] = (_A*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))) / _convFactor;

        return true;
      };
  };

template <typename T>
  class AnalyticalPoissonStress : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _mu;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalPoissonStress(T ts, T mu, T convFactor = T(1), T A = T(1), T k = T(2 * M_PI))
        : AnalyticalF2D<T, T>(4),
          _ts(ts),
          _mu(mu),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override
      {
        T x = input[0];
        T y = input[1];

        output[0] = (-_mu*(T(6)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) + _A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)))) / _convFactor;
        output[1] = (-_mu*(-_A*_k*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000))) - T(2)*_A*_k*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000))))) / _convFactor;
        output[2] = (-_mu*(-_A*_k*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000))) - T(2)*_A*_k*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000))))) / _convFactor;
        output[3] = (-_mu*(T(2)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) + T(3)*_A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)))) / _convFactor;

        return true;
      };
  };


  template <typename T>
  class AnalyticalCauchyStress : public AnalyticalF2D<T, T> {
    protected:
      T _ts;
      T _lambda;
      T _mu;
      T _convFactor;
      T _A;
      T _k;

    public:
      AnalyticalCauchyStress(T ts, T lambda, T mu, T convFactor = T(1), T A = T(1), T k = T(2 * M_PI))
        : AnalyticalF2D<T, T>(4),
          _ts(ts),
          _lambda(lambda),
          _mu(mu),
          _convFactor(convFactor),
          _A(A),
          _k(k) {};

      bool operator()(T output[], const T input[]) override
      {
        T x = input[0];
        T y = input[1];

        output[0] = (T(4)*_A*_k*_mu*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) + _lambda*(T(2)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) + _A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)))) / _convFactor;
        output[1] = (T(2)*_mu*(-T(0.50000)*_A*_k*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000))) - T(1.00000000000000)*_A*_k*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000))))) / _convFactor;
        output[2] = (T(2)*_mu*(-T(0.50000)*_A*_k*util::sin(_k*(-T(0.80000)*_ts + y))*util::sin(T(2)*_k*(-T(0.30000)*_ts + x))*util::sin(T(2)*_k*(_ts + T(-0.10000))) - T(1.00000000000000)*_A*_k*util::sin(T(2)*_k*(-T(0.70000)*_ts + x))*util::sin(_k*(-T(0.10000)*_ts + y))*util::sin(T(2)*_k*(_ts + T(0.40000))))) / _convFactor;
        output[3] = (T(2)*_A*_k*_mu*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)) + _lambda*(T(2)*_A*_k*util::sin(T(2)*_k*(_ts + T(-0.10000)))*util::cos(_k*(-T(0.80000)*_ts + y))*util::cos(T(2)*_k*(-T(0.30000)*_ts + x)) + _A*_k*util::sin(T(2)*_k*(_ts + T(0.40000)))*util::cos(T(2)*_k*(-T(0.70000)*_ts + x))*util::cos(_k*(-T(0.10000)*_ts + y)))) / _convFactor;

        return true;
      };
  };

  template <typename T>
  class AnalyticalGradRho : public AnalyticalF2D<T,T> {
    T _ts, _rho0, _mu, _conv, _A, _k, _h;
  public:
    AnalyticalGradRho(T ts, T rho0, T mu, T conv = T(1),
                      T A = T(1), T k = T(2*M_PI), T h = T(1e-6))
      : AnalyticalF2D<T,T>(2), _ts(ts), _rho0(rho0), _mu(mu),
        _conv(conv), _A(A), _k(k), _h(h) {}

    T rhoDev(T x, T y) {
      AnalyticalPoissonStress<T> pi(_ts, _mu, T(1), _A, _k);
      T out[4]; T in[2] = {x, y};
      pi(out, in);
      return _rho0 * (out[0] + out[3]) / (T(4) * _mu);
    }

    bool operator()(T output[], const T input[]) override {
      const T x = input[0], y = input[1];
      output[0] = (rhoDev(x + _h, y) - rhoDev(x - _h, y)) / (T(2) * _h) / _conv;
      output[1] = (rhoDev(x, y + _h) - rhoDev(x, y - _h)) / (T(2) * _h) / _conv;
      return true;
    }
  };

} // namespace olb

#endif
