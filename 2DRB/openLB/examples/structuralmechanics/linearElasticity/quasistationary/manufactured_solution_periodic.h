#ifndef MANUFACTURED_SOLUTION_PERIODIC_H
#define MANUFACTURED_SOLUTION_PERIODIC_H

namespace olb {

  template <typename T>
  class ForceField2D : public AnalyticalF2D<T, T> {
    private:
      T _lambda;
      T _mu;
      T _conv;

    public:
      ForceField2D(T lambda, T mu, T conv) : AnalyticalF2D<T, T>(2),
        _lambda(lambda),
        _mu(mu),
        _conv(conv) {}

    bool operator()(T output[], const T input[]) override {
      T x = input[0];
      T y = input[1];
      T pi = T(M_PI);
      T pi2 = pi * pi;

      output[0] = (T(4) / T(625) * (_lambda + _mu) * pi2 * util::cos(T(2) * pi * x) * util::sin(T(2) * pi * y) + T(9) / T(1250) * _mu * pi2 * util::cos(T(2) * pi * x) * util::sin(T(2) * pi * y)) / _conv;

      output[1] = (T(4) / T(625) * (_lambda + _mu) * pi2 * util::sin(T(2) * pi * x) * util::cos(T(2) * pi * y) + T(7) / T(1250) * _mu * pi2 * util::sin(T(2) * pi * x) * util::cos(T(2) * pi * y)) / _conv;

      return true;
    }
  };

  template <typename T>
  class ManufacturedSolutionU2D : public AnalyticalF2D<T, T> {
    private:
      bool _do;

    public:
      ManufacturedSolutionU2D(bool run) : AnalyticalF2D<T, T>(2),
      _do(run) {}

      bool operator()(T output[], const T input[]) override {
        T x = input[0];
        T y = input[1];
        T pi = T(M_PI);

        output[0] = T(9) / T(10000) * util::cos(T(2) * pi * x) * util::sin(T(2) * pi * y);
        output[1] = T(7) / T(10000) * util::sin(T(2) * pi * x) * util::cos(T(2) * pi * y);

        return true;
      }
  };

  template <typename T>
  class ManufacturedSolutionStress2D : public AnalyticalF2D<T, T> {
    private:
      T _lambda;
      T _mu;

    public:
      ManufacturedSolutionStress2D(T lambda, T mu) : AnalyticalF2D<T, T>(4),
      _lambda(lambda),
      _mu(mu) {}

      bool operator()(T output[], const T input[]) override {
        T x = input[0];
        T y = input[1];
        T pi = T(M_PI);

        output[0] = -T(2) / T(625) * (_lambda + _mu) * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) - _mu * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) / T(2500);

        output[1] = T(2) / T(625) * _mu * util::cos(T(2) * pi * x) * pi * util::cos(T(2) * pi * y);

        output[2] = T(2) / T(625) * _mu * util::cos(T(2) * pi * x) * pi * util::cos(T(2) * pi * y);

        output[3] = -T(2) / T(625) * (_lambda + _mu) * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) + _mu * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) / T(2500);

        return true;
      }
  };


}

#endif
