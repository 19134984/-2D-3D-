#ifndef MANUFACTURED_SOLUTION_BOUNDARIES_H
#define MANUFACTURED_SOLUTION_BOUNDARIES_H

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

    bool operator()(T output[], const T input[]) override
    {
      T x = input[0];
      T y = input[1];
      T pi = T(M_PI);
      T pi2 = pi * pi;

      output[0] = (-(_lambda + _mu) * (-T(9) / T(2500) * pi2 * (y * y) * util::sin(T(2) * pi * x * y) - T(7) / T(2500) * pi * util::sin(T(2) * pi * y) * x) - _mu * (-T(9) / T(2500) * pi2 * (y * y) * util::sin(T(2) * pi * x * y) - T(9) / T(2500) * pi2 * (x * x) * util::sin(T(2) * pi * x * y))) / _conv;

      output[1] = (-(_lambda + _mu) * (T(9) / T(5000) * pi * util::cos(T(2) * pi * x * y) - T(9) / T(2500) * pi2 * y * x * util::sin(T(2) * pi * x * y) - T(7) / T(2500) * pi2 * util::cos(T(2) * pi * y) * (x * x + T(1))) - _mu * (T(7) / T(5000) * util::cos(T(2) * pi * y) - T(7) / T(2500) * pi2 * util::cos(T(2) * pi * y) * (x * x + T(1)))) / _conv;

      return true;
    }
  };

  template<typename T>
  class ManufacturedSolutionU2D : public AnalyticalF2D<T, T> {
    protected:
      bool _do;

    public:
      ManufacturedSolutionU2D(bool run) : AnalyticalF2D<T, T>(2),
      _do(run) {}

      bool operator()(T output[], const T input[]) override {
        T x = input[0];
        T y = input[1];
        T pi = T(M_PI);

        output[0] = T(9) / T(10000) * util::sin(T(2) * pi * x * y);
        output[1] = T(7) / T(10000) * util::cos(T(2) * pi * y) * (x * x + T(1));

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

        output[0] = (_lambda + _mu) * (T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y) - T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1))) + _mu * (T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y) + T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1)));

        output[1] = _mu * (T(7) / T(5000) * util::cos(T(2) * pi * y) * x + T(9) / T(5000) * pi * x * util::cos(T(2) * pi * x * y));

        output[2] = _mu * (T(7) / T(5000) * util::cos(T(2) * pi * y) * x + T(9) / T(5000) * pi * x * util::cos(T(2) * pi * x * y));

        output[3] = (_lambda + _mu) * (T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y) - T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1))) + _mu * (-T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1)) - T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y));

        return true;
      }
  };

} // namespace olb

#endif
