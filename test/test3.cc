#include "MatrixVector.hpp"
#include "matrix_vector/MatrixVectorOperations.hpp"
#include "expressions/MathOperations.hpp"
// #include "expressions/MatrixVectorTerms.hpp"

using namespace AMDiS;

template <class T0, class T1 = T0>
struct Mult : FunctorBase
{
  static constexpr auto eval(T0 const& v0, T1 const& v1) -> decltype(v0 * v1) { return v0 * v1; }
//   constexpr auto operator() (T0 const& v0, T1 const& v1) const RETURNS_CONST_REF( eval(v0, v1) )
};

int main(int argc, char** argv)
{
    {
      using V = WorldVector<double>;
      std::cout << "1)\n";
      V v{1.0, 2.0, 3.0};

//       std::cout << "2)\n";
//       using F = Mult<double, V>;
//
//       std::cout << "3)\n";
//       auto&& erg = F::eval(2.0, v);
//       std::cout << "4)\n";
//       V w = erg;
//
//       std::cout << "5)\n";
//       std::cout << "v: " << size(v) << "\n";
//       std::cout << "6)\n";
//       std::cout << "erg: " << size(erg) << "\n";
//       std::cout << "7)\n";
//       std::cout << "w: " << size(w) << "\n";
//
// //       print_type( erg );
//
      std::cout << "8)\n";
      auto&& g = constant(2.0) * X();
//       std::cout << "9)\n";
//       std::cout << "g(v): " << size(g(v)) << "\n";
//
//       std::cout << "10)\n";
//       std::cout << "g(w): " << size(g(w)) << "\n";

      std::cout << "11)\n";
      print_type(g);
      print_type(v);
      V w2 = g(v);
//       std::cout << "12)\n";
//       V w3 = g(w);
//       print_type( g(v) );
    }
}