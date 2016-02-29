#include "AMDiS.h"
#include "expressions/MatrixVectorTerms.hpp"

int main(int argc, char** argv)
{
    using namespace AMDiS;
    
    AMDiS::init(argc, argv);
    
    // Test 1 (boolean expressions)
    {
      auto f1 = (X(0) > 0.2) && (X(0) > 0.8);
      auto f2 = (X(0) == X(1)) || (X(0) != 0.8);
      
      WorldVector<double> x1 = {0.3, 0.7};
      WorldVector<double> x2{0.1, 0.1};
      WorldVector<double> x3; x3 = {0.2, 0.9};
      
      std::cout << "f1(" << x1 << ") = " << f1(x1) << "\n";
      std::cout << "f1(" << x2 << ") = " << f1(x2) << "\n";
      std::cout << "f1(" << x3 << ") = " << f1(x3) << "\n";
      
      std::cout << "f2(" << x1 << ") = " << f2(x1) << "\n";
      std::cout << "f2(" << x2 << ") = " << f2(x2) << "\n";
      std::cout << "f2(" << x3 << ") = " << f2(x3) << "\n";
    }
    
    // Test 2 (conditonals in expressions)
    {
      auto g1 = if_( X(0) > 0.2, X(0), X(0)+1.0 );

      WorldVector<double> x1; x1 = 0.3;
      WorldVector<double> x2; x2 = 0.1;
      
      std::cout << "g1(" << x1 << ") = " << g1(x1) << "\n";
      std::cout << "g1(" << x2 << ") = " << g1(x2) << "\n";
      
      auto g2 = if_( true, X(0), X(1) );
      
      WorldVector<double> x3; x3[0] = 1.0; x3[1] = 0.0;
      
      std::cout << "g2(" << x3 << ") = " << g2(x3) << "\n";
    }
    
    // Test 3 (assign an expression-term to std::function)
    {
      std::function<bool(WorldVector<double>)> f = (X(0) > 0.2) && (X(0) < 0.8);
      WorldVector<double> x1; x1 = 0.3;
      std::cout << "f(" << x1 << ") = " << f(x1) << "\n";
      
      std::function<bool(WorldVector<double>)> g = f; 
      std::cout << "g(" << x1 << ") = " << g(x1) << "\n";
    }
    
    // Test 4
    {
      using V = Vector<double>;
      V v{0.1, 0.2};
      auto f = 2.0*X();
      auto g = two_norm(2.0*X());
      auto h = f[0] + f[1];
      
      
      WorldVector<double> x{1.0, 1.0};
      // V y = f(x);
      // std::cout << "f(x) = " << y << "\n";
      std::cout << "g(x) = " << g(x) << "\n";
      std::cout << "h(x) = " << h(x) << "\n";
    }
    
    
    AMDiS::finalize();
}