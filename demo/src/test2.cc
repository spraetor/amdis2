#include "AMDiS.h"

int main(int argc, char** argv)
{
    using namespace AMDiS;
    
    AMDiS::init(argc, argv);
    
    // Test 1
    {
      auto f = (X(0) > 0.2) && (X(0) < 0.8);
      WorldVector<double> x1; x1 = 0.3;
      std::cout << "f(" << x1 << ") = " << f(x1) << "\n";
      WorldVector<double> x2; x2 = 0.1;
      std::cout << "f(" << x2 << ") = " << f(x2) << "\n";
      WorldVector<double> x3; x3 = 0.9;
      std::cout << "f(" << x3 << ") = " << f(x3) << "\n";
    }
    
    // Test 2
    {
      auto g1 = if_( X(0) > 0.2, X(0), X(0)+1.0 );

      WorldVector<double> x1; x1 = 0.3;
      std::cout << "g1(" << x1 << ") = " << g1(x1) << "\n";
      WorldVector<double> x2; x2 = 0.1;
      std::cout << "g1(" << x2 << ") = " << g1(x2) << "\n";
      
      auto g2 = if_( true, X(0), X(1) );
      WorldVector<double> x3; x3[0] = 1.0; x3[1] = 0.0;
      std::cout << "g2(" << x3 << ") = " << g2(x3) << "\n";
    }
    
    // Test 3
    {
      std::function<bool(WorldVector<double>)> f = (X(0) > 0.2) && (X(0) < 0.8);
      WorldVector<double> x1; x1 = 0.3;
      std::cout << "f(" << x1 << ") = " << f(x1) << "\n";
      
      std::function<bool(WorldVector<double>)> g = f; 
      std::cout << "g(" << x1 << ") = " << g(x1) << "\n";
    }
    
    AMDiS::finalize();
}