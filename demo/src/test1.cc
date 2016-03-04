#include "AMDiS.h"

int main(int argc, char** argv)
{
    using namespace AMDiS;
    
    AMDiS::init(argc, argv);
    
    // Test 1
    {
        WorldVector<double> w1;
        WorldVector<double> w2(Global::getGeo(WORLD));
        WorldVector<double> w3(DEFAULT_SIZE);
        
        WorldVector<double> w4(DEFAULT_SIZE, 4.0);
        
        std::cout << "w4 = " << w4 << "\n";
    }
    
    // Test 2
    {
        int dim = 2;
        
        FixVec<double, VERTEX> w5(dim);
        FixVec<double, VERTEX> w6(dim, 6.0);
        
        std::cout << "w6 = " << w6 << "\n";
        
        FixVec<WorldVector<double>, VERTEX> w7(dim);
        
        WorldVector<double> _w(DEFAULT_SIZE, 8.0);
        FixVec<WorldVector<double>, VERTEX> w8(dim, _w);
        
        std::cout << "w8 = " << w8 << "\n";
    }
    
    // Test 3
    {
        int dim = 2;
        
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> w9(dim);
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> w10(dim, FixVec<WorldVector<double>, VERTEX>(dim));
        
        WorldVector<double> _w0(DEFAULT_SIZE, 11.0);
        FixVec<WorldVector<double>, VERTEX> _w1(dim, _w0);
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> w11(dim, _w1);
        std::cout << "w11 = " << w11 << "\n";
    }
    
    // Test 4
    {
        int dim = 2;
        
        WorldVector<double> _w0(DEFAULT_SIZE, 1.0);
        WorldVector<double> _w1(DEFAULT_SIZE, 2.0);
        FixVec<WorldVector<double>, VERTEX> _w2(dim, _w0);
        FixVec<WorldVector<double>, VERTEX> _w3(dim, _w1);
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> v1(dim, _w2), v2(dim, _w3);
        
        std::cout << "dot(v1, v2) = " << dot(v1, v2) << "\n";
    }
    
    // Test 5
    {
        int dim = 2;
        
        DimVec<double> v1(dim-1);     v1 = 1.0;
        std::cout << "v1 = " << v1 << "\n";
        StaticVector<double, 2> v2; v2 = 2.0;
        std::cout << "v2 = " << v2 << "\n";
        Vector<double> v3(dim);     v3 = 3.0;
        std::cout << "v3 = " << v3 << "\n";
        
        std::cout << "v1*v2 = " << (v1 * v2) << "\n";
        std::cout << "v1*v3 = " << (v1 * v3) << "\n";
        std::cout << "v2*v3 = " << (v2 * v3) << "\n";
        
        FixMat<double, WORLD> m1; m1 = 1.0;
        WorldMatrix<double>   m2; m2 = 2.0;
        DimMat<double>        m3(dim-1); m3 = 3.0;
        StaticMatrix<double, 2, 2> m4(2,2); m4 = 4.0;
        Matrix<double> m5(dim, dim); m5 = 5.0;
        
        std::cout << "m1 = " << m1 << "\n";
        std::cout << "m2 = " << m2 << "\n";
        std::cout << "m3 = " << m3 << "\n";
        std::cout << "m4 = " << m4 << "\n";
        std::cout << "m5 = " << m5 << "\n";
    }
    
    // Test 6
    {
      int dim = 2;
        
      WorldVector<double> v1; v1 = 1.0;
      Vector<double> v2(dim); v2 = 1.0;
      
      std::cout << "v1 == v2: " << (v1 == v2) << "\n";
      std::cout << "v1 != v2: " << (v1 != v2) << "\n";
      std::cout << "v1 < v2: " << (v1 < v2) << "\n";
      std::cout << "v1 > v2: " << (v1 > v2) << "\n";
    }
    
    
    AMDiS::finalize();
}