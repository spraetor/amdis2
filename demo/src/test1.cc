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
        
        std::cout << "w1 = " << w1 << "\n";
        std::cout << "w2 = " << w2 << "\n";
        std::cout << "w3 = " << w3 << "\n";
        
        WorldVector<double> w4(DEFAULT_SIZE, 1.0);
        
        std::cout << "w4 = " << w4 << "\n";
    }
    
    // Test 2
    {
        int dim = 2;
        
        FixVec<double, VERTEX> w5(dim);
        FixVec<double, VERTEX> w6(dim, 0.0);
        
        std::cout << "w5 = " << w5 << "\n";
        std::cout << "w6 = " << w6 << "\n";
        
        FixVec<WorldVector<double>, VERTEX> w7(dim);
        
        std::cout << "w7 = " << w7 << "\n";
        
        WorldVector<double> _w(DEFAULT_SIZE, 2.0);
        FixVec<WorldVector<double>, VERTEX> w8(dim, _w);
        
        std::cout << "w8 = " << w8 << "\n";
    }
    
    // Test 3
    {
        int dim = 2;
        
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> w9(dim);
        //w9.set(0.0);
        std::cout << "w9 = " << w9 << "\n";
        
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> w10(dim, FixVec<WorldVector<double>, VERTEX>(dim));
        //w10.set(0.0);
        
        std::cout << "w10 = " << w10 << "\n";
        
        WorldVector<double> _w0(DEFAULT_SIZE, 0.0);
        FixVec<WorldVector<double>, VERTEX> _w1(dim, _w0);
        FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> w11(dim, _w1);
        
        std::cout << "w11 = " << w11 << "\n";
        
    }

    // neighbourCoord(mesh->getDim(), FixVec<WorldVector<double>, VERTEX>(mesh->getDim())),
    
    
    
    AMDiS::finalize();
}