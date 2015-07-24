#include <vector>

#include "Global.h"
#include "Initfile.h"
#include "Element.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"

namespace AMDiS 
{
  int Global::dimOfWorld = 0;
  std::vector<std::vector<int> > Global::geoIndexTable;

  Element *Global::referenceElement[4] = 
    { NULL, 
      new Line(NULL), 
      new Triangle(NULL), 
      new Tetrahedron(NULL) 
    };

  void Global::init()
  {
    int d = -1;

    // get dimension
    TEST_EXIT(Parameters::initialized())("Parameters not initialized!\n");
    Parameters::get("dimension of world",d,0);
    TEST_EXIT(d > 0)("Cannot initialize dimension!\n");
    TEST_EXIT((d == 1) || (d == 2) || (d == 3))("Invalid world dimension %d!\n",d);

    // set dimension
    dimOfWorld = d;

    // prepare geoIndex-Table
    int geoTableSize = abs(static_cast<int>(MINPART)) + MAXPART + 1;
    geoIndexTable.resize(4);
    for (int i = 0; i < 4; i++) {
      geoIndexTable[i].resize(geoTableSize);
      for (int j = 0; j < geoTableSize; j++)
        geoIndexTable[i][j] = 0;      
    }

    geoIndexTable[0][PARTS - MINPART] = 1;
    geoIndexTable[0][VERTEX - MINPART] = 1;
    geoIndexTable[0][EDGE - MINPART] = 0;
    geoIndexTable[0][FACE - MINPART] = 0;
    geoIndexTable[0][WORLD - MINPART] = dimOfWorld;

    for (int i = 1; i < 4; i++) {
      geoIndexTable[i][CENTER - MINPART] = referenceElement[i]->getGeo(CENTER);
      geoIndexTable[i][VERTEX - MINPART] = referenceElement[i]->getGeo(VERTEX);
      geoIndexTable[i][EDGE - MINPART] = referenceElement[i]->getGeo(EDGE);
      geoIndexTable[i][FACE - MINPART] = referenceElement[i]->getGeo(FACE);
      geoIndexTable[i][DIMEN - MINPART] = referenceElement[i]->getGeo(DIMEN);
      geoIndexTable[i][PARTS - MINPART] = referenceElement[i]->getGeo(PARTS);
      geoIndexTable[i][NEIGH - MINPART] = referenceElement[i]->getGeo(NEIGH);
      geoIndexTable[i][WORLD - MINPART] = dimOfWorld;
      geoIndexTable[i][BOUNDARY - MINPART] = referenceElement[i]->getGeo(BOUNDARY);
      geoIndexTable[i][PROJECTION - MINPART] = referenceElement[i]->getGeo(PROJECTION);
    }

    // set msgWait
    Msg::setMsgWait(!(Parameters::getMsgWait() == 0));
  }


  void Global::clear()
  {
    delete referenceElement[1];
    delete referenceElement[2];
    delete referenceElement[3];
  }

} // end namespace AMDiS
