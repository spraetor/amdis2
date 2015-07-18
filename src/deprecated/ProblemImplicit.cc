/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


#include "ProblemImplicit.h"
#include "MathFunctions.h"
#include "io/Reader.h"
using namespace std;

namespace AMDiS {

  void ProblemImplicit::readDofVec(const std::string& filename, DOFVector<double>* vec, 
                                       Mesh* mesh) 
  {
#if 0
    long size = mesh->getNumberOfVertices();

    TEST_EXIT(vec->getSize()>=size)("dofvector smaller than vertex data vector\n");
    std::string buf;
    getline(in, buf);
    while (!in.eof() && buf != "vertex values: solution") 
      getline(in, buf);

    TEST_EXIT(!in.eof())("no vertex data\n");

    for (long i = 0; i < size ; ++i) {
      in >> buf;
      (*vec)[i] = atof(buf.c_str());
    }
#endif
    io::readFile(filename, vec);
  }


  void ProblemImplicit::readR(const std::string& inStream, double eps, 
	                          Mesh* mesh, int implMesh , int comp) 
  {
    FUNCNAME("ProblemImplicit::readR()");

    DOFVector<double>* r = getSignedDistance(implMesh, comp);
    DOFVector<double>* phi1 = getPhi1(implMesh, comp);
    DOFVector<double>* phi2 = getPhi2(implMesh, comp);
    DOFVector<double>* levelSet = getLevelset(implMesh, comp);

    TEST_EXIT(r != NULL)("no signed distance vector\n");
    TEST_EXIT(phi1 != NULL)("no phasefield1 vector\n");
    TEST_EXIT(phi2 != NULL)("no phasefield2 vector\n");
    TEST_EXIT(levelSet != NULL)("no levelSet vector\n");

    bool checkSize = r->getSize() == phi1->getSize() && 
      r->getSize() == phi2->getSize();

    TEST_EXIT(checkSize)("something went wrong\n");
    
    readDofVec(inStream, r, mesh);

    for (int i = r->getSize() - 1; i >= 0; --i) {
      (*phi1)[i] = Phi1((*r)[i], eps);
      (*phi2)[i] = Phi2((*r)[i], eps);
      (*levelSet)[i] = LevelSet((*r)[i]);
    }
  }


  void ProblemImplicit::readPhi1(const std::string& inStream, double eps, 
                                     Mesh* mesh, int implMesh, int comp)
  {
    DOFVector<double>* r = getSignedDistance(implMesh, comp);
    DOFVector<double>* phi1 = getPhi1(implMesh, comp);
    DOFVector<double>* phi2 = getPhi2(implMesh, comp);
    DOFVector<double>* levelSet = getLevelset(implMesh, comp);

    TEST_EXIT(r != NULL)("no signed distance vector\n");
    TEST_EXIT(phi1 != NULL)("no phasefield1 vector\n");
    TEST_EXIT(phi2 != NULL)("no phasefield2 vector\n");
    TEST_EXIT(levelSet != NULL)("no levelSet vector\n");

    bool checkSize = r->getSize() == phi1->getSize() && 
      r->getSize() == phi2->getSize();

    TEST_EXIT(checkSize)("something went wrong\n");

    readDofVec(inStream, phi1, mesh);

    for (int i = r->getSize() - 1; i>= 0; --i) {
      (*r)[i] = Phi1ToR((*phi1)[i], eps);
      (*phi2)[i] = 1 - (*phi1)[i];
      (*levelSet)[i] = LevelSet((*r)[i]);
    }
  }


  void ProblemImplicit::readPhi2(const std::string& inStream, double eps, 
	                             Mesh* mesh, int implMesh, int comp)
  { 
    DOFVector<double>* r = getSignedDistance(implMesh, comp);
    DOFVector<double>* phi1 = getPhi1(implMesh, comp);
    DOFVector<double>* phi2 = getPhi2(implMesh, comp);
    DOFVector<double>* levelSet = getLevelset(implMesh, comp);

    TEST_EXIT(r != NULL)("no signed distance vector\n");
    TEST_EXIT(phi1 != NULL)("no phasefield1 vector\n");
    TEST_EXIT(phi2 != NULL)("no phasefield2 vector\n");
    TEST_EXIT(levelSet != NULL)("no levelSet vector\n");

    bool checkSize = r->getSize() == phi1->getSize() &&
      r->getSize() == phi2->getSize();
    TEST_EXIT(checkSize)("something went wrong\n");

    readDofVec(inStream, phi2, mesh);

    for (int i = r->getSize() - 1; i >= 0; --i) {
      (*r)[i] = Phi2ToR((*r)[i], eps);
      (*phi1)[i] = 1 - (*phi2)[i];
      (*levelSet)[i] = LevelSet((*r)[i]);
    }
  }
  

  std::string ProblemImplicit::getDofFilename(std::string path, int implMesh)
  {
    FUNCNAME("ProblemImplicit::getDofFilename()");
    std::string suffix = "[" + 
      boost::lexical_cast< std::string >(implMesh) +
      "]";

    std::string dofFilename("");
    Parameters::get(path + "dof file" + suffix, dofFilename);
    TEST_EXIT(dofFilename.length() == 0)("dof file was changed to \"arh file\" and reads an arh file now");
    Parameters::get(path + "arh file" + suffix, dofFilename);
    if (implMesh == 0 && dofFilename.length() == 0) {
      Parameters::get(path + "dof file", dofFilename);
      TEST_EXIT(dofFilename.length() == 0)("dof file was changed to \"arh file\" and reads an arh file now");
      Parameters::get(path + "arh file", dofFilename);
    }
    return dofFilename;    
  }


  double ProblemImplicit::getEpsilon(std::string path, int implMesh)
  {
    std::string suffix = "[" + 
      boost::lexical_cast< std::string >(implMesh) +
      "]";

    double eps(-1);
    Parameters::get(path + "eps" + suffix, eps);
    if (implMesh == 0 && eps < 0)
      Parameters::get(path + "eps", eps);
    return eps;
  }


  int ProblemImplicit::getType(std::string path, int implMesh)
  {
    std::string suffix = "[" + 
      boost::lexical_cast< std::string >(implMesh) +
      "]";
    int serType(-1);
    Parameters::get(path + "type" + suffix, serType);
    if (implMesh == 0 && serType < 0)
      Parameters::get(path + "type", serType);
    return serType;
  }


  DOFVector<double>* ProblemImplicit::getSignedDistance(unsigned int im , unsigned int m) 
  { 
    if (m >= r.size() || im >= r[m].size())
      return NULL;
    return (r[m])[im]; 
  }


  DOFVector<double>* ProblemImplicit::getPhi1(unsigned int im, unsigned int m)
  {
    if (m >= phi1.size() || im >= phi1[m].size())
      return NULL;

    return (phi1[m])[im];
  }


  DOFVector<double>* ProblemImplicit::getPhi2(unsigned int im, unsigned int m)
  {
    if (m >= phi2.size() || im >= phi2[m].size())
      return NULL;

    return (phi2[m])[im];
  }


  DOFVector<double>* ProblemImplicit::getLevelset(unsigned int im, unsigned int m)
  {
    if (m >= levelSet.size() || im >= levelSet[m].size())
      return NULL;

    return (levelSet[m])[im];
  }


  bool ProblemImplicit::createImplicitMesh()
  {
    //check each mesh if it's an implicit one
    r.resize(meshes.size());
    phi1.resize(meshes.size());
    phi2.resize(meshes.size());
    levelSet.resize(meshes.size());
    for (unsigned int i = 0; i < meshes.size(); ++i)
      createImplicitMesh(i);
    return true;
  }


  bool ProblemImplicit::createImplicitMesh(int p) 
  {
    std::string path = name + "->implicit mesh[" 
      + boost::lexical_cast< std::string >(p) + "]->";
    int nImplMeshes(0);
    Parameters::get(path + "nr meshes", nImplMeshes);
    if (nImplMeshes == 0)
      return false;
    r[p].resize(nImplMeshes, NULL);
    phi1[p].resize(nImplMeshes, NULL);
    phi2[p].resize(nImplMeshes, NULL);
    levelSet[p].resize(nImplMeshes, NULL);

    for ( int i = 0; i < nImplMeshes ; ++i ) {
      (r[p])[i] = new DOFVector< double >(getFeSpace(p), "r");
      (phi1[p])[i] = new DOFVector< double >(getFeSpace(p), "phi1");
      (phi2[p])[i] = new DOFVector< double >(getFeSpace(p), "phi2");
      (levelSet[p])[i] = new DOFVector< double >(getFeSpace(p), "levelSet");
      createImplicitMesh(path, i, p);
    }
    return true;
  }
  

  bool ProblemImplicit::createImplicitMesh(std::string path, int implMesh, 
					   int comp)
  {
    FUNCNAME("ProblemImplicit::createImplicitMesh()");

    std::string dofFilename = getDofFilename(path, implMesh);
    if (dofFilename.length() == 0)
      return false;

    double eps = getEpsilon(path, implMesh);
    if (eps < 0.0)
      return false;
    int serType = getType(path, implMesh);
    if (serType < 0)
      return false;

    TEST_EXIT(meshes[comp] != NULL)("the mesh was not created\n");
    
    switch (serType) {
    case 0:
      readR(dofFilename, eps, meshes[comp], implMesh, comp);
      break;
    case 1:
      readPhi1(dofFilename, eps, meshes[comp], implMesh, comp);
      break;
    case 2:
      readPhi2(dofFilename, eps, meshes[comp], implMesh, comp);
      break;
    default:
      WARNING("unknown implicit mesh type\n");
    }
    return true;
  }


  void ProblemImplicit::createMesh()
  {
    ProblemStat::createMesh();
#if 0
    implMesh.resize(meshes.size());
    for (unsigned int i = 0 ; i < meshes.size() ; ++i ) {
      std::string path = name + "->implicit mesh[" +
	boost::lexical_cast< std::string >(i) + "]";
      std::string meshFilename("");
      Parameters::get(path + "->macro file name", meshFilename);
      std::string serFilename("");
      Parameters::get(path + "->serialization file name", serFilename);
      implMesh[i] = true;
      if (meshFilename.length() > 0)
	meshes[i]->setName(path);
      else if (serFilename.length() > 0) {
	std::ifstream inStream(serFilename.c_str());
	meshes[i]->deserialize(inStream);
	inStream.close();
      } else
	implMesh[i] = false;
    }
#endif    
  }


  void ProblemImplicit::initialize(Flag initFlag, 
				   ProblemStatSeq* adaptProblem, 
				   Flag adoptFlag)
  {
    ProblemStat::initialize(initFlag);
    if ( initFlag.isSet(INIT_IMPLICIT_MESH) )
      createImplicitMesh();
  }

  ProblemImplicit::~ProblemImplicit()
  {
    for ( unsigned int p(0); p < meshes.size(); ++p ) {
      for ( unsigned int i = 0; i < r[p].size() ; ++i )
	  if ( r[p][i] != NULL)
	    delete r[p][i];
       for ( unsigned int i(0); i < phi1[p].size(); ++i )
	 if ( phi1[p][i] != NULL)
	    delete phi1[p][i];
       for ( unsigned int i(0); i < phi2[p].size(); ++i )
	 if ( phi2[p][i] != NULL)
	    delete phi2[p][i];
       for ( unsigned int i(0); i < levelSet[p].size(); ++i )
	 if ( levelSet[p][i] != NULL)
	    delete levelSet[p][i];
    }

  }
}
