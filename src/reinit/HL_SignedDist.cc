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


#include "HL_SignedDist.h"

namespace reinit
{

  using namespace AMDiS;

  const Flag HL_SignedDist::VEL_EXT                = 0X01L;
  const Flag HL_SignedDist::VEL_EXT_FROM_VEL_FIELD = 0X02L;

  void HL_SignedDist::calcSignedDistFct(AdaptInfo* adaptInfo_,
                                        const DOFVector<double>* lS_DOF_,
                                        DOFVector<double>* sD_DOF_,
                                        ElementFunction<double>* elFct)
  {
    adaptInfo = adaptInfo_;

    if (elFct == NULL)
    {
      TEST_EXIT(lS_DOF_)("illegal level set function !\n");
      TEST_EXIT(lS_DOF_->getFeSpace() == sD_DOF_->getFeSpace())("DOF vectors do not have the same FE space!\n");
    }
    TEST_EXIT(sD_DOF_)("illegal DOF vector for signed distance function !\n");
    TEST_EXIT(sD_DOF_->getFeSpace()->getBasisFcts()->getDegree() == 1)("does it work for non-linear finite elements ?\n");

    lS_DOF = lS_DOF_;
    sD_DOF = sD_DOF_;
    feSpace = sD_DOF->getFeSpace();

    // ===== Initialization. =====
    initialize(elFct);

    // ===== Print level set function. =====
    if (lS_DOF)
      printLevelSetFct();

    // ===== Create DOF vector to mark boundary vertices. =====
    if (bound_DOF)
      delete bound_DOF;
    bound_DOF = new DOFVector<double>(feSpace, "bound_DOF");
    bound_DOF->set(0.0);

    // ===== Boundary vertex and boundary value initialization. =====
    initializeBoundary();

    // ===== Print boundary initialization. =====
    printBoundInitFct();

    // ===== Calculate distance function. =====
    HL_updateIteration();

    // ===== Transformation to signed distance function. =====
    if (elFct == NULL)
      setSign();

    // ===== Print calculated signed distance function. =====
    printSignedDistFct();
  }


  void HL_SignedDist::calcSignedDistFct(AdaptInfo* adaptInfo_, DOFVector<double>* lS_DOF_)
  {
    TEST_EXIT(lS_DOF_)("illegal level set function lS_DOF_ !\n");

    sD_DOF = new DOFVector<double>(lS_DOF_->getFeSpace(), "sD_DOF");

    calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, NULL);

    *lS_DOF_ = *sD_DOF;

    delete sD_DOF;
  }


  void HL_SignedDist::calcVelocityExt(AdaptInfo* adaptInfo_,
                                      DOFVector<double>* origVel_DOF_,
                                      DOFVector<double>* vel_DOF_,
                                      const DOFVector<double>* lS_DOF_,
                                      DOFVector<double>* sD_DOF_,
                                      bool calcSDFct,
                                      ElementFunction<double>* elFct)
  {
    TEST_EXIT(velExtType.isSet(VEL_EXT))("illegal velocity extension type !\n");
    TEST_EXIT(velExt)("velExt not defined !\n");

    velExt->setVelocity(origVel_DOF_, vel_DOF_);

    velExt->printOrigVelDOF(adaptInfo_);

    if (calcSDFct || sD_DOF_ != NULL)
    {
      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF_, elFct);
    }
    else
    {

      TEST_EXIT(vel_DOF_)("illegal level set function lS_DOF_ !\n");

      sD_DOF = new DOFVector<double>(vel_DOF_->getFeSpace(), "sD_DOF");

      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, elFct);

      delete sD_DOF;
    }

    velExt->printVelDOF(adaptInfo_);
  }


  void HL_SignedDist::calcVelocityExt(AdaptInfo* adaptInfo_,
                                      DOFVector<double>* vel_DOF_,
                                      DOFVector<double>* lS_DOF_,
                                      bool calcSDFct,
                                      ElementFunction<double>* elFct)
  {
    TEST_EXIT(velExtType.isSet(VEL_EXT))("illegal velocity extension type !\n");
    TEST_EXIT(velExt)("velExt not defined !\n");
    TEST_EXIT(vel_DOF_)("velocity vector vel_DOF_ not defined !\n");

    if (calcSDFct)
    {
      TEST_EXIT(lS_DOF_)("illegal level set function lS_DOF_ !\n");
    }

    DOFVector<double>* newVel_DOF_ =
      new DOFVector<double>(vel_DOF_->getFeSpace(), "vel_DOF_");

    velExt->setVelocity(vel_DOF_, newVel_DOF_);

    velExt->printOrigVelDOF(adaptInfo_);

    if (calcSDFct && elFct == NULL)
    {
      calcSignedDistFct(adaptInfo_, lS_DOF_);
    }
    else
    {
      sD_DOF = new DOFVector<double>(vel_DOF_->getFeSpace(), "sD_DOF");

      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, elFct);

      if (calcSDFct)
        lS_DOF_->copy(*sD_DOF);

      delete sD_DOF;
    }

    vel_DOF_->copy(*newVel_DOF_);

    velExt->printVelDOF(adaptInfo_);

    delete newVel_DOF_;
  }


  void HL_SignedDist::calcVelocityExt(AdaptInfo* adaptInfo_,
                                      std::vector<DOFVector<double> *> origVel_DOF_,
                                      std::vector<DOFVector<double> *> vel_DOF_,
                                      const DOFVector<double>* lS_DOF_,
                                      DOFVector<double>* sD_DOF_,
                                      bool calcSDFct,
                                      ElementFunction<double>* elFct)
  {
    TEST_EXIT(velExtType.isSet(VEL_EXT))("illegal velocity extension type !\n");
    TEST_EXIT(velExt)("velExt not defined !\n");

    velExt->setVelocity(origVel_DOF_, vel_DOF_);

    velExt->printOrigVelDOF(adaptInfo_);

    if (calcSDFct || sD_DOF_ != NULL)
    {
      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF_, elFct);
    }
    else
    {

      TEST_EXIT(vel_DOF_[0])("illegal level set function lS_DOF_ !\n");

      sD_DOF = new DOFVector<double>(vel_DOF_[0]->getFeSpace(), "sD_DOF");

      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, elFct);

      delete sD_DOF;
    }

    velExt->printVelDOF(adaptInfo_);
  }


  void HL_SignedDist::calcVelocityExt(AdaptInfo* adaptInfo_,
                                      std::vector<DOFVector<double> *> vel_DOF_,
                                      DOFVector<double>* lS_DOF_,
                                      bool calcSDFct,
                                      ElementFunction<double>* elFct)
  {
    TEST_EXIT(velExtType.isSet(VEL_EXT))("illegal velocity extension type !\n");
    TEST_EXIT(velExt)("velExt not defined !\n");

    int nVelDOFs = vel_DOF_.size();

    if (calcSDFct)
    {
      TEST_EXIT(lS_DOF_)("illegal level set function lS_DOF_ !\n");
    }

    std::vector<DOFVector<double> *> newVel_DOF_(nVelDOFs);
    for (int nV = 0; nV < nVelDOFs; ++nV)
    {
      TEST_EXIT(vel_DOF_[nV])("velocity vector vel_DOF_ not defined !\n");

      newVel_DOF_[nV] = new DOFVector<double>((vel_DOF_[nV])->getFeSpace(), "vel_DOF_");
    }

    velExt->setVelocity(vel_DOF_, newVel_DOF_);

    velExt->printOrigVelDOF(adaptInfo_);

    if (calcSDFct && elFct == NULL)
    {
      calcSignedDistFct(adaptInfo_, lS_DOF_);
    }
    else
    {
      sD_DOF = new DOFVector<double>(vel_DOF_[0]->getFeSpace(), "sD_DOF");

      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, elFct);

      if (calcSDFct)
        lS_DOF_->copy(*sD_DOF);

      delete sD_DOF;
    }

    for (int nV=0; nV<nVelDOFs; ++nV)
      (vel_DOF_[nV])->copy(*(newVel_DOF_[nV]));

    velExt->printVelDOF(adaptInfo_);

    for (int nV = 0; nV < nVelDOFs; ++nV)
      delete newVel_DOF_[nV];
  }


  void HL_SignedDist::calcVelocityExtFromVelocityField(AdaptInfo* adaptInfo_,
      std::vector<DOFVector<double> *>& velField_,
      DOFVector<double>* vel_DOF_,
      const DOFVector<double>* lS_DOF_,
      DOFVector<double>* sD_DOF_,
      bool calcSDFct,
      ElementFunction<double>* elFct)
  {
    TEST_EXIT(velExtType.isSet(VEL_EXT_FROM_VEL_FIELD))
    ("illegal velocity extension type !\n");
    TEST_EXIT(velExt)("velExt not defined !\n");
    TEST_EXIT(elFct == NULL)("not implemented yet for elFct != NULL !\n");

    ((VelocityExtFromVelocityField*)(velExt))->setVelocityField(velField_, lS_DOF_, vel_DOF_);

    if (calcSDFct)
    {
      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF_, elFct);
    }
    else
    {

      TEST_EXIT(vel_DOF_)("illegal velocity vector vel_DOF_ !\n");

      sD_DOF = new DOFVector<double>(vel_DOF_->getFeSpace(), "sD_DOF");

      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, elFct);

      delete sD_DOF;
    }

    velExt->printVelDOF(adaptInfo_);
  }


  void HL_SignedDist::calcVelocityExtFromVelocityField(AdaptInfo* adaptInfo_,
      std::vector<DOFVector<double> *>& velField_,
      DOFVector<double>* vel_DOF_,
      DOFVector<double>* lS_DOF_,
      bool calcSDFct,
      ElementFunction<double>* elFct)
  {
    TEST_EXIT(velExtType.isSet(VEL_EXT_FROM_VEL_FIELD))
    ("illegal velocity extension type !\n");
    TEST_EXIT(velExt)("velExt not defined !\n");
    TEST_EXIT(elFct == NULL)("not implemented yet for elFct != NULL !\n");

    if (calcSDFct)
    {
      TEST_EXIT(lS_DOF_)("illegal level set function lS_DOF_ !\n");
    }

    ((VelocityExtFromVelocityField*)(velExt))->setVelocityField(velField_, lS_DOF_, vel_DOF_);

    if (calcSDFct && elFct == NULL)
    {
      calcSignedDistFct(adaptInfo_, lS_DOF_);
    }
    else
    {

      TEST_EXIT(vel_DOF_)("illegal velocity vector vel_DOF_ !\n");

      sD_DOF = new DOFVector<double>(vel_DOF_->getFeSpace(), "sD_DOF");

      calcSignedDistFct(adaptInfo_, lS_DOF_, sD_DOF, elFct);

      if (calcSDFct)
        lS_DOF_->copy(*sD_DOF);

      delete sD_DOF;
    }

    velExt->printVelDOF(adaptInfo_);
  }


  void HL_SignedDist::initialize(ElementFunction<double>* elFct)
  {
    if (elUpdate)
      delete elUpdate;
    if (bndElDist)
      delete bndElDist;

    if (elLS)
      delete elLS;
    if (phi)
      delete phi;

    // ===== Create ElementLevelSet. =====
    if (elFct == NULL)
    {
      phi = new ElementFunctionDOFVec<double>(lS_DOF);
      elLS = new ElementLevelSet("ElementLevelSet", phi, feSpace->getMesh());
    }
    else
    {
      elLS = new ElementLevelSet("ElementLevelSet", elFct, feSpace->getMesh());
    }


    // ===== Define boundary initialization strategy. =====
    int boundInitFlag;
    Parameters::get(name + "->boundary initialization", boundInitFlag);
    switch (boundInitFlag)
    {
    case 0:
      bndElDist = new BoundaryElementLevelSetDist(elLS, dim);
      break;
    case 1:
      bndElDist = new BoundaryElementTopDist(elLS, dim, velExt);
      break;
    case 2:
      bndElDist = new BoundaryElementEdgeDist(elLS, dim);
      break;
    case 3:
      bndElDist = new BoundaryElementNormalDist(elLS, dim);
      break;
    default:
      ERROR_EXIT("illegal boundary initialization !\n");
      break;
    }

    // ===== Create elUpdate. =====
    TEST_EXIT(velExt == NULL || boundInitFlag == 1)
    ("velocity extension only works with topological boundary element initialization method !\n");

    switch (dim)
    {
    case 2:
      elUpdate = new ElementUpdate_2d(velExt);
      break;
    case 3:
      elUpdate = new ElementUpdate_3d(velExt);
      break;
    default:
      ERROR_EXIT("illegal dimension !\n");
      break;
    }
  }


  void HL_SignedDist::setSign()
  {
    DOFVector<double>::Iterator it_sD(const_cast<DOFVector<double> *>(sD_DOF), USED_DOFS);
    DOFVector<double>::Iterator it_lS(const_cast<DOFVector<double> *>(lS_DOF), USED_DOFS);

    // Vertex lies inside the zero level set ?
    for (it_sD.reset(), it_lS.reset(); !it_sD.end(); ++it_sD, ++it_lS)
      if ((*it_lS) < 0)
        (*it_sD) *= -1;
  }

}
