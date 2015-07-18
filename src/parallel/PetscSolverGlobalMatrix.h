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



/** \file PetscSolverGlobalMatrix.h */

#ifndef AMDIS_PETSC_SOLVER_GLOBAL_MATRIX_H
#define AMDIS_PETSC_SOLVER_GLOBAL_MATRIX_H

#include <boost/tuple/tuple.hpp>
#include "AMDiS_fwd.h"
#include "parallel/PetscSolver.h"

namespace AMDiS 
{
  namespace Parallel
  {

    /** \ingroup Solver
     * 
     * \brief
     * PETSc solver which creates a globally distributed matrix. Supports also
     * coarse space assembling.
     */
    class PetscSolverGlobalMatrix : public PetscSolver
    {
    public:
      /// Creator class
      class Creator : public LinearSolverCreator
      {
      public:
	virtual ~Creator() {}

	/// Returns a new PetscSolver object.
	LinearSolverInterface* create() 
	{ 
	  return new PetscSolverGlobalMatrix(this->name); 
	}
      };
      
      PetscSolverGlobalMatrix(std::string name, bool setOptions = true);

      void fillPetscMatrix(Matrix<DOFMatrix*> *mat);

      void fillPetscMatrixWithCoarseSpace(Matrix<DOFMatrix*> *mat);

      void fillPetscRhs(SystemVector *vec);

      void solvePetscMatrix(SystemVector &vec, AdaptInfo *adaptInfo);    

      void solveGlobal(Vec &rhs, Vec &sol);
      
      void extractVectorComponent(Vec input, int i, Vec *output, int numberOfComponents=1);
      
      void extractMatrixComponent(Mat input, int startRow, int numberOfRows, int startCol, int numberOfCols, Mat *output);

      void destroyMatrixData();

      void destroyVectorData();

    protected:
      void removeDirichletRows(Matrix<DOFMatrix*> *seqMat);

      void removeDirichletRows(SystemVector *seqVec);

      /// Reads field split information and creats a splitting based on 
      /// component numbers.
      void createFieldSplit(PC pc);

      /** \brief 
      * Creates a new field split for a preconditioner object.
      *
      * \param[in] pc          PETSc preconditioner object, must be of
      *                        type PCFIELDSPLIT
      * \param[in] splitName   Name of the field split, can be used to set other
      *                        parameters of the PCFIELDSPLIT preconditioner on
      *                        the command line.
      * \param[in] components  System component numbers of the field split. At
      *                        the moment only continuous splits are allowed.
      */
      void createFieldSplit(PC pc, const char* splitName, std::vector<int> &components);

      /// Wrapper to create field split from only one component.
      void createFieldSplit(PC pc, const char* splitName, int component)
      {
	std::vector<int> components;
	components.push_back(component);
	createFieldSplit(pc, splitName, components);
      }
      
      

      virtual void initSolver(KSP &ksp);

      virtual void exitSolver(KSP &ksp);

      virtual void initPreconditioner(PC pc);
      virtual void initPreconditioner(const Matrix<DOFMatrix*>& A, const Mat& fullMatrix)
      {
	initPreconditioner(getPc());
      }

      virtual void exitPreconditioner(PC pc);

      /// Takes a DOF matrix and sends the values to the global PETSc matrix.
      void setDofMatrix(DOFMatrix* mat, int rowComp = 0, int colComp = 0);

      /// Takes a DOF vector and sends its values to a given PETSc vector.
      void setDofVector(Vec vecInterior,
			Vec vecCoarse,
			DOFVector<double>* vec, 
			int rowCompxo, 
			bool rankOnly = false);

      inline void setDofVector(Vec vecInterior,
			      DOFVector<double>* vec, 
			      int rowComp, 
			      bool rankOnly = false)
      {
	setDofVector(vecInterior, PETSC_NULL, vec, rowComp, rankOnly);
      }

      inline void setDofVector(Vec vecInterior, 
			      SystemVector &vec, 
			      bool rankOnly = false)
      {
	for (int i = 0; i < vec.getSize(); i++)
	  setDofVector(vecInterior, PETSC_NULL, vec.getDOFVector(i), i, rankOnly);
      }

      PetscSolver* createSubSolver(int component, std::string kspPrefix);

      void setConstantNullSpace(KSP ksp, int constFeSpace, bool test = false);

      void setConstantNullSpace(KSP ksp);

    protected:
      bool zeroStartVector;

      /// If true, after parallel assembling, information about the matrix
      /// are printed.
      bool printMatInfo;
    };

  } // end namespace Parallel
  
} // end namespace AMDiS

#endif

