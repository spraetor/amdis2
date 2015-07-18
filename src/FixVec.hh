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


namespace AMDiS {

  template<typename T>
  void DimVec<T>::multMatrixVec(DimMat<T> &m, DimVec<T> &v)
  {
    T *mIt, *thisIt;
    for (thisIt = this->begin(), mIt = m.begin();
	 thisIt != this->end();
	 thisIt++) {
      *thisIt = 0;
      for (T* vIt = v.begin(); vIt != v.end(); vIt++, mIt++)
	*thisIt += *vIt * *mIt;
    }
  }

  template<typename T>
  void WorldVector<T>::multMatrixVec(const WorldMatrix<T> &m, const WorldVector<T> &v)
  {
    FUNCNAME_DBG("WorldVector<T>::multMatrix()");
    TEST_EXIT_DBG(m.getNumRows() == this->getSize())("invalide size\n");
    TEST_EXIT_DBG(v.getSize() == this->getSize())("invalide size\n");

    T const* mIt;
    T* thisIt;
    for (thisIt = this->begin(), mIt = m.begin();
	 thisIt != this->end(); 
	 thisIt++) {
      *thisIt = 0;
      
      for (T const* vIt = v.begin(); vIt != v.end(); vIt++, mIt++) 
	*thisIt += *vIt * *mIt;      
    }
  }

  template<typename T>
  bool WorldMatrix<T>::isDiagMatrix() const
  {
    for (int i = 0; i < this->getSize(); i++)
      for (int j = i + 1; j < this->getSize(); j++)
	if (abs((*this)[i][j]) > DBL_TOL || abs((*this)[j][i]) > DBL_TOL) 
	  return(false);

    return(true);
  }

  template<typename T>
  bool WorldMatrix<T>::isSymmetric() const
  {
    for (int i = 0; i < this->getSize(); i++)
      for (int j = i + 1; j < this->getSize(); j++)
	if (abs((*this)[i][j] - (*this)[j][i]) > DBL_TOL) 
	  return false;

    return true;
  }
 
  template<typename T>
  void WorldMatrix<T>::setDiag(T value)
  {
    for (int i = 0; i < this->rows; i++)
      this->valArray[i * this->cols + i] = value;
  }

  template<typename T>
  void WorldMatrix<T>::vecProduct(const WorldVector<T>& v1, const WorldVector<T>& v2)
  {
    FUNCNAME_DBG("WorldMatrix<T>::vecProduct()");

    TEST_EXIT_DBG(v1.getSize() == v2.getSize())("size(v1) != size(v2), %d != %d\n", v1.getSize(), v2.getSize());
    TEST_EXIT_DBG(v1.getSize() == this->getNumRows())("size(v1) != num_rows(this), %d != %d\n", v1.getSize(), this->getNumRows());

    T *thisIt = this->begin();

    for (T const* v1It = v1.begin(); v1It != v1.end(); v1It++)
      for (T const* v2It = v2.begin(); v2It != v2.end(); v2It++, thisIt++)
	*thisIt = *v1It * *v2It;
  }
  
}
