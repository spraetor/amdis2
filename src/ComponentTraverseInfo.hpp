#pragma once

// std c++ headers
#include <set>

// AMDiS includes
#include "Global.hpp"

namespace AMDiS
{
  //forward declaration
  class FiniteElemSpace;


  class SingleComponentInfo
  {
  public:
    SingleComponentInfo()
      : rowFeSpace(NULL),
        colFeSpace(NULL),
        status(0)
    {}

    void setFeSpace(const FiniteElemSpace* row, const FiniteElemSpace* col = NULL)
    {
      rowFeSpace = row;
      colFeSpace = col;
    }

    void setAuxFeSpaces(std::set<const FiniteElemSpace*> feSpaces)
    {
      auxFeSpaces = feSpaces;
    }

    void addAuxFeSpace(const FiniteElemSpace* fe)
    {
      auxFeSpaces.insert(fe);
    }

    bool hasFeSpace() const
    {
      return rowFeSpace != NULL;
    }

    void updateStatus();

    size_t getNumAuxFeSpaces() const
    {
      return auxFeSpaces.size();
    }

    const FiniteElemSpace* getRowFeSpace() const
    {
      return rowFeSpace;
    }

    const FiniteElemSpace* getColFeSpace() const
    {
      return colFeSpace;
    }

    const FiniteElemSpace* getAuxFeSpace() const
    {
      FUNCNAME_DBG("SingleComponentInfo::getAuxFeSpace()");

      TEST_EXIT_DBG(auxFeSpaces.size() <= 1)("More than one aux FE space!\n");

      if (auxFeSpaces.size() == 1)
        return (*(auxFeSpaces.begin()));

      return NULL;
    }

    int getStatus() const
    {
      return status;
    }

  protected:
    const FiniteElemSpace* rowFeSpace;

    const FiniteElemSpace* colFeSpace;

    std::set<const FiniteElemSpace*> auxFeSpaces;

    /// Status of the component, see the possible status flags below.
    int status;

  public:
    /// Single component status flag: empty component, no fe spaces
    static constexpr int EMPTY = 0;

    /// Single component status flag: row = col, no aux
    static constexpr int EQ_SPACES_NO_AUX = 1;

    /// Single component status flag: row = col = aux
    static constexpr int EQ_SPACES_WITH_AUX = 2;

    /// Single component status flag: row = col, different aux
    static constexpr int EQ_SPACES_WITH_DIF_AUX = 3;

    /// Single component status flag: row, col, no aux
    static constexpr int DIF_SPACES_NO_AUX = 4;

    /// Single component status flag: row, col, aux either equal to row or to col
    static constexpr int DIF_SPACES_WITH_AUX = 5;

    /// Single component status flag: row, col, aux (at least 3 different fe spaces)
    static constexpr int DIF_SPACES_WITH_DIF_AUX = 6;
  };


  class ComponentTraverseInfo
  {
  public:
    ComponentTraverseInfo(int n)
      : nComponents(n)
    {
      resize(n);
    }

    void resize(int n)
    {
      nComponents = n;

      matrixComponents.resize(n);
      vectorComponents.resize(n);

      for (int i = 0; i < n; i++)
        matrixComponents[i].resize(n);
    }

    void updateStatus()
    {
      for (int i = 0; i < nComponents; i++)
      {
        for (int j = 0; j < nComponents; j++)
          matrixComponents[i][j].updateStatus();

        vectorComponents[i].updateStatus();
      }
    }

    SingleComponentInfo& getMatrix(int row, int col)
    {
      return matrixComponents[row][col];
    }

    SingleComponentInfo& getVector(int row)
    {
      return vectorComponents[row];
    }

    int getStatus(int row, int col) const
    {
      return matrixComponents[row][col].getStatus();
    }

    int getStatus(int row) const
    {
      return vectorComponents[row].getStatus();
    }

    /// Returns true, if for the given matrix component the row and the col FE spaces
    /// are equal. Note that this is also the case, if there is another aux FE space,
    /// that may be different from both, the row and the col FE spaces.
    bool eqSpaces(int row, int col) const
    {
      int status = matrixComponents[row][col].getStatus();

      return (status == SingleComponentInfo::EQ_SPACES_NO_AUX ||
              status == SingleComponentInfo::EQ_SPACES_WITH_AUX ||
              status == SingleComponentInfo::EQ_SPACES_WITH_DIF_AUX);
    }

    const FiniteElemSpace* getAuxFeSpace(int row, int col) const
    {
      return matrixComponents[row][col].getAuxFeSpace();
    }

    const FiniteElemSpace* getAuxFeSpace(int row) const
    {
      return vectorComponents[row].getAuxFeSpace();
    }

    /// Returns true if there is an aux FE that is different from row and col FE spaces.
    bool difAuxSpace(int row, int col) const
    {
      int status = matrixComponents[row][col].getStatus();

      return (status == SingleComponentInfo::EQ_SPACES_WITH_DIF_AUX ||
              status == SingleComponentInfo::DIF_SPACES_WITH_DIF_AUX);
    }

    /// Returns true if there is an aux FE that is different from row and col FE spaces.
    bool difAuxSpace(int row) const
    {
      int status = vectorComponents[row].getStatus();

      return (status == SingleComponentInfo::EQ_SPACES_WITH_DIF_AUX ||
              status == SingleComponentInfo::DIF_SPACES_WITH_DIF_AUX);
    }

    /** \brief
     * Returns the row FE space for a given row number, i.e., the FE space
     * of the diagonal matrix.
     *
     * \param[in]  row   Row number of the matrix line for which the FE space
     *                   should be returned.
     */
    const FiniteElemSpace* getRowFeSpace(int row) const;


    /** \brief
     * Returns the non row FE space for a given row number. This is either the
     * col FE space of an off diagonal matrix or the aux fe space of another
     * matrix in the row or of the right hand side vector. If there is no non row
     * FE space, this function returns a null pointer.
     *
     * \param[in]  row   Row number of the matrix line for which the non FE space
     *                   should be returned.
     */
    const FiniteElemSpace* getNonRowFeSpace(int row) const;

  protected:
    int nComponents;

    std::vector<std::vector<SingleComponentInfo>> matrixComponents;

    std::vector<SingleComponentInfo> vectorComponents;
  };

} // end namespace AMDiS
