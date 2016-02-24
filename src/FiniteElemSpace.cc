#include "FiniteElemSpace.h"
#include "DOFAdmin.h"
#include "DOFVector.h"
#include "BasisFunction.h"
#include "Mesh.h"
#include <algorithm>

// using namespace std;

namespace AMDiS
{
  std::vector<FiniteElemSpace*> FiniteElemSpace::feSpaces;

  FiniteElemSpace::FiniteElemSpace(DOFAdmin* admin_,
                                   const BasisFunction* bas_fcts,
                                   Mesh* aMesh,
                                   std::string aString)
    : name(aString),
      admin(admin_),
      basFcts(bas_fcts),
      mesh(aMesh)
  {
    FUNCNAME("FiniteElemSpace::FiniteElemSpace()");

    TEST_EXIT(mesh)("No mesh!\n");
    TEST_EXIT(basFcts)("No basis functions!\n");

    if (!admin)
    {
      const DOFAdmin* admin_local = NULL;
      const DimVec<int>* ndof = NULL;

      ndof = basFcts->getNumberOfDofs();
      TEST_EXIT(ndof)("no n_dof or basFcts->n_dof\n");

      for (int i = 0; i < mesh->getNumberOfDOFAdmin(); i++)
      {
        admin_local = &(mesh->getDofAdmin(i));
        int j = 0;
        for (; j <= mesh->getDim(); j++)
          if (admin_local->getNumberOfDofs(j) != (*ndof)[j])
            break;
        if (j > mesh->getDim())
          break;
        admin_local = NULL;
      }

      if (!admin_local)
        admin_local = mesh->createDOFAdmin(name, *ndof);

      admin = const_cast<DOFAdmin*>(admin_local);
    }

    feSpaces.push_back(this);
  }


  FiniteElemSpace& FiniteElemSpace::operator=(const FiniteElemSpace& feSpace)
  {
    if (&feSpace == this)
      return *this;

    mesh = new Mesh(feSpace.mesh->getName(), feSpace.mesh->getDim());
    *mesh = *(feSpace.mesh);
    admin = &(const_cast<DOFAdmin&>(mesh->getDofAdmin(0)));
    basFcts = feSpace.basFcts;

    // TODO: remove nonsense error-messages
    TEST_EXIT(feSpace.admin == &(feSpace.mesh->getDofAdmin(0)))
    ("Gut, das muss ich mir noch mal ueberlegen!\n");

    return *this;
  }


  FiniteElemSpace* FiniteElemSpace::provideFeSpace(DOFAdmin* admin,
      const BasisFunction* basFcts,
      Mesh* mesh,
      std::string name_)
  {
    for (size_t i = 0; i < feSpaces.size(); i++)
      if (feSpaces[i]->basFcts == basFcts &&
          feSpaces[i]->mesh == mesh &&
          (!admin || (admin && feSpaces[i]->admin == admin)))
        return feSpaces[i];

    return new FiniteElemSpace(admin, basFcts, mesh, name_);
  }


  void FiniteElemSpace::destroyFeSpaces()
  {
    for (size_t i = 0; i < feSpaces.size(); i++)
      delete feSpaces[i];

    feSpaces.resize(0);
  }

#if DEBUG
  FiniteElemSpace* FiniteElemSpace::provideFeSpace(Mesh* mesh)
  {
    FUNCNAME("FiniteElemSpace::provideFeSpace()");

    for (size_t i = 0; i < feSpaces.size(); i++)
      if (feSpaces[i]->mesh == mesh)
        return feSpaces[i];

    ERROR_EXIT("FE space not found!\n");

    return NULL;
  }
#endif


  int FiniteElemSpace::calcMemoryUsage()
  {
    int result = sizeof(FiniteElemSpace);
    result += mesh->calcMemoryUsage();
    return result;
  }


  void FiniteElemSpace::clear()
  {
    for (size_t i = 0; i < feSpaces.size(); i++)
    {
      if (feSpaces[i])
      {
        delete feSpaces[i];
        feSpaces[i] = NULL;
      }
    }
  }


  const FiniteElemSpace*
  FiniteElemSpace::getHighest(std::vector<const FiniteElemSpace*>& feSpaces)
  {
    const FiniteElemSpace* feSpace = feSpaces[0];
    for (size_t i = 1; i < feSpaces.size(); i++)
      if (feSpaces[i]->getBasisFcts()->getDegree() >
          feSpace->getBasisFcts()->getDegree())
        feSpace = feSpaces[i];

    return feSpace;
  }

} // end namespace AMDiS
