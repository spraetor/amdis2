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




#ifndef HL_SIGNEDDISTBORNEMANN
#define HL_SIGNEDDISTBORNEMANN

#include <queue>
#include <time.h>
#include "ElInfo.h"
#include "FixVec.h"
#include "Traverse.h"
#include "ElementLevelSet.h"
#include "BoundaryElementDist.h"
#include "ElementUpdate.h"
#include "ElementUpdate_2d.h"
#include "ElementUpdate_3d.h"
#include "HL_SignedDist.h"

#include "SMIAdapter.h"
#include "smi.h"

namespace reinit {

using namespace AMDiS;

typedef struct
{
  int ElNum;
  int VertNum;
} Vert_Struct;

class HL_SignedDistBornemann : public HL_SignedDist
{
public:
  HL_SignedDistBornemann(const char *name_,int dim_)
    : HL_SignedDist(name_, dim_),
      smiAdapter(NULL)
  {
    FUNCNAME("HL_SignedDistBornemann::HL_SignedDistBornemann");

    // ===== Read parameters from init file. =====
    Parameters::get(name + "->tolerance", "%f", &tol);
    Parameters::get(name + "->count_how_often_saved_in_list", "%d", &count_in_list);
    Parameters::get(name + "->save_in_list->the ..th", "%d", &print_in_list);

    TEST_EXIT(tol > 0)("illegal tolerance !\n");
  }

 protected:
  /**
   * Initializes the boundary: calculation of the distance of boundary 
   * vertices to the interface. 
   * Interface is given by \ref lS_DOF and result is stored in 
   * \ref sD_DOF.
   */
  void initializeBoundary() 
  {
    FUNCNAME("HL_SignedDistBornemann::initializeBoundary()");

    TraverseStack stack;
    FixVec<double, VERTEX> distVec(dim);
    int elStatus;
    const int *elVertStatusVec;
    int numIntersecPoints;
    int NumVertIntPoints;
    int ElNum;
    double TIME;

    //=== transvering Mesh to SMI and add quantities
    Timer t;
    Mesh_to_SMI_and_quantity ();
    TIME = t.elapsed();

    cout <<"Zeit zum Transformieren nach SMI: "<<TIME<<" sec.\n";

    // ===== All non-boundary vertices are initialized with "infinity". =====
    sD_DOF->set(inftyValue);

    // ===== Traverse mesh and initialize boundary elements. =====
    const int nBasFcts = feSpace->getBasisFcts()->getNumber();
    DegreeOfFreedom *locInd = new DegreeOfFreedom[nBasFcts];

    ElInfo *elInfo = stack.traverseFirst(feSpace->getMesh(),
					 -1, 
					 Mesh::CALL_LEAF_EL | 
					 Mesh::FILL_BOUND |
					 Mesh::FILL_COORDS);
    while(elInfo) {

      // Get local indices of vertices.
      feSpace->getBasisFcts()->getLocalIndices(
	  const_cast<Element *>(elInfo->getElement()),
	  const_cast<DOFAdmin *>(feSpace->getAdmin()),
	  locInd); 

      // Get element status.
      elStatus = ElementLevelSet::createElementLevelSet(elInfo);

      //Get some information for creating the first list
      
      elVertStatusVec = ElementLevelSet::getElVertStatusVec();
      numIntersecPoints = ElementLevelSet::getNumElIntersecPoints();
      
      //Beginn creating the first list
      
      NumVertIntPoints = ElementLevelSet::getNumVertIntPoints();
      ElNum = elInfo->getElement()->getIndex();
      
      createLists( elStatus, elVertStatusVec, numIntersecPoints, ElNum, NumVertIntPoints, locInd);
            
      //end creating two lists
      
      // Is element cut by the interface ?
      if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

	// Reset element distance vector.
	for (int i=0; i<=dim; ++i) {
	  distVec[i] = inftyValue;
	}

	// Mark all vertices as boundary vertices.
	for (int i=0; i<=dim; ++i) {
	  (*bound_DOF)[locInd[i]] = 1.0;
	}



	// Calculate distance for all vertices.
	if (bndElDist->calcDistOnBoundaryElement(elInfo, distVec) !=
	    ElementLevelSet::LEVEL_SET_BOUNDARY) {
	  ERROR_EXIT("error in distance calculation !\n");
	}
	else {

	  // If distance is smaller, correct to new distance.
	  for (int i=0; i<=dim; ++i) {
	    if ((*sD_DOF)[locInd[i]] > distVec[i]) {
	      (*sD_DOF)[locInd[i]] = distVec[i];
	    }
	  }
	}
      }  // end of: elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY

      else if (ElementLevelSet::getNumVertIntPoints() != 0){
	// Interface cuts element only in vertices.
	elVertStatusVec = ElementLevelSet::getElVertStatusVec();

	for (int i=0; i<=dim; ++i) {
	  if (elVertStatusVec[i] == ElementLevelSet::LEVEL_SET_BOUNDARY) {
	    (*bound_DOF)[locInd[i]] = 1.0;
	    (*sD_DOF)[locInd[i]] = 0.0;
	  }
	}
      } 

      elInfo = stack.traverseNext(elInfo);
    }  // end of: mesh traverse

    FREE_MEMORY(locInd, DegreeOfFreedom, nBasFcts);

    return;
  };

  /**
   * Calculates the distance function and stores result in sD_DOF.
   * Requirement: The boundary values are already set in \ref sD_DOF.
   */
  void HL_updateIteration()
  {
    int numNodes;
    int *nodeIndices;
    int max_q3 = 0;
    int max_q4 = 0;
    int sum_q3 = 0;
    int sum_q4 = 0;
    double TIME;
    
    SMI_Get_all_nodes(1,1,const_cast<DegreeOfFreedom*>(&numNodes),&nodeIndices);
    int *values_q3 = new int [numNodes];
    int *values_q4 = new int [numNodes];
    
    Parameters::get("SignedDist->count updates", count_updates);
    
    Timer t;
    traverseListElement();
    TIME = t.elapsed();

    cout<<"Zeit zum durchlaufen der ersten Liste: "<<TIME<<" sec.\n\n";
    t.reset();

    traversingListELVert ( sD_DOF );
    TIME = t.elapsed();

    cout<<"Zeit zum Durchlaufen der zweiten Liste: "<<TIME<<" sec.\n\n";
    
    std::string smiOutFile;
    Parameters::get("SignedDist->count updates->output->filename", &smiOutFile);
    cout << "count updates Ausgabe-Datei:  " << smiOutFile.c_str() << "\n\n";
    
    ofstream out (smiOutFile.c_str());
    
    SMI_Get_quantity_values(1,1,3,SMI_TYPE_INT,1,numNodes,nodeIndices,values_q3);
    SMI_Get_quantity_values(1,1,4,SMI_TYPE_INT,1,numNodes,nodeIndices,values_q4);  

    for (int i=0; i<numNodes; i++)
      {
	out<<nodeIndices[i]<<" "<<values_q3[i]<<" "<<values_q4[i]<<"\n";

	if(max_q3 < values_q3[i])
	  {
	    max_q3 = values_q3[i];
	  }
	if(max_q4 < values_q4[i])
	  {
	    max_q4 = values_q4[i];
	  }
	sum_q3 = sum_q3 + values_q3[i];
	sum_q4 = sum_q4 + values_q4[i];
      }

    out<<"\n\n maximale Anzahl an versuchten Updates auf einem Knoten: "<<max_q3<<"\n maximale Anzahl an durchgefuehrten Updates auf einem Knoten:  "<<max_q4<<"\n";
    out<<"\n Summe aller versuchten Updates: "<<sum_q3<<"\n Summe aller durchgefuehrten updates: "<<sum_q4<<"\n";

    out<<"Anzahl an Knoten: "<<numNodes<<"\n";

    out.close();

 
  

    SMI_End_write_transaction(1,1);

    return;
  };

  /**
   *function for transvering Mesh to SMi an Adding the quantities
   **/

  void Mesh_to_SMI_and_quantity ()
    {
      if(!smiAdapter)
	{
	  smiAdapter = new SMIAdapter(1,1,const_cast<FiniteElemSpace*>(feSpace),1,-1);
	}
      
      // cout << "\n\n\tSMI-Adapter angelegt !\n\n";
      
      //====== transfer Mesh to SMI ======================
      
      smiAdapter->addNeighbourInfo();
      
      // cout << "\n\n\tNachbarschafts-Infos in SMI ergaenzt !\n\n";
      
      smiAdapter-> transferMeshToSMI();
      
      //   cout << "\n\n\tGitter nach SMI geschrieben !\n\n";
    
      SMI_Begin_write_transaction(1,1);
    
      int nul=0;
      //which pair of element and node is saved in list
      SMI_Add_quantity(1, 1, 1, SMI_LOC_ELEM, SMI_TYPE_INT,dim+1, &nul);
      //which node is a boundary-node
      SMI_Add_quantity(1,1,2,SMI_LOC_NODE, SMI_TYPE_INT,1,&nul);
      //saves the number of tried updates
      SMI_Add_quantity(1,1,3,SMI_LOC_NODE,SMI_TYPE_INT,1,&nul);
      //saves the number of updates
      SMI_Add_quantity(1,1,4,SMI_LOC_NODE,SMI_TYPE_INT,1,&nul);
      //how often a node is saved in the second list
      SMI_Add_quantity(1,1,5,SMI_LOC_NODE,SMI_TYPE_INT,1,&nul);
      
      // cout << "\n\n\tQuantities in SMI ergaenzt !\n\n";
    }

  //Begin fuction for creating two lists
  //*
    
    void createLists(const int elStatus, 
		     const int *elVertStatusVec, 
		     const int numIntersecPoints, 
		     const int ElNum, 
		     const int NumVertIntPoints,
		     const DegreeOfFreedom *locInd)
      {
	Vert_Struct vertStruct;
	int *sv = new int [dim+1];
	int value=0;
	int value_q5 = 0;
	int * nodeIndices;
	/*    int *nodes = new int [dim+1]; */

    if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY)
      {
	List_Element.push( ElNum );

	//saving which nodes are boundary-nodes
	for (int i=0; i<=dim; i++)
	  {
	    sv[i] = 1;
	  }
	SMI_Get_elems(1,1,1,const_cast<DegreeOfFreedom*>(&ElNum),NULL,&nodeIndices,NULL,NULL);
	SMI_Set_quantity_values(1,1,2,SMI_TYPE_INT,1,3,nodeIndices,sv);
      }
    //if the elemet isn't a boundary-element, but the interface cuts the FE in two nodes
    else if (NumVertIntPoints == dim)
      {
	for (int i=0; i<=dim; i++)
	  {
	    if (elVertStatusVec[i] == ElementLevelSet::LEVEL_SET_BOUNDARY && elVertStatusVec[(i+1)%3] ==ElementLevelSet::LEVEL_SET_BOUNDARY)
	      {
		value = 1;
		SMI_Get_elems(1,1,1,const_cast<DegreeOfFreedom*>(&ElNum),NULL,&nodeIndices,NULL,NULL);
		SMI_Set_quantity_values(1,1,2,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&nodeIndices[i]),&value);
		SMI_Set_quantity_values(1,1,2,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&nodeIndices[(i+1)%3]),&value);

		//saving which node of this element will be included into the second list 
		//part 1
		for (int j=0; j<=dim; j++)         //ACHTUNG: aus < ein <= gemacht
		  { 
		    sv[j] = 0; 
		  } 
		sv[(i+2)%3] = 1;

		//include pair into the second list
		vertStruct.ElNum = ElNum ; 
		vertStruct.VertNum = locInd[(i+2)%3]; 
 		List_El_Vert.push( vertStruct ); 

		//saving which node of this element will be included into the second list 
		//part 2
 		SMI_Set_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,const_cast<DegreeOfFreedom*>(&ElNum),sv); 

		//counts how often a node is saved in the second list
		if(count_in_list == 1)
		  {
		    SMI_Get_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct.VertNum),&value_q5);
		    value_q5 = value_q5 + 1;
		    SMI_Set_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct.VertNum),&value_q5);
		  }
	      }
	  }
      }
    else 
      {
	//saving that none of the nodes of this element is saved in the second list
	for (int i=0; i<=dim; i++)
	  {
	    sv[i] = 0;
	  }
	SMI_Set_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,const_cast<DegreeOfFreedom*>(&ElNum),sv);
      }
    delete [] sv;
      };
 
//End function for creating two lists
//
 
//function for traversing the list "List_Element"
//
 
 void traverseListElement()
   {
     Vert_Struct vertStruct;
     
     int Element;
     int *neighbour=new int[dim+1];
     int *oppVertices=new int[dim+1];
     int *values=new int[dim*(dim+1)];
     int *value = new int [dim+1];
     int value_q2;
     int value_q5;
     int *nodeIndices;
     
     while (List_Element.size() != 0)
       {
	 Element = List_Element.front();
	 smiAdapter->getNeighbourInfo(Element, neighbour, oppVertices);
	 
	  for ( int i=0; i<=dim; i++)
	    {
	      SMI_Get_elems(1,1,1,&neighbour[i],NULL,&nodeIndices, NULL,NULL);
	      SMI_Get_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,&neighbour[i],value);
	      SMI_Get_quantity_values(1,1,2,SMI_TYPE_INT,1,1,&nodeIndices[oppVertices[i]],&value_q2);
	      if (value[oppVertices[i]] == 0 && value_q2 == 0)
		{
		  vertStruct.ElNum = neighbour[i] ;
		  vertStruct.VertNum = nodeIndices[oppVertices[i]] ;                           
		  List_El_Vert.push( vertStruct );
		  value[oppVertices[i]] = 1;
		  
		  //counts how often a node is saved in the second list
		  if(count_in_list == 1)
		    {
		      SMI_Get_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct.VertNum),&value_q5);
		      value_q5 = value_q5 + 1;
		      SMI_Set_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct.VertNum),&value_q5);
		    }
		}
	      SMI_Set_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,&neighbour[i],value);  
	    }
	  List_Element.pop();
       }
     delete [] neighbour;
     delete [] oppVertices;
     delete [] values;
     delete [] value;
   };
 
 //end of the function for traversing the list "List_Element"
 //
 
 
 /*
  *gets the neighbour according to the node "VertNum_in" an two of its nodes
  *returns "false" if the called neighbour exists and "true" if not
  */
 
 bool getNextNeighbour (const int ElNum_in, const int Vert_Up_in, const int VertNum_in, int &ElNum_out, int &VertNum_1_out, int &VertNum_2_out)
   {
     int *neighbour = new int[dim+1];
     int *oppVertices = new int[dim+1];
     int *nodeIndicesOfElem;
     int locVertNum;
     
     smiAdapter->getNeighbourInfo(ElNum_in, neighbour,oppVertices);
     
     //which local node is the node "VertNum_in"?
      SMI_Get_elems(1,1,1,const_cast<DegreeOfFreedom*>(&(ElNum_in)),NULL,&nodeIndicesOfElem, NULL, NULL);
      for (int i=0; i<=dim; i++)
	{
	  if(nodeIndicesOfElem[i] == VertNum_in)
	    {
	      locVertNum = i;
	    }
	}
      
      ElNum_out = neighbour [locVertNum];
      
      
      
      SMI_Get_elems(1,1,1,const_cast<DegreeOfFreedom*>(&(ElNum_out)),NULL,&nodeIndicesOfElem, NULL, NULL);
      VertNum_1_out = nodeIndicesOfElem [oppVertices[locVertNum]];
      for (int i=0; i<=dim; i++)
	{
	  if (nodeIndicesOfElem[i] != Vert_Up_in && nodeIndicesOfElem[i] != VertNum_1_out)
	    {  
	      VertNum_2_out = nodeIndicesOfElem[i];
	    }
	}
      
      delete [] neighbour;
      delete [] oppVertices;
           
      if (ElNum_out == -1)
	{
	  return true;
	}
      
      return false;
      
   }
 
 
 /*
  * checking whetherthe element "ElNum_in" has to be included into the second list
  *if yes it will be included
  */
 
 void includeIntoList (int ElNum_in, int VertNum_1_in, int VertNum_2_in)
   {
     Vert_Struct vertStruct_0;
     Vert_Struct vertStruct_1;
     
     int *nodeIndices;
     int *valuesINT = new int[dim+1];
     int *locVert = new int[dim];
     int value_q2;
     int value_q5;
     
     if(ElNum_in != -1)
       {
	 
	  SMI_Get_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,const_cast<DegreeOfFreedom*>(&ElNum_in),valuesINT);
	  
	  //which locla node is the node VerNum_1_in, which one is the node VertNum_2_in?
	  SMI_Get_elems (1, 1, 1, &ElNum_in, NULL, &nodeIndices, NULL, NULL);
	  for (int i=0; i<=dim; i++)
	    {
	      if(nodeIndices[i] == VertNum_1_in)
		{
		  locVert[0] = i;
		}
	      if(nodeIndices[i] == VertNum_2_in)
		{
		  locVert[1] = i;
		}
	    }
	  
	  //if the pair of element and node isn't saved in the second List an if it isn't a boundary-node include into second list
	  SMI_Get_quantity_values(1,1,2,SMI_TYPE_INT,1,1,&VertNum_1_in, &value_q2);
	  if(valuesINT[locVert[0]] == 0 && value_q2 == 0)
	    {
	      vertStruct_0.ElNum = ElNum_in;
	      vertStruct_0.VertNum = VertNum_1_in;
	      List_El_Vert.push( vertStruct_0 );
	      
	      valuesINT[locVert[0]] = 1;
	      SMI_Set_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,&ElNum_in,&valuesINT);

	      //counts how often a node is saved in the second list
	      if(count_in_list == 1)
		{
		  SMI_Get_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct_0.VertNum),&value_q5);
		  value_q5 = value_q5 + 1;
		  SMI_Set_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct_0.VertNum),&value_q5);
		}
	    }
	  
	  //if the pair of element and node isn't saved in the second List an if it isn't a boundary-node include into second list
	  SMI_Get_quantity_values(1,1,2,SMI_TYPE_INT,1,1,&VertNum_2_in, &value_q2);
	  if(valuesINT[locVert[1]] == 0 && value_q2 == 0)
	    {
	      vertStruct_1.ElNum = ElNum_in;
	      vertStruct_1.VertNum = VertNum_2_in;
	      List_El_Vert.push( vertStruct_1 );
	      
	      valuesINT[locVert[1]] = 1;
	      SMI_Set_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,&ElNum_in,&valuesINT);

	      //counts how often a node is saved in the second list
	      if(count_in_list == 1)
		{
		  SMI_Get_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct_1.VertNum),&value_q5);
		  value_q5 = value_q5 + 1;
		  SMI_Set_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&vertStruct_1.VertNum),&value_q5);
		}
	    }
       }
     
     delete [] valuesINT;
     delete [] locVert;
   }
 

 /**
  *function needed in the function "traversingListELVert"
  **/
 int getNext_node_l_r (int elem_l_r_in, int node_l_r_in, const int Vert)
   {
     int *nodeIndicesOfElem;
  
     SMI_Get_elems(1,1,1,const_cast<DegreeOfFreedom*>(&elem_l_r_in),NULL,&nodeIndicesOfElem, NULL, NULL);
     for (int i=0; i<=dim; i++)
       {
	 if(nodeIndicesOfElem[i] != Vert && nodeIndicesOfElem[i] != node_l_r_in)
	   {
	     return nodeIndicesOfElem[i];
	   }
       }
   }
 
 //Beginn function for traversing the list "List_El_Vert"
 //
 
 void traversingListELVert ( DOFVector<double> *boundVal_DOF )
   {
     Vert_Struct vertStruct_0;
     Vert_Struct vertStruct_1;
     
     Vert_Struct El_Vert;
     
     int Vert;
     
    
    
     //int * new_neighbour_l = new int [dim+1];
     //int *new_neighbour_r = new int [dim+1];
     //  int *oppVertices_l = new int[dim+1];
     //int *oppVertices_r = new int[dim+1];
    
     int locVert;
     int value_q3;
     int value_q4;
     int value_q5;
     int counter = 0;
     int *nodeIndices;
     double *coords = new double [(dim+1)*dim];
     int *valuesINT=new int[dim+1];
     double update;
     int value_q2;
     int counter_list = 0;
     FixVec<WorldVector<double> *, VERTEX> coordsOfNodes(dim);
     for ( int i=0; i<=dim; i++)
       {
	 coordsOfNodes[i] = new WorldVector<double>;
       }
     FixVec<double,VERTEX> boundVal(dim);
     WorldVector<double> helpVec;
      
      

     while (List_El_Vert.size() != 0)
       {
	 El_Vert = List_El_Vert.front();
	 Vert = El_Vert.VertNum;
	 
	 SMI_Get_elems (1, 1, 1, &(El_Vert.ElNum), NULL, &(nodeIndices), NULL, NULL);
	 SMI_Get_nodes (1, 1, 3, dim, nodeIndices, coords);
	 
	 counter = 0;
	 for (int i=0; i<=dim; i++)
	   {
	     if (nodeIndices[i] == El_Vert.VertNum)
	       {
		 (*(coordsOfNodes[dim]))[0]=coords[i*2];
		 (*(coordsOfNodes[dim]))[1]=coords[i*2+1];
		 boundVal[dim] = (*boundVal_DOF)[nodeIndices[i]];
		 locVert = i;
	       }
	     else
	       {
		 (*(coordsOfNodes[counter]))[0]=coords[i*2];
		 (*(coordsOfNodes[counter]))[1]=coords[i*2+1];
		 boundVal[counter] = (*boundVal_DOF)[nodeIndices[i]];
		 counter = counter+1;
	       }
	   }
	 
	 update = elUpdate->calcElementUpdate( coordsOfNodes, boundVal );
	 
	  if(count_updates == 1)
	    {
	      SMI_Get_quantity_values(1,1,3,SMI_TYPE_INT,1,1,&nodeIndices[locVert],&value_q3); 
	      value_q3 = value_q3 + 1;
	      SMI_Set_quantity_values(1,1,3,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&nodeIndices[locVert]),&value_q3);
	    }
	  
	  if (update < (*boundVal_DOF)[El_Vert.VertNum] && ((*boundVal_DOF)[El_Vert.VertNum]-update) > tol)
	    {
	      (*boundVal_DOF)[El_Vert.VertNum] = update;
	      
	      if(count_updates == 1)
		{
		  SMI_Get_quantity_values(1,1,4,SMI_TYPE_INT,1,1,&nodeIndices[locVert],&value_q4); 
		  value_q4 = value_q4 + 1;
		  SMI_Set_quantity_values(1,1,4,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(&nodeIndices[locVert]),&value_q4);
		}
	      search_and_include_comb(El_Vert.ElNum, Vert, nodeIndices);
	    } 
 
	  SMI_Get_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,&El_Vert.ElNum, valuesINT);           
	  List_El_Vert.pop();
	  valuesINT[locVert] = 0;
	  SMI_Set_quantity_values(1,1,1,SMI_TYPE_INT,dim+1,1,&El_Vert.ElNum, valuesINT);

	  //counts how often a node is saved in the second list
	  if(count_in_list == 1)
	    {
	      SMI_Get_quantity_values(1,1,5,SMI_TYPE_INT,1,1,const_cast<DegreeOfFreedom*>(& El_Vert.VertNum),&value_q5);
	      value_q5 = value_q5 - 1;
	      SMI_Set_quantity_values(1,1,5,SMI_TYPE_INT,dim+1,1,const_cast<DegreeOfFreedom*>(& El_Vert.VertNum),&value_q5);
	    }

	  if(count_in_list == 1)
	    {
	      counter_list++;
	      if(counter_list % print_in_list == 0)
		{
		  print_quantity_5(counter_list);
		}
	    }	 
       }

// delete [] new_neighbour_l;
//   delete [] new_neighbour_r;
// delete [] oppVertices_l;
//     delete [] oppVertices_r;
     delete [] coords;
     delete [] valuesINT;
     for ( int i=0; i<Global::getGeo(VERTEX,dim); i++)
       {
	 delete coordsOfNodes[i];
       }
   };
 
 //End of the function for traversing the list "List_El_Vert"
 //

 void search_and_include_comb(int ElNum, int Vert, int *nodeIndices)
   {
     bool stop_l = false;
     bool stop_r = false;
     int neighbour_l;
     int neighbour_r;
     int Vert_1_l;
     int Vert_2_l;
     int Vert_1_r;
     int Vert_2_r;
     int locVert;
     int counter;
     int elem_l;
     int elem_r;
     int node_l;
     int node_r;

     //here the process of including elements/nodes into the second list
     stop_l = false;
     stop_r = false;	      
     
     elem_l = ElNum;
     if (elem_l == -1)
       {
	 stop_l = true;
       }
     elem_r = ElNum;
     if (elem_r == -1)
       {
	 stop_r = true;
       }
     
     counter = 0;
     for (int i=0; i<=dim; i++)
       {
	 if(nodeIndices[i] != Vert && counter == 1)
	   {
	     node_r = nodeIndices[i];
	     counter = counter + 1;
	   }
	 if(nodeIndices[i] != Vert && counter == 0)
	   {
	     node_l = nodeIndices[i];
	     counter = counter + 1;
	   }
	 
	 if(nodeIndices[i] == Vert)
	   {
	     locVert = i;
	   }
       }
     
     while (!stop_l || !stop_r)
       {
	 //get next neighbour
	 if (!stop_l)
	   {
	     stop_l = getNextNeighbour (elem_l, Vert, node_l, neighbour_l, Vert_1_l, Vert_2_l);
	   }
	 if(!stop_r)
	   {
	     stop_r = getNextNeighbour (elem_r, Vert, node_r, neighbour_r, Vert_1_r, Vert_2_r);
	   }
	 
	 if(neighbour_l == elem_r)
	   {
	     break;
	   }
	 
	 //indclude into the second list (only if possible)
	 if (!stop_l)
	   {
	     includeIntoList (neighbour_l, Vert_1_l, Vert_2_l);
	     if(neighbour_l == neighbour_r)
	       {
		 break;
	       }
	     //"elem_l", "node_l", "elem_r", "node_r" have to be set on the next elements
	     node_l = getNext_node_l_r(elem_l, node_l,Vert);
	     elem_l = neighbour_l;
	   }
	 if (!stop_r && neighbour_l != neighbour_r)
	   {
	     includeIntoList (neighbour_r, Vert_1_r, Vert_2_r);
	     //"elem_l", "node_l", "elem_r", "node_r" have to be set on the next elements
	     node_r = getNext_node_l_r(elem_r, node_r, Vert);
	     elem_r = neighbour_r;
	   }
       }
   }

 //==================================
 void print_quantity_5 (int cntr)
   {
     int numNodes;
     int*nodeIndices;
     double *coords = new double[dim];
     SMI_Get_all_nodes(1,1,const_cast<DegreeOfFreedom*>(&numNodes),&nodeIndices);
     int *values_q5 = new int [numNodes];

      std::string q5_OutFile_0;
      std::string q5_OutFile_1;
      std::string q5_OutFile_2;
      std::string q5_OutFile_3;
      std::string q5_OutFile_4;

    Parameters::get("SignedDist->count saving->output->filename_0", q5_OutFile_0);
    Parameters::get("SignedDist->count saving->output->filename_1", q5_OutFile_1);
    Parameters::get("SignedDist->count saving->output->filename_2", q5_OutFile_2);
    Parameters::get("SignedDist->count saving->output->filename_3", q5_OutFile_3);
    Parameters::get("SignedDist->count saving->output->filename_4", q5_OutFile_4);
    
    
    char cntrStr[20];
    sprintf(cntrStr, "%d", cntr);
    q5_OutFile_0 += cntrStr;
    q5_OutFile_1 += cntrStr;
    q5_OutFile_2 += cntrStr;
    q5_OutFile_3 += cntrStr;
    q5_OutFile_4 += cntrStr;

    cout << "count saving Ausgabe-Datei_0:  " << q5_OutFile_0.c_str() << "\n\n";
    cout << "count saving Ausgabe-Datei_1:  " << q5_OutFile_1.c_str() << "\n\n";
    cout << "count saving Ausgabe-Datei_2:  " << q5_OutFile_2.c_str() << "\n\n";
    cout << "count saving Ausgabe-Datei_3:  " << q5_OutFile_3.c_str() << "\n\n";
    cout << "count saving Ausgabe-Datei_4:  " << q5_OutFile_4.c_str() << "\n\n";

    ofstream out_0 (q5_OutFile_0.c_str());
    ofstream out_1 (q5_OutFile_1.c_str());
    ofstream out_2 (q5_OutFile_2.c_str());
    ofstream out_3 (q5_OutFile_3.c_str());
    ofstream out_4 (q5_OutFile_4.c_str());
    
    SMI_Get_quantity_values(1,1,5,SMI_TYPE_INT,1,numNodes,nodeIndices,values_q5);
   
    for (int i=0; i<numNodes; i++)
      {

	SMI_Get_nodes (1,1,1,dim,const_cast<DegreeOfFreedom*>(&nodeIndices[i]),coords);

	if(values_q5[i] == 0)
	  {
	    out_0<<coords[0]<<" "<<coords[1]<<"\n";
	  }
	if(values_q5[i] == 1)
	  {
	    out_1<<coords[0]<<" "<<coords[1]<<"\n";
	  }
	if(values_q5[i] == 2)
	  {
	    out_2<<coords[0]<<" "<<coords[1]<<"\n";
	  }
	if(values_q5[i] == 3)
	  {
	    out_3<<coords[0]<<" "<<coords[1]<<"\n";
	  }
	if(values_q5[i] >= 4)
	  {
	    out_4<<coords[0]<<" "<<coords[1]<<"\n";
	  }


      }

    out_0.close();
    out_1.close();
    out_2.close();
    out_3.close();
    out_4.close();

    delete [] coords;
    delete [] values_q5;
   }
 //====================================

 protected:
 /**
  * Tolerance for Hopf-Lax update iteration loop.
  */
 double tol;

 /**
  *is needed for transfering the mesh to SMI
  */
 SMIAdapter *smiAdapter;

 /**
  * in this list boundary-elements are saved;
  * is needed for creating the list "List_El_Vert"
  */
 queue<int> List_Element;
  
 /**
  *in this list structs filled with element-number and vertex-number are saved;
  *is needed for traversing the mesh efficiently
  */
 queue<Vert_Struct> List_El_Vert;

 /**
  *0 ->do not count updates
  *1 ->count updates
  **/
 int count_updates; 

 /**
  *1 -> it will be count how often a node is saved in the second list
  *0 ->it will not be count
  **/
 int count_in_list; 

 int print_in_list;
};

} // end namespace reinit

using reinit::HL_SignedDistBornemann;

#endif  // HL_SIGNEDDISTBORNEMANN
