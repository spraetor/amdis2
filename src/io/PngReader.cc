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

#if defined HAVE_PNG

#include "PngReader.h"
#include "png.h"

#include "detail/PngReader.h"

using namespace std;

namespace AMDiS { namespace io {
  
  namespace PngReader
  {

    /** \brief
    * Copies the values of a value file to a DOF vector.
    */
    void readFile(std::string filename, DOFVector<double> *vec)
    {
      FUNCNAME("PngReader::readFile()");

      TEST_EXIT(filename != "")("Filename not specified!\n");
      TEST_EXIT(vec)("no DOF vector specified\n");

      png_structp png_ptr;
      png_infop info_ptr;
      FILE *fp;
      unsigned int sig_read = 0;
      int row, col;
      int bytesPerPixel = 0;

      // Open files and create the png data structures.
      if ((fp = fopen(filename.c_str(), "rb")) == NULL) {
	TEST_EXIT(0)("ERROR: file can not >be opened\n");
      }

      png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
      if (png_ptr == NULL) {
	TEST_EXIT(0)("ERROR in png_create_read_struct\n");
      }

      info_ptr = png_create_info_struct(png_ptr);
      if (info_ptr == NULL) {
	TEST_EXIT(0)("ERROR in png_create_info_struct\n");
      }

      if (setjmp(png_jmpbuf(png_ptr))) {
	TEST_EXIT(0)("ERROR in png_jmpbuf\n");
      }

      png_init_io(png_ptr, fp);
      png_set_sig_bytes(png_ptr, sig_read);

      // Read the whole png at once to the pointer info_ptr.
      png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

      bytesPerPixel = info_ptr->rowbytes /  info_ptr->width;

      cout << "Read image: " << filename << endl;
      cout << "Size: " << info_ptr->width << " x " << info_ptr->height << " pixel" << endl;
      cout << "Bytes per pixel: " << bytesPerPixel << endl;

      double value = 0;

      const BasisFunction *basFcts = vec->getFeSpace()->getBasisFcts();
      int numBasFcts = basFcts->getNumber();
      std::vector<DegreeOfFreedom> localIndices(numBasFcts);

      double xMin = 0.0, xMax = 1.0, yMin = 0.0, yMax = 1.0;
      detail::getMeshDimension(vec->getFeSpace()->getMesh(), xMin, xMax, yMin, yMax);
      
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(vec->getFeSpace()->getMesh(), -1,
					  Mesh::CALL_LEAF_EL |
					  Mesh::FILL_COORDS);
      while (elInfo) {
	Element *el = elInfo->getElement();
	basFcts->getLocalIndices(el, vec->getFeSpace()->getAdmin(), localIndices);
	
	for (int i = 0; i < numBasFcts; i++) {
	  double lambdaX = ((elInfo->getCoords())[i][0] - xMin)/(xMax - xMin);
	  double lambdaY = ((elInfo->getCoords())[i][1] - yMin)/(yMax - yMin);
	  col = static_cast<int>(lambdaX*(info_ptr->width-1));
	  row = static_cast<int>((1.0-lambdaY)*(info_ptr->height-1));
	  switch (bytesPerPixel) {
	    case 1:
	      value = static_cast<double>(info_ptr->row_pointers[row][col]);
	      break;
	    case 3:
	      value =
		1.0/3.0/255.0 * (static_cast<double>(info_ptr->row_pointers[row][3*col]) +
		  static_cast<double>(info_ptr->row_pointers[row][col*3 + 1]) +
		  static_cast<double>(info_ptr->row_pointers[row][col*3 + 2]));
	      break;
	    default:
	      TEST_EXIT(false)("ERROR: bytesPerPixel=%d is unknown case!\n",bytesPerPixel);
	  }
	  (*vec)[localIndices[i]] = value;
	}
	elInfo = stack.traverseNext(elInfo);
      }

      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(fp);
    }
  
  } // end namespace PngReader
} } // end namespace io, AMDiS

#endif
