#ifdef HAVE_PNG

#include "PngWriter.hpp"

#include <float.h>
#include <png.h>

#include "DataCollector.hpp"
#include "DOFVector.hpp"
#include "Traverse.hpp"

namespace AMDiS
{
  namespace io
  {

    int PngWriter::writeFile(std::string filename, int imageType)
    {
      FUNCNAME("PngWriter::writeFile()");

      double minX = DBL_MAX, minY = DBL_MAX;
      double maxX = DBL_MIN, maxY = DBL_MIN;

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(dataCollector->getMesh(), -1,
                                           Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
      double pointdist = std::min(distance(elInfo->getCoord(0), elInfo->getCoord(1)),
                                  std::min(distance(elInfo->getCoord(1), elInfo->getCoord(2)),
                                           distance(elInfo->getCoord(2), elInfo->getCoord(0))));
      while (elInfo)
      {
        for (int i = 0; i < 3; i++)
        {
          double x = (elInfo->getCoord(i))[0];
          double y = (elInfo->getCoord(i))[1];
          if (x < minX) minX = x;
          if (x > maxX) maxX = x;
          if (y < minY) minY = y;
          if (y > maxY) maxY = y;
        }

        elInfo = stack.traverseNext(elInfo);
      }

      TEST_EXIT(minX == 0.0 && minY == 0.0)("Only supported for minX = minY = 0.0!\n");
      TEST_EXIT(pointdist > 0.0)("This should not happen!\n");

      int imageX = static_cast<int>(maxX / pointdist) + 1;
      int imageY = static_cast<int>(maxY / pointdist) + 1;

      png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                            NULL, NULL, NULL);
      if (!png_ptr)
        return 0;

      png_bytep* rowPointers = new png_bytep[imageY];
      for (int i = 0; i < imageY; i++)
      {
        //       rowPointers[i] = (png_byte*)png_malloc(png_ptr,
        // 					     (imageType == 0 ? imageX : imageX * 3));
        rowPointers[i] = (png_byte*)malloc(sizeof(png_byte) * imageX * 3);
      }

      const BasisFunction* basisFcts = dataCollector->getFeSpace()->getBasisFcts();
      std::vector<DegreeOfFreedom> localDofs(3);
      DOFVector<double>* dofvalues = dataCollector->getValues();

      elInfo = stack.traverseFirst(dataCollector->getMesh(), -1,
                                   Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
      while (elInfo)
      {
        basisFcts->getLocalIndices(elInfo->getElement(),
                                   dataCollector->getFeSpace()->getAdmin(),
                                   localDofs);

        for (int i = 0; i < 3; i++)
        {
          if (imageType == 0)
          {
            int indexX = static_cast<int>((elInfo->getCoord(i))[0] / pointdist);
            int indexY = static_cast<int>((elInfo->getCoord(i))[1] / pointdist);
            rowPointers[indexY][indexX] =
              static_cast<unsigned char>((*dofvalues)[localDofs[i]]);
          }
          else
          {
            int indexX = static_cast<int>((elInfo->getCoord(i))[0] / pointdist);
            int indexY = static_cast<int>((elInfo->getCoord(i))[1] / pointdist);

            TEST_EXIT(indexX >= 0 && indexX < imageX)("X-index out of range!");
            TEST_EXIT(indexY >= 0 && indexY < imageY)("Y-index out of range!");

            int value = static_cast<int>((*dofvalues)[localDofs[i]]);

            unsigned char r = value % 256;
            unsigned char g = (value - r % (256 * 256)) / 256;
            unsigned char b = (value - r - g) / (256 * 256);
            rowPointers[indexY][indexX * 3] = r;
            rowPointers[indexY][indexX * 3 + 1] = g;
            rowPointers[indexY][indexX * 3 + 2] = b;
          }

        }

        elInfo = stack.traverseNext(elInfo);
      }

      FILE* fp = fopen(filename.c_str(), "wb");
      TEST_EXIT(fp)("Cannot open file for writing matrix picture file!\n");

      png_infop info_ptr = png_create_info_struct(png_ptr);

      if (!info_ptr)
      {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        return 0;
      }

      if (setjmp(png_jmpbuf(png_ptr)))
      {
        return 0;
      }

      png_init_io(png_ptr, fp);

      png_set_IHDR(png_ptr, info_ptr, imageX, imageY, 8,
                   (imageType == 0 ? PNG_COLOR_TYPE_GRAY : PNG_COLOR_TYPE_RGB),
                   PNG_INTERLACE_NONE,
                   PNG_COMPRESSION_TYPE_DEFAULT,
                   PNG_FILTER_TYPE_DEFAULT);

      png_set_rows(png_ptr, info_ptr, rowPointers);

      png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, png_voidp_NULL);

      png_destroy_write_struct(&png_ptr, &info_ptr);

      fclose(fp);

      delete [] rowPointers;
      return 0;
    }

  }
} // end namespace io, AMDiS

#endif
