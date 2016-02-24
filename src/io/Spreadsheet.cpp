#include <algorithm>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "Spreadsheet.h"
#include "Global.h"

namespace AMDiS
{
  namespace io
  {

    using namespace std;
    using namespace boost;

    void Spreadsheet::write(string filename)
    {
      int nRows = static_cast<int>(data.size());
      int maxDataLength = 0;

      for (int i = 0; i < nRows; i++)
        maxDataLength = std::max(maxDataLength, static_cast<int>(data[i].size()));

      for (int i = 0; i < nRows; i++)
        data[i].resize(maxDataLength, 0.0);

      ofstream file;
      file.open(filename.c_str());

      for (int i = 0; i < nRows; i++)
      {
        for (int j = 0; j < maxDataLength; j++)
        {
          file << data[i][j];
          if (j + 1 < maxDataLength)
            file << " ";
          else
            file << "\n";
        }
      }

      file.close();
    }


    void Spreadsheet::read(string filename)
    {
      data.clear();

      string line;
      vector<string> lineSplit;

      ifstream file;
      file.open(filename.c_str());

      while (!file.eof())
      {
        getline(file, line);
        split(lineSplit, line, is_any_of(" ,"));

        if (lineSplit.size() == 0)
          continue;

        if ((lineSplit.size() == 1) && (lineSplit[0] == ""))
          continue;

        if (lineSplit[0].find("#") == 0)
          continue;

        vector<double> lineData(lineSplit.size());
        for (unsigned int i = 0; i < lineSplit.size(); i++)
          lineData[i] = boost::lexical_cast<double>(lineSplit[i]);

        data.push_back(lineData);
      }

      file.close();
    }

  }
} // end namespace io, AMDiS
