#pragma once

// std c++ headers
#include <string>

// Boost includes
#include <boost/property_tree/detail/file_parser_error.hpp>

namespace boost { namespace property_tree { namespace amdis_parser
{

    //! Json parser error
    class amdis_parser_error: public file_parser_error
    {
    public:
        amdis_parser_error(const std::string &message, 
                          const std::string &filename, 
                          unsigned long line): 
            file_parser_error(message, filename, line)
        { 
        }
    };

} } }
