#pragma once

// std c++ headers
#include <fstream>
#include <locale>
#include <string>

// Boost includes
#include <boost/property_tree/ptree.hpp>

// AMDiS includes
#include "utility/amdis_parser_read.hpp"
// #include "utility/amdis_parser_write.hpp"
#include "utility/amdis_parser_error.hpp"


namespace boost { namespace property_tree { namespace amdis_parser
{

    /**
     * Read JSON from a the given stream and translate it to a property tree.
     * @note Clears existing contents of property tree.  In case of error the
     *       property tree unmodified.
     * @note Items of JSON arrays are translated into ptree keys with empty
     *       names. Members of objects are translated into named keys.
     * @note JSON data can be a string, a numeric value, or one of literals
     *       "null", "true" and "false". During parse, any of the above is
     *       copied verbatim into ptree data string.
     * @throw amdis_parser_error In case of error deserializing the property
     *                          tree.
     * @param stream Stream from which to read in the property tree.
     * @param[out] pt The property tree to populate.
     */
    template<class Ptree>
    void read_amdis(std::basic_istream<
                       typename Ptree::key_type::value_type
                   > &stream,
                   Ptree &pt)
    {
        read_amdis_internal(stream, pt, std::string());
    }

    /**
     * Read JSON from a the given file and translate it to a property tree.
     * @note Clears existing contents of property tree.  In case of error the
     *       property tree unmodified.
     * @note Items of JSON arrays are translated into ptree keys with empty
     *       names. Members of objects are translated into named keys.
     * @note JSON data can be a string, a numeric value, or one of literals
     *       "null", "true" and "false". During parse, any of the above is
     *       copied verbatim into ptree data string.
     * @throw amdis_parser_error In case of error deserializing the property
     *                          tree.
     * @param filename Name of file from which to read in the property tree.
     * @param[out] pt The property tree to populate.
     * @param loc The locale to use when reading in the file contents.
     */

    template<class Ptree>
    void read_amdis(const std::string &filename,
                   Ptree &pt,
                   const std::locale &loc = std::locale())
    {
        std::basic_ifstream<typename Ptree::key_type::value_type>
            stream(filename.c_str());
        if (!stream)
            BOOST_PROPERTY_TREE_THROW(amdis_parser_error(
                "cannot open file", filename, 0));
        stream.imbue(loc);
        read_amdis_internal(stream, pt, filename);
    }

    /**
     * Translates the property tree to JSON and writes it the given output
     * stream.
     * @note Any property tree key containing only unnamed subkeys will be
     *       rendered as JSON arrays.
     * @pre @e pt cannot contain keys that have both subkeys and non-empty data.
     * @throw amdis_parser_error In case of error translating the property tree
     *                          to JSON or writing to the output stream.
     * @param stream The stream to which to write the JSON representation of the
     *               property tree.
     * @param pt The property tree to tranlsate to JSON and output.
     * @param pretty Whether to pretty-print. Defaults to true for backward
     *               compatibility.
     */
#if 0
    template<class Ptree>
    void write_amdis(std::basic_ostream<
                        typename Ptree::key_type::value_type
                    > &stream,
                    const Ptree &pt,
                    bool pretty = true)
    {
        write_amdis_internal(stream, pt, std::string(), pretty);
    }
#endif
    /**
     * Translates the property tree to JSON and writes it the given file.
     * @note Any property tree key containing only unnamed subkeys will be
     *       rendered as JSON arrays.
     * @pre @e pt cannot contain keys that have both subkeys and non-empty data.
     * @throw amdis_parser_error In case of error translating the property tree
     *                          to JSON or writing to the file.
     * @param filename The name of the file to which to write the JSON
     *                 representation of the property tree.
     * @param pt The property tree to translate to JSON and output.
     * @param loc The locale to use when writing out to the output file.
     * @param pretty Whether to pretty-print. Defaults to true and last place
     *               for backward compatibility.
     */
    template<class Ptree>
    void write_amdis(const std::string &filename,
                    const Ptree &pt,
                    const std::locale &loc = std::locale(),
                    bool pretty = true)
    {
        std::basic_ofstream<typename Ptree::key_type::value_type>
            stream(filename.c_str());
        if (!stream)
            BOOST_PROPERTY_TREE_THROW(amdis_parser_error(
                "cannot open file", filename, 0));
        stream.imbue(loc);
        write_amdis_internal(stream, pt, filename, pretty);
    }

} } }

namespace boost { namespace property_tree
{
    using amdis_parser::read_amdis;
//     using amdis_parser::write_amdis;
    using amdis_parser::amdis_parser_error;
} }
