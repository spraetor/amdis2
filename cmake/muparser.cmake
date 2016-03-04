# specify the target and requirements for the muparser library

set(MUPARSER_SOURCE_DIR ${BASE_DIR}/lib/muparser_v2_2_4/src)
set(MUPARSER_INCLUDE_DIR ${BASE_DIR}/lib/muparser_v2_2_4/include)
add_library(muparser
    ${MUPARSER_SOURCE_DIR}/muParser.cpp
    ${MUPARSER_SOURCE_DIR}/muParserBase.cpp
    ${MUPARSER_SOURCE_DIR}/muParserBytecode.cpp
    ${MUPARSER_SOURCE_DIR}/muParserCallback.cpp
    ${MUPARSER_SOURCE_DIR}/muParserDLL.cpp
    ${MUPARSER_SOURCE_DIR}/muParserError.cpp
    ${MUPARSER_SOURCE_DIR}/muParserInt.cpp
    ${MUPARSER_SOURCE_DIR}/muParserTest.cpp
    ${MUPARSER_SOURCE_DIR}/muParserTokenReader.cpp)

target_include_directories(muparser PUBLIC ${MUPARSER_INCLUDE_DIR})

set_property(TARGET muparser PROPERTY CXX_STANDARD 11)
set_property(TARGET muparser PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(muparser PRIVATE HAS_CXX11=1)

# specify how to install this target:
# -----------------------------------

file(GLOB MUPARSER_HEADERS "${MUPARSER_SOURCE_DIR}/include/*.h")
install(FILES ${MUPARSER_HEADERS}
	DESTINATION include/amdis/muparser)
install(TARGETS muparser
	DESTINATION lib/amdis/ )