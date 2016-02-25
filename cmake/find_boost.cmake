if (BOOST_ROOT)
    file(TO_CMAKE_PATH ${BOOST_ROOT} BOOST_ROOT)
endif (BOOST_ROOT)
if (BOOST_LIBRARYDIR)
    file(TO_CMAKE_PATH ${BOOST_LIBRARYDIR} BOOST_LIBRARYDIR)
endif (BOOST_LIBRARYDIR)

set(BOOST_VERSION "1.48")
set(BOOST_LIBS_REQUIRED system iostreams filesystem program_options date_time)
if (WIN32)
    list(APPEND BOOST_LIBS_REQUIRED zlib)
    if (ENABLE_COMPRESSION OR AMDIS_NEED_COMPRESSION)
	   list(APPEND BOOST_LIBS_REQUIRED bzip2)
    endif (ENABLE_COMPRESSION OR AMDIS_NEED_COMPRESSION)
endif (WIN32)

if (NOT BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
endif (NOT BUILD_SHARED_LIBS)
set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost ${BOOST_VERSION} REQUIRED ${BOOST_LIBS_REQUIRED})
if (Boost_FOUND)
    add_library(boost INTERFACE)
    target_include_directories(boost INTERFACE ${Boost_INCLUDE_DIR})
    target_link_libraries(boost INTERFACE  ${Boost_LIBRARIES})
#     target_link_libraries(boost INTERFACE Boost::boost Boost::system Boost::iostreams Boost::filesystem Boost::program_options Boost::date_time)
    target_link_libraries(amdis_base INTERFACE boost)

    if (MSVC_SHARED_LIBS)
        link_directories(${Boost_LIBRARY_DIRS})
        target_compile_definitions(amdis_base INTERFACE ${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
    endif (MSVC_SHARED_LIBS)
endif (Boost_FOUND)
