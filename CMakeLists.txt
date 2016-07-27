cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#---------------------------
set( PackageName DDKalTest )
#---------------------------

project(${PackageName})

# project version
SET( ${PackageName}_VERSION_MAJOR 1 )
SET( ${PackageName}_VERSION_MINOR 0 )
SET( ${PackageName}_VERSION_PATCH 0 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${CMAKE_SOURCE_DIR}/cmake ) 
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

#------------- set the default installation directory to be the source directory

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET( CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH  
    "install prefix path  - overwrite with -D CMAKE_INSTALL_PREFIX = ..." 
    FORCE )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} - overwrite with -D CMAKE_INSTALL_PREFIX" )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

find_package( DD4hep REQUIRED COMPONENTS DDRec )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake ) 
include( DD4hep )

find_package( ROOT REQUIRED )
if ( ${ROOT_VERSION_MAJOR} GREATER 5 )
  set( ROOT_COMPONENT_LIBRARIES Geom )
else()
  set( ROOT_COMPONENT_LIBRARIES Geom Reflex)
endif()


if(DD4HEP_USE_XERCESC)
  find_package(XercesC)
endif()

include(DD4hep_XML_setup)
 
#fg: fixme: for now need to link against geant4 libraries:
find_package(Geant4 REQUIRED gdml ui_all vis_all)
INCLUDE(${Geant4_USE_FILE})   # this also takes care of geant 4 definitions and include dirs
LINK_LIBRARIES( ${Geant4_LIBRARIES} )



#-------------------------------------------------------------
# add additional packages here

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES streamlog REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( KalTest 1.5 REQUIRED )
FIND_PACKAGE( streamlog REQUIRED )
FIND_PACKAGE( aidaTT REQUIRED )
FIND_PACKAGE( GSL REQUIRED )

#NB: need ROOT before KalTest (needed on case sensitive systems and classes that differ only in capitalization)
FOREACH( pkg KalTest LCIO streamlog aidaTT GSL)

  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
  ENDIF()

ENDFOREACH()

add_definitions( -D AIDATT_USE_DD4HEP ) 
#-------------------------------------------------------------

include_directories( ${CMAKE_SOURCE_DIR}/include
  ${DD4hep_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIR}
 )

file(GLOB sources ./src/*cc )

add_library(${PackageName} SHARED ${sources})

target_link_libraries(${PackageName} ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}
  ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} 
  )


SET_TARGET_PROPERTIES( ${PackageName} PROPERTIES VERSION ${${PackageName}_VERSION} SOVERSION ${${PackageName}_SOVERSION} )



add_executable( printSurfaces ./examples/printSurfaces.cc )
target_link_libraries( printSurfaces ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} ${PackageName})


#---Rootmap generation--------------------------------------------------------------

#if(APPLE)
#  dd4hep_generate_rootmap_apple(${PackageName} )
#else()
#  dd4hep_generate_rootmap(${PackageName} )
#endif()


#---- configure run environment ---------------
configure_file( ${DD4hep_ROOT}/cmake/thisdd4hep_package.sh.in  this${PackageName}.sh @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/this${PackageName}.sh
  DESTINATION bin
  )


#--- install target-------------------------------------


FILE(GLOB hfiles "ILD/include/*.h")
INSTALL(FILES ${hfiles} 
  DESTINATION include/${PackageName} 
  )


install(TARGETS ${PackageName}  printSurfaces
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  )
# to do: add corresponding uninstall...
#-------------------------------------------------------

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()


# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( DDKalTestConfig.cmake DDKalTestConfigVersion.cmake DDKalTestLibDeps.cmake )
