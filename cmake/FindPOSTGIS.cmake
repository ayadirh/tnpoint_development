# (c) 2015  pgRouting colaborators
#
# Finds the most recent postGIS for a particular postgreSQL
# We need this for the tests
#
# Usage:
# find_package(POSTGIS)
#
# The following variables are set if PostGIS is found:
# POSTGIS_FOUND             - Set to true when PostGIS is found.
# POSTGIS_LIBRARY           - if we ever need to link it
# POSTGIS_CONTROL           - if we ever need to see the contents
# POSTGIS_VERSION           - The full numbers
# POSTGIS_VERSION_STR       - The th PostGIS prefix

if (POSTGIS_FOUND)
  return()
endif()

if (NOT POSTGRESQL_FOUND)
  find_package(POSTGRESQL REQUIRED)
endif()

# Find PostGIS library
file(GLOB POSTGIS_LIBRARY "${POSTGRESQL_DYNLIB_DIR}/postgis-*.so")
if(POSTGIS_LIBRARY STREQUAL "")
  message(FATAL_ERROR "No PostGIS library have been found")
endif()
list(LENGTH POSTGIS_LIBRARY NO_POSTGIS_LIBRARIES)
if(NO_POSTGIS_LIBRARIES GREATER 1)
  message(FATAL_ERROR "Several versions of the PostGIS library have been found")
endif()

find_file(POSTGIS_CONTROL postgis.control
  PATHS "${POSTGRESQL_SHARE_DIR}/extension")

if (POSTGIS_CONTROL)
  file(READ ${POSTGIS_CONTROL} control_contents)
  string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" POSTGIS_VERSION ${control_contents})
  set(POSTGIS_VERSION_STR "PostGIS ${POSTGIS_VERSION}")
  string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\1" POSTGIS_MAJOR_VERSION ${POSTGIS_VERSION})
  string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\2" POSTGIS_MINOR_VERSION ${POSTGIS_VERSION})
  string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" "\\3" POSTGIS_MICRO_VERSION ${POSTGIS_VERSION})
  math(EXPR POSTGIS_VERSION_NUMBER "${POSTGIS_MAJOR_VERSION} * 10000 + ${POSTGIS_MINOR_VERSION} * 100 + ${POSTGIS_MICRO_VERSION}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(POSTGIS
  FOUND_VAR POSTGIS_FOUND
  REQUIRED_VARS POSTGIS_LIBRARY POSTGIS_CONTROL
  VERSION_VAR POSTGIS_VERSION
  FAIL_MESSAGE "Could NOT find PostGIS")

if (POSTGIS_FOUND)
  mark_as_advanced(POSTGIS_LIBRARY POSTGIS_CONTROL POSTGIS_VERSION POSTGIS_VERSION_STR)
endif()

message(STATUS "POSTGIS_LIBRARY: ${POSTGIS_LIBRARY}")
message(STATUS "POSTGIS_CONTROL: ${POSTGIS_CONTROL}")
message(STATUS "POSTGIS_VERSION: ${POSTGIS_VERSION}")
message(STATUS "POSTGIS_VERSION_STR: ${POSTGIS_VERSION_STR}")
