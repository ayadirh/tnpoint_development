/**
 * @brief GIN index for temporal points.
 */

#include "pg_point/tpoint_gin.h"

/* C */
#include <assert.h>
#include <float.h>
/* PostgreSQL */
#include <postgres.h>
#include <access/gist.h>
#include <utils/float.h>
#include <utils/timestamp.h>
/* MEOS */
#include <meos.h>
#include <meos_internal.h>
#include "general/temporal_util.h"
/* MobilityDB */
#include "pg_general/temporal.h"
#include "pg_general/temporal_catalog.h"

/*****************************************************************************
 * GIN compare methods
 *****************************************************************************/
PG_FUNCTION_INFO_V1(tnpoint_gin_compare);
/**
 * GIN compare method for temporal points
 */
PGDLLEXPORT int
tnpoint_gin_compare(PG_FUNCTION_ARGS)
{
    Datum a  = PG_GETARG_POINTER(0);
    Datum b  = PG_GETARG_POINTER(1);
    int32      result = NULL;

  PG_RETURN_INT32(result);
}

/*****************************************************************************
 * GIN extract value methods
 *****************************************************************************/
PG_FUNCTION_INFO_V1(tnpoint_gin_extract_value);
/**
 * GIN extract value method for temporal points
 */
PGDLLEXPORT Datum
tnpoint_gin_extract_value(PG_FUNCTION_ARGS)
{
    // Datum itemValue  = PG_GETARG_POINTER(0);
    // Datum itemValue  = PG_GETARG_DATUM(0);
    int32      *nkeys = (int32 *) PG_GETARG_POINTER(1);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(2);
    Datum      *result = NULL;

    // itemValue = PointerGetDatum(PG_DETOAST_DATUM(itemValue));
    Temporal *temp = PG_GETARG_TEMPORAL_P(0);
    int count;
    int64 *routes = tnpoint_routes(temp, &count);
    *nkeys = count;

    // for(int i =0; i<count; i++){
    //     result[i]=routes[i];
    // }
    // ArrayType *result_array = int64arr_array(routes, count);

    // pfree(routes);
    PG_FREE_IF_COPY(temp, 0);

  // PG_RETURN_POINTER(result);
    PG_RETURN_POINTER(routes);
}

// #define ginOverlapsStrategyNumber 1
// #define ginContainsStrategyNumber 2
// #define ginIsContainedByStrategyNumber 3
// #define ginEqualStrategyNumber 4
// #define     GinOverlapStrategy   1 
// #define     GinContainsStrategy   2 
// #define     GinContainedStrategy   3 
// #define     GinEqualStrategy   4
/*****************************************************************************
 * GIN extract query methods
 *****************************************************************************/
PG_FUNCTION_INFO_V1(tnpoint_gin_extract_query);
/**
 * GIN extract value method for temporal points
 */
PGDLLEXPORT Datum
tnpoint_gin_extract_query(PG_FUNCTION_ARGS)
{
    Datum query  = PG_GETARG_POINTER(0);
    int32      *nentries = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Datum extra_data ** = PG_GETARG_POINTER(4);
    bool **nullFlags = PG_GETARG_POINTER(5);
    int32      *searchMode = (int32 *) PG_GETARG_POINTER(6);
    Datum      *result = NULL;

    // Temporal *temp = PG_GETARG_TEMPORAL_P(0);
    // int count;
    // int64 *routes = tnpoint_routes(temp, &count);
    // *nkeys = count;

    switch(strategy){
        case GinOverlapStrategy :
            //need to match the sequence of routeid
            elog(ERROR, "to work on it");
            result = NULL;
            break;
        case GinContainsStrategy:
            //if type = tnpoint
            //tnpoint_routes();
            break;
        case GinContainedStrategy:
            result = NULL;
            break;
        case GinEqualStrategy:
            elog(ERROR, "to work on it");
            result = NULL;
            break;
        default:
          elog(ERROR, "unrecognized strategy number: %d", strategy);
          result = NULL;    /* keep compiler quiet */
          break;
    }

  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * GIN consistent methods
 *****************************************************************************/
PG_FUNCTION_INFO_V1(tnpoint_gin_consistent);
/**
 * GIN consistent method for temporal points
 */
PGDLLEXPORT Datum
tnpoint_gin_consistent(PG_FUNCTION_ARGS)
{
    bool *check = (bool *) PG_GETARG_POINTER(0);
    StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(1);
    Datum query = PG_GETARG_POINTER(2);
    int32 nkeys = (int32) PG_GETARG_POINTER(3);
    Datum *extra_data = PG_GETARG_POINTER(4);
    bool *recheck = (bool *) PG_GETARG_POINTER(5);
    Datum *queryKeys = PG_GETARG_POINTER(6);
    bool *nullFlags = (bool *) PG_GETARG_POINTER(7);

    Datum      *result = NULL;

  PG_RETURN_BOOL(result);
}