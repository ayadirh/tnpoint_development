/*****************************************************************************
 *
 * This MobilityDB code is provided under The PostgreSQL License.
 * Copyright (c) 2016-2022, Université libre de Bruxelles and MobilityDB
 * contributors
 *
 * MobilityDB includes portions of PostGIS version 3 source code released
 * under the GNU General Public License (GPLv2 or later).
 * Copyright (c) 2001-2022, PostGIS contributors
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without a written
 * agreement is hereby granted, provided that the above copyright notice and
 * this paragraph and the following two paragraphs appear in all copies.
 *
 * IN NO EVENT SHALL UNIVERSITE LIBRE DE BRUXELLES BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
 * EVEN IF UNIVERSITE LIBRE DE BRUXELLES HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 * UNIVERSITE LIBRE DE BRUXELLES SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON
 * AN "AS IS" BASIS, AND UNIVERSITE LIBRE DE BRUXELLES HAS NO OBLIGATIONS TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
 *
 *****************************************************************************/

/**
 * @brief Basic functions for temporal network points.
 */

#include "npoint/tnpoint.h"

/* MEOS */
#include <meos.h>
#include "general/temporal_parser.h"
#include "general/temporal_util.h"
#include "general/lifting.h"
#include "npoint/tnpoint_static.h"
#include "npoint/tnpoint_parser.h"
/* MobilityDB */
#include "pg_general/temporal_catalog.h"
#include "pg_npoint/tnpoint_static.h"


/*****************************************************************************
 * General functions
 *****************************************************************************/

/**
 * @brief Convert a C array of int64 values into a PostgreSQL array
 */
ArrayType *
int64arr_array(const int64 *int64arr, int count)
{
  return construct_array((Datum *)int64arr, count, INT8OID, 8, true, 'd');
}

#if 0 /* not used */
/**
 * @brief Convert a C array of network point values into a PostgreSQL array
 */
ArrayType *
npointarr_array(Npoint **npointarr, int count)
{
  return construct_array((Datum *)npointarr, count, type_oid(T_NPOINT),
    sizeof(Npoint), false, 'd');
}
#endif

/**
 * @brief Convert a C array of network segment values into a PostgreSQL array
 */
ArrayType *
nsegmentarr_array(Nsegment **nsegmentarr, int count)
{
  return construct_array((Datum *)nsegmentarr, count, type_oid(T_NSEGMENT),
    sizeof(Nsegment), false, 'd');
}

/*****************************************************************************
 * Input/output functions
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tnpoint_in);
/**
 * @ingroup mobilitydb_temporal_in_out
 * @brief Input function for temporal network points
 * @sqlfunc tnpoint_in()
 */
PGDLLEXPORT Datum
Tnpoint_in(PG_FUNCTION_ARGS)
{
  char *input = PG_GETARG_CSTRING(0);
  Oid temptypid = PG_GETARG_OID(1);
  Temporal *result = temporal_parse(&input, oid_type(temptypid));
  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * Cast functions
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tnpoint_to_tgeompoint);
/**
 * @ingroup mobilitydb_temporal_cast
 * @brief Cast a temporal network point as a temporal geometric point
 * @sqlfunc tgeompoint()
 */
PGDLLEXPORT Datum
Tnpoint_to_tgeompoint(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Temporal *result = tnpoint_tgeompoint(temp);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Tgeompoint_to_tnpoint);
/**
 * @ingroup mobilitydb_temporal_cast
 * @brief Cast a temporal geometric point as a temporal network point
 * @sqlfunc tnpoint()
 */
PGDLLEXPORT Datum
Tgeompoint_to_tnpoint(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Temporal *result = tgeompoint_tnpoint(temp);
  PG_FREE_IF_COPY(temp, 0);
  if (result == NULL)
    PG_RETURN_NULL();
  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * Transformation functions
 *****************************************************************************/

/**
 * @brief Set the precision of the fraction of the temporal network point to the
 * number of decimal places.
 */
Temporal *
tnpoint_round(const Temporal *temp, Datum size)
{
  /* We only need to fill these parameters for tfunc_temporal */
  LiftedFunctionInfo lfinfo;
  memset(&lfinfo, 0, sizeof(LiftedFunctionInfo));
  lfinfo.func = (varfunc) &datum_npoint_round;
  lfinfo.numparam = 1;
  lfinfo.param[0] = size;
  lfinfo.restype = temp->temptype;
  lfinfo.tpfunc_base = NULL;
  lfinfo.tpfunc = NULL;
  Temporal *result = tfunc_temporal(temp, &lfinfo);
  return result;
}

PG_FUNCTION_INFO_V1(Tnpoint_round);
/**
 * @ingroup mobilitydb_temporal_transf
 * @brief Set the precision of the fraction of the temporal network point to the
 * number of decimal places.
 * @sqlfunc round()
 */
PGDLLEXPORT Datum
Tnpoint_round(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Datum size = PG_GETARG_DATUM(1);
  Temporal *result = tnpoint_round(temp, size);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_POINTER(result);
}

/*****************************************************************************
 * Accessor functions
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tnpoint_positions);
/**
 * @ingroup mobilitydb_temporal_accessor
 * @brief Return the network segments covered by the temporal network point
 * @sqlfunc positions()
 */
PGDLLEXPORT Datum
Tnpoint_positions(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  int count;
  Nsegment **segments = tnpoint_positions(temp, &count);
  ArrayType *result = nsegmentarr_array(segments, count);
  pfree_array((void **) segments, count);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_POINTER(result);
}

PG_FUNCTION_INFO_V1(Tnpoint_route);
/**
 * @ingroup mobilitydb_temporal_accessor
 * @brief Return the route of a temporal network point
 * @sqlfunc route()
 */
PGDLLEXPORT Datum
Tnpoint_route(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  int64 result = tnpoint_route(temp);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_INT64(result);
}

PG_FUNCTION_INFO_V1(Tnpoint_routes);
/**
 * @ingroup mobilitydb_temporal_accessor
 * @brief Return the array of routes of a temporal network point
 * @sqlfunc routes()
 */
PGDLLEXPORT Datum
Tnpoint_routes(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  int count;
  int64 *routes = tnpoint_routes(temp, &count);
  ArrayType *result = int64arr_array(routes, count);
  pfree(routes);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_POINTER(result);
}

/*****************************************************************************/

/*****************************************************************************
 * TnPoint Complex
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Tgeompoint_to_TnpointComplex);
/**
 * Cast a temporal geometric point as a temporal network point
 */
PGDLLEXPORT Datum
Tgeompoint_to_TnpointComplex(PG_FUNCTION_ARGS)
{
    Temporal *temp = PG_GETARG_TEMPORAL_P(0);
    Datum result = tgeompoint_tnpointComplex(fcinfo, temp);
    PG_FREE_IF_COPY(temp, 0);
    // if (result == NULL)
    //     PG_RETURN_NULL();
    PG_RETURN_DATUM(result);
    // PG_RETURN_POINTER(DatumGetPointer(result));
    // PG_RETURN_POINTER(result);
    //PG_RETURN_NULL();
}

// /**
//  * @brief Parse a temporal value from the buffer (dispatch function).
//  *
//  * @param[in] str Input string
//  */
// Temporal *
// temporal_parse_c(char **str)
// {
//   //p_whitespace(str);
//     Temporal *result_tn = NULL;  /* keep compiler quiet */
//     Temporal *result_tg = NULL;  /* keep compiler quiet */
    
//   // /* Starts with "Interp=Stepwise;" */
//   // if (strncasecmp(*str, "Interp=Stepwise;", 16) == 0)
//   // {
//   //   /* Move str after the semicolon */
//   //   *str += 16;
//   //   linear = false;
//   // }
//   //find the comma
//   int comma_position = NULL;
//   int counter = 0;
//   for(i=0;i++){
//     if(**str[i]==NULL)
//         break;
//     else if(**str[i]=='","'){
//         comma_postion=i;
//     }
//     counter++;
//   }

//   //get from postion 0 to comma_position as tnpoints
//   //"{[NPoint(1,0)@2000-01-01 00:00:00+01]}"
  
//   if (**str != '{' && **str != '[' && **str != '(')
//     //tinstant_parse

//   else if (**str == '{')
//   {
//     // next_position_on_str
//     if (**str == '[' || **str == '(')
//     {
//       // *str = bak;
//       result = (Temporal *) tsequenceset_parse(str, temptype, linear);
//     }
//     else
//     {
//       // *str = bak;
//       result = (Temporal *) tinstantset_parse(str, temptype);
//     }
//   }
//   //tinstantset_parse
//   //tsequenceset_parse

//   //get from comma_position as tgeompoints
//   //"{[0101000020E61000005A17128B7337114047D4FA1175654940@2000-01-01 00:00:00+01, 0101000020E6100000680355489E581140F2333F4860744A40@2000-01-02 00:00:00+01)}"
//     if (**str != '{' && **str != '[' && **str != '(')
//     //tinstant_parse

//   else if (**str == '{')
//   {
//     // next_position_on_str
//     if (**str == '[' || **str == '(')
//     {
//       // *str = bak;
//       result = (Temporal *) tsequenceset_parse(str, temptype, linear);
//     }
//     else
//     {
//       // *str = bak;
//       result = (Temporal *) tinstantset_parse(str, temptype);
//     }
//   }

//   // if (**str != '{' && **str != '[' && **str != '(')
//   //   result = (Temporal *) tinstant_parse(str, temptype, true, true);
//   // else if (**str == '[' || **str == '(')
//   //   result = (Temporal *) tsequence_parse(str, temptype, linear, true, true);
//   // else if (**str == '{')
//   // {
//   //   char *bak = *str;
//   //   p_obrace(str);
//   //   p_whitespace(str);
//   //   if (**str == '[' || **str == '(')
//   //   {
//   //     *str = bak;
//   //     result = (Temporal *) tsequenceset_parse(str, temptype, linear);
//   //   }
//   //   else
//   //   {
//   //     *str = bak;
//   //     result = (Temporal *) tinstantset_parse(str, temptype);
//   //   }
//   // }
//   // return result;
// }

/**
 * @brief Cast a temporal network point as a temporal geometric point.
 */
Temporal *
tnpoint_tgeompoint_c(const Temporal *temp)
{
  Temporal *result;
  ensure_valid_tempsubtype(temp->subtype);
  printf("ok\n");
  fflush(stdout);
  if (temp->subtype == TINSTANT)
    result = (Temporal *) tnpointinst_tgeompointinst((TInstant *) temp);
  else if (temp->subtype == TINSTANTSET)
    result = (Temporal *) tnpointinstset_tgeompointinstset((TInstantSet *) temp);
  else if (temp->subtype == TSEQUENCE)
    result = (Temporal *) tnpointseq_tgeompointseq((TSequence *) temp);
  else /* temp->subtype == TSEQUENCESET */
    result = (Temporal *) tnpointseqset_tgeompointseqset((TSequenceSet *) temp);
  return result;
}

 Temporal* datum_from_heap_tupple(HeapTupleHeader hth, char *attribute, bool *flag){
    return (Temporal *)GetAttributeByName(hth, attribute, &flag);
}

PG_FUNCTION_INFO_V1(TnpointComplex_to_tgeompoint);
/**
 * @ingroup mobilitydb_temporal_cast
 * @brief Cast a temporal complex network point as a temporal geometric point
 * @sqlfunc tgeompoint()
 */
PGDLLEXPORT Datum
TnpointComplex_to_tgeompoint(PG_FUNCTION_ARGS)
{
  // Temporal *temp = PG_GETARG_TEMPORAL_P(0);
    // Temporal *result = tnpoint_tgeompoint(temp);
  Temporal *result = NULL;
  //PG_FREE_IF_COPY(temp, 0);
  //PG_RETURN_POINTER(result);

    // HeapTupleHeader input = PG_GETARG_HEAPTUPLEHEADER(0);
    // HeapTupleHeader input1 = PG_GETARG_HEAPTUPLEHEADER(0);
    bool isnull01;
    bool isnull02;

    // Datum tnp = GetAttributeByName(input, "tnpoints", &isnull01);
    // Datum tgm = GetAttributeByName(input1, "tgeompoints", &isnull02);

    Temporal* tnp = datum_from_heap_tupple(PG_GETARG_HEAPTUPLEHEADER(0), "tnpoints", &isnull01);
    //Datum tgm = datum_from_heap_tupple(PG_GETARG_HEAPTUPLEHEADER(0), "tgeompoints", &isnull02);

    // tn01 = text_to_cstring(PG_GETARG_TEXT_P(0));
    // tg01 = text_to_cstring(PG_GETARG_TEXT_P(1));
    if(!isnull01){
        printf("tnpoints is not null \n");
        // printf(": %lu \n", tnp);
        // fflush(stdout);
        printf("ok1\n");
        fflush(stdout);
        //Temporal *tnp_temp = (Temporal *) PG_DETOAST_DATUM(tnp);
        // printf("ok2\n");
        // fflush(stdout);
        //printf("%d\n", VARSIZE(tnp_temp));
        //fflush(stdout);
        //printf("%d\n", tnp_temp->subtype);
        //fflush(stdout);
        //cast it to tnpoints
        result = tnpoint_tgeompoint_c(tnp);
        PG_RETURN_POINTER(result);
    }else{
        printf("tnpoints is null \n");
        fflush(stdout);
        PG_RETURN_NULL();
    }
    //


    // TnpointComplex tnc = (TnpointComplex *)PG_DETOAST_DATUM(0);
    
    // HeapTupleHeader tg = PG_GETARG_HEAPTUPLEHEADER(1);
    // HeapTuple *tnc = (HeapTuple *)PG_GETARG_DATUM(0);
    //Datum *tg = PG_GETARG_DATUM(1);
    //Temporal *result_01 = tnpoint_tgeompoint((const Temporal **)tn);
    PG_RETURN_NULL();
}

PG_FUNCTION_INFO_V1(intComplex_to_int);
/**
 * demo
 */
PGDLLEXPORT Datum
intComplex_to_int(PG_FUNCTION_ARGS){

    HeapTupleHeader input = PG_GETARG_HEAPTUPLEHEADER(0);
    bool isnull01;
    bool isnull02;

    int i = GetAttributeByName(input, "numbers", &isnull01);
    char c = GetAttributeByName(input, "character", &isnull02);
    printf("results: %d and %c \n", i, c);
    fflush(stdout);

    if(isnull01){
        PG_RETURN_NULL();
    }
    PG_RETURN_INT32(i);

    // Temporal* tnp = datum_from_heap_tupple(PG_GETARG_HEAPTUPLEHEADER(0), "tnpoints", &isnull01);
    //Datum tgm = datum_from_heap_tupple(PG_GETARG_HEAPTUPLEHEADER(0), "tgeompoints", &isnull02);
}

PG_FUNCTION_INFO_V1(intPlus_to_int);
/**
 * demo
 */
PGDLLEXPORT Datum
intPlus_to_int(PG_FUNCTION_ARGS){

    HeapTupleHeader input = PG_GETARG_HEAPTUPLEHEADER(0);
    bool isnull01;
    bool isnull02;

    int i = GetAttributeByName(input, "numbers", &isnull01);
    Datum c = GetAttributeByName(input, "tgeompoints", &isnull02);

    Temporal *d = DatumGetPointer(GetAttributeByName(input, "tgeompoints", &isnull02));
    Datum b = PointerGetDatum((Temporal *)PG_DETOAST_DATUM(c)); 
    Temporal *a = (Temporal *)PG_DETOAST_DATUM(c);
    if(isnull02){
        PG_RETURN_NULL();
    }
    PG_RETURN_DATUM(b);

    // Temporal* tnp = datum_from_heap_tupple(PG_GETARG_HEAPTUPLEHEADER(0), "tnpoints", &isnull01);
    //Datum tgm = datum_from_heap_tupple(PG_GETARG_HEAPTUPLEHEADER(0), "tgeompoints", &isnull02);
}


PG_FUNCTION_INFO_V1(tnpointandtgeompoint_to_tgeompoint);
/**
 * demo
 */
PGDLLEXPORT Datum
tnpointandtgeompoint_to_tgeompoint(PG_FUNCTION_ARGS){

    Temporal *tnpoints = PG_GETARG_TEMPORAL_P(0);
    Temporal *tgeompoints = PG_GETARG_TEMPORAL_P(0);
    bool isnull01;
    bool isnull02;

    Temporal *result01 = tnpoint_tgeompoint(tnpoints);
    //Temporal merge two tgeompoints
    //result01 and tgeompoints as Temporal *result

    PG_FREE_IF_COPY(tnpoints, 0);
    // PG_RETURN_NULL();
    PG_RETURN_DATUM(PointerGetDatum(result01));
}