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

/* C */
#include <assert.h>
/* MobilityDB */
#include <meos.h>
#include <meos_internal.h>
#include "general/temporal_parser.h"
#include "general/temporal_util.h"
#include "general/lifting.h"
#include "point/tpoint_spatialfuncs.h"
#include "npoint/tnpoint_static.h"
#include "npoint/tnpoint_parser.h"

/* Heap Tuple */
#include <funcapi.h>

/*****************************************************************************
 * Input/output functions
 *****************************************************************************/


/*****************************************************************************
 * Cast functions
 *****************************************************************************/

/**
 * @brief Cast a temporal network point as a temporal geometric point.
 */
TInstant *
tnpointinst_tgeompointinst(const TInstant *inst)
{
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  GSERIALIZED *geom = npoint_geom(np);
  TInstant *result = 
    tinstant_make(PointerGetDatum(geom), T_TGEOMPOINT, inst->t);
  pfree(geom);
  return result;
}

/**
 * @brief Cast a temporal network point as a temporal geometric point.
 */
TInstantSet *
tnpointinstset_tgeompointinstset(const TInstantSet *is)
{
  TInstant **instants = palloc(sizeof(TInstant *) * is->count);
  for (int i = 0; i < is->count; i++)
  {
    const TInstant *inst = tinstantset_inst_n(is, i);
    instants[i] = tnpointinst_tgeompointinst(inst);
  }
  TInstantSet *result = tinstantset_make_free(instants, is->count, MERGE_NO);
  return result;
}

/**
 * @brief Cast a temporal network point as a temporal geometric point.
 */
TSequence *
tnpointseq_tgeompointseq(const TSequence *seq)
{
  TInstant **instants = palloc(sizeof(TInstant *) * seq->count);
  Npoint *np = DatumGetNpointP(tinstant_value(tsequence_inst_n(seq, 0)));
  /* We are sure line is not empty */
  GSERIALIZED *line = route_geom(np->rid);
  int srid = gserialized_get_srid(line);
  LWLINE *lwline = (LWLINE *) lwgeom_from_gserialized(line);
  for (int i = 0; i < seq->count; i++)
  {
    const TInstant *inst = tsequence_inst_n(seq, i);
    np = DatumGetNpointP(tinstant_value(inst));
    POINTARRAY *opa = lwline_interpolate_points(lwline, np->pos, 0);
    LWGEOM *lwpoint;
    assert(opa->npoints <= 1);
    lwpoint = lwpoint_as_lwgeom(lwpoint_construct(srid, NULL, opa));
    Datum point = PointerGetDatum(geo_serialize(lwpoint));
    instants[i] = tinstant_make(point, T_TGEOMPOINT, inst->t);
    pfree(DatumGetPointer(point));
  }
  TSequence *result = tsequence_make_free(instants, seq->count,
    seq->period.lower_inc, seq->period.upper_inc,
    MOBDB_FLAGS_GET_LINEAR(seq->flags), false);
  pfree(DatumGetPointer(line));
  return result;
}

/**
 * @brief Cast a temporal network point as a temporal geometric point.
 */
TSequenceSet *
tnpointseqset_tgeompointseqset(const TSequenceSet *ss)
{
  TSequence **sequences = palloc(sizeof(TSequence *) * ss->count);
  for (int i = 0; i < ss->count; i++)
  {
    const TSequence *seq = tsequenceset_seq_n(ss, i);
    sequences[i] = tnpointseq_tgeompointseq(seq);
  }
  TSequenceSet *result = tsequenceset_make_free(sequences, ss->count, false);
  return result;
}

/**
 * @brief Cast a temporal network point as a temporal geometric point.
 */
Temporal *
tnpoint_tgeompoint(const Temporal *temp)
{
  Temporal *result;
  ensure_valid_tempsubtype(temp->subtype);
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

/*****************************************************************************/

/**
 * @brief Cast a temporal geometric point as a temporal network point.
 */
TInstant *
tgeompointinst_tnpointinst(const TInstant *inst)
{
  GSERIALIZED *gs = DatumGetGserializedP(tinstant_value(inst));
  Npoint *np = geom_npoint(gs);
  if (np == NULL)
    return NULL;
  TInstant *result = tinstant_make(PointerGetDatum(np), T_TNPOINT, inst->t);
  pfree(np);
  return result;
}

/**
 * @brief Cast a temporal geometric point as a temporal network point.
 */
TInstantSet *
tgeompointinstset_tnpointinstset(const TInstantSet *is)
{
  TInstant **instants = palloc(sizeof(TInstant *) * is->count);
  for (int i = 0; i < is->count; i++)
  {
    const TInstant *inst = tinstantset_inst_n(is, i);
    TInstant *inst1 = tgeompointinst_tnpointinst(inst);
    if (inst1 == NULL)
    {
      pfree_array((void **) instants, i);
      return NULL;
    }
    instants[i] = inst1;
  }
  TInstantSet *result = tinstantset_make_free(instants, is->count, MERGE_NO);
  return result;
}

/**
 * @brief Cast a temporal geometric point as a temporal network point.
 */
TSequence *
tgeompointseq_tnpointseq(const TSequence *seq)
{
  TInstant **instants = palloc(sizeof(TInstant *) * seq->count);
  for (int i = 0; i < seq->count; i++)
  {
    const TInstant *inst = tsequence_inst_n(seq, i);
    TInstant *inst1 = tgeompointinst_tnpointinst(inst);
    if (inst1 == NULL)
    {
      pfree_array((void **) instants, i);
      return NULL;
    }
    instants[i] = inst1;
  }
  TSequence *result = tsequence_make_free(instants, seq->count,
    seq->period.lower_inc, seq->period.upper_inc,
    MOBDB_FLAGS_GET_LINEAR(seq->flags), true);
  return result;
}

/**
 * @brief Cast a temporal geometric point as a temporal network point.
 */
TSequenceSet *
tgeompointseqset_tnpointseqset(const TSequenceSet *ss)
{
  TSequence **sequences = palloc(sizeof(TSequence *) * ss->count);
  for (int i = 0; i < ss->count; i++)
  {
    const TSequence *seq = tsequenceset_seq_n(ss, i);
    TSequence *seq1 = tgeompointseq_tnpointseq(seq);
    if (seq1 == NULL)
    {
      pfree_array((void **) sequences, i);
      return NULL;
    }
    sequences[i] = seq1;
  }
  TSequenceSet *result = tsequenceset_make_free(sequences, ss->count, true);
  return result;
}

/**
 * @brief Cast a temporal geometric point as a temporal network point.
 */
Temporal *
tgeompoint_tnpoint(const Temporal *temp)
{
  int32_t srid_tpoint = tpoint_srid(temp);
  int32_t srid_ways = get_srid_ways();
  ensure_same_srid(srid_tpoint, srid_ways);
  Temporal *result;
  ensure_valid_tempsubtype(temp->subtype);
  if (temp->subtype == TINSTANT)
    result = (Temporal *) tgeompointinst_tnpointinst((TInstant *) temp);
  else if (temp->subtype == TINSTANTSET)
    result = (Temporal *) tgeompointinstset_tnpointinstset((TInstantSet *) temp);
  else if (temp->subtype == TSEQUENCE)
    result = (Temporal *) tgeompointseq_tnpointseq((TSequence *) temp);
  else /* temp->subtype == TSEQUENCESET */
    result = (Temporal *) tgeompointseqset_tnpointseqset((TSequenceSet *) temp);
  return result;
}

/*****************************************************************************
 * Accessor functions
 *****************************************************************************/

/**
 * @brief Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointinst_positions(const TInstant *inst)
{
  Nsegment **result = palloc(sizeof(Nsegment *));
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  result[0] = nsegment_make(np->rid, np->pos, np->pos);
  return result;
}

/**
 * @brief Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointinstset_positions(const TInstantSet *is, int *count)
{
  int count1;
  /* The following function removes duplicate values */
  Datum *values = tinstantset_values(is, &count1);
  Nsegment **result = palloc(sizeof(Nsegment *) * count1);
  for (int i = 0; i < count1; i++)
  {
    Npoint *np = DatumGetNpointP(values[i]);
    result[i] = nsegment_make(np->rid, np->pos, np->pos);
  }
  pfree(values);
  *count = count1;
  return result;
}

/**
 * Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointseq_step_positions(const TSequence *seq, int *count)
{
  int count1;
  /* The following function removes duplicate values */
  Datum *values = tsequence_values(seq, &count1);
  Nsegment **result = palloc(sizeof(Nsegment *) * count1);
  for (int i = 0; i < count1; i++)
  {
    Npoint *np = DatumGetNpointP(values[i]);
    result[i] = nsegment_make(np->rid, np->pos, np->pos);
  }
  pfree(values);
  *count = count1;
  return result;
}

/**
 * Return the network segments covered by the temporal network point.
 */
Nsegment *
tnpointseq_linear_positions(const TSequence *seq)
{
  const TInstant *inst = tsequence_inst_n(seq, 0);
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  int64 rid = np->rid;
  double minPos = np->pos, maxPos = np->pos;
  for (int i = 1; i < seq->count; i++)
  {
    inst = tsequence_inst_n(seq, i);
    np = DatumGetNpointP(tinstant_value(inst));
    minPos = Min(minPos, np->pos);
    maxPos = Max(maxPos, np->pos);
  }
  return nsegment_make(rid, minPos, maxPos);
}

/**
 * @brief Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointseq_positions(const TSequence *seq, int *count)
{
  if (MOBDB_FLAGS_GET_LINEAR(seq->flags))
  {
    Nsegment **result = palloc(sizeof(Nsegment *));
    result[0] = tnpointseq_linear_positions(seq);
    *count = 1;
    return result;
  }
  else
    return tnpointseq_step_positions(seq, count);
}

/**
 * Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointseqset_linear_positions(const TSequenceSet *ss, int *count)
{
  Nsegment **segments = palloc(sizeof(Nsegment *) * ss->count);
  for (int i = 0; i < ss->count; i++)
  {
    const TSequence *seq = tsequenceset_seq_n(ss, i);
    segments[i] = tnpointseq_linear_positions(seq);
  }
  Nsegment **result = segments;
  int count1 = ss->count;
  if (count1 > 1)
    result = nsegmentarr_normalize(segments, &count1);
  *count = count1;
  return result;
}

/**
 * Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointseqset_step_positions(const TSequenceSet *ss, int *count)
{
  /* The following function removes duplicate values */
  int newcount;
  Datum *values = tsequenceset_values(ss, &newcount);
  Nsegment **result = palloc(sizeof(Nsegment *) * newcount);
  for (int i = 0; i < newcount; i++)
  {
    Npoint *np = DatumGetNpointP(values[i]);
    result[i] = nsegment_make(np->rid, np->pos, np->pos);
  }
  pfree(values);
  *count = newcount;
  return result;
}

/**
 * @brief Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpointseqset_positions(const TSequenceSet *ss, int *count)
{
  Nsegment **result;
  if (MOBDB_FLAGS_GET_LINEAR(ss->flags))
    result = tnpointseqset_linear_positions(ss, count);
  else
    result = tnpointseqset_step_positions(ss, count);
  return result;
}

/**
 * @brief Return the network segments covered by the temporal network point.
 */
Nsegment **
tnpoint_positions(const Temporal *temp, int *count)
{
  Nsegment **result;
  ensure_valid_tempsubtype(temp->subtype);
  if (temp->subtype == TINSTANT)
  {
    result = tnpointinst_positions((TInstant *) temp);
    *count = 1;
  }
  else if (temp->subtype == TINSTANTSET)
    result = tnpointinstset_positions((TInstantSet *) temp, count);
  else if (temp->subtype == TSEQUENCE)
    result = tnpointseq_positions((TSequence *) temp, count);
  else /* temp->subtype == TSEQUENCESET */
    result = tnpointseqset_positions((TSequenceSet *) temp, count);
  return result;
}

/*****************************************************************************/

/**
 * @brief Return the route of the temporal network point.
 */
int64
tnpointinst_route(const TInstant *inst)
{
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  return np->rid;
}

/**
 * @brief Return the route of a temporal network point.
 */
int64
tnpoint_route(const Temporal *temp)
{
  if (temp->subtype != TINSTANT && temp->subtype != TSEQUENCE)
    elog(ERROR, "Input must be a temporal instant or a temporal sequence");

  const TInstant *inst = (temp->subtype == TINSTANT) ?
    (const TInstant *) temp : tsequence_inst_n((const TSequence *) temp, 0);
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  return np->rid;
}

/**
 * @brief Return the array of routes of a temporal network point
 */
int64 *
tnpointinst_routes(const TInstant *inst)
{
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  int64 *result = palloc(sizeof(int64));
  result[0]= np->rid;
  return result;
}

/**
 * @brief Return the array of routes of a temporal network point
 */
int64 *
tnpointinstset_routes(const TInstantSet *is)
{
  int64 *result = palloc(sizeof(int64) * is->count);
  for (int i = 0; i < is->count; i++)
  {
    const TInstant *inst = tinstantset_inst_n(is, i);
    Npoint *np = DatumGetNpointP(tinstant_value(inst));
    result[i] = np->rid;
  }
  return result;
}

/**
 * @brief Return the array of routes of a temporal network point
 */
int64 *
tnpointseq_routes(const TSequence *seq)
{
  const TInstant *inst = tsequence_inst_n(seq, 0);
  Npoint *np = DatumGetNpointP(tinstant_value(inst));
  int64 *result = palloc(sizeof(int64));
  result[0]= np->rid;
  return result;
}

/**
 * @brief Return the array of routes of a temporal network point
 */
int64 *
tnpointseqset_routes(const TSequenceSet *ss)
{
  int64 *result = palloc(sizeof(int64) * ss->count);
  for (int i = 0; i < ss->count; i++)
  {
    const TSequence *seq = tsequenceset_seq_n(ss, i);
    const TInstant *inst = tsequence_inst_n(seq, 0);
    Npoint *np = DatumGetNpointP(tinstant_value(inst));
    result[i] = np->rid;
  }
  return result;
}

/**
 * @brief Return the array of routes of a temporal network point
 */
int64 *
tnpoint_routes(const Temporal *temp, int *count)
{
  int64 *result;
  ensure_valid_tempsubtype(temp->subtype);
  if (temp->subtype == TINSTANT)
  {
    result = tnpointinst_routes((TInstant *) temp);
    *count = 1;
  }
  else if (temp->subtype == TINSTANTSET)
  {
    result = tnpointinstset_routes((TInstantSet *) temp);
    *count = ((TInstantSet *) temp)->count;
  }
  else if (temp->subtype == TSEQUENCE)
  {
    result = tnpointseq_routes((TSequence *) temp);
    *count = 1;
  }
  else /* temp->subtype == TSEQUENCESET */
  {
    result = tnpointseqset_routes((TSequenceSet *) temp);
    *count = ((TSequenceSet *) temp)->count;
  }
  return result;
}

/*****************************************************************************/

/*****************************************************************************
 * TnPoint Complex
 *****************************************************************************/
/**
 * Free a C array of pointers
 */
void
pfree_array_c(void **array, int count)
{
  for (int i = 0; i < count; i++)
    pfree(array[i]);
  return;
}


/**
 * @ingroup
 * @brief Return a Heap Tuple with temporal network point and tgeompoint tuple.
 */
Datum
tgeompoint_tnpointHeapTupleDatum(FunctionCallInfo fcinfo,const Temporal *result_tnpoints, const Temporal *result_tgeompoints)
{
    /* Build a tuple description for the function output */
    TupleDesc resultTupleDesc;
    //FunctionCallInfo fcinfo
    get_call_result_type(fcinfo, NULL, &resultTupleDesc);
    BlessTupleDesc(resultTupleDesc);

    /* Construct the result */
    HeapTuple resultTuple;
    bool result_is_null[2]; /* needed to say no value is null */
    Datum result_values[2]; /* used to construct the composite return value */
    Datum result_combined; /* the actual composite return value */

    /* Store Tnpoints */
    result_values[0] = PointerGetDatum(result_tnpoints);
    /* Store Exceptions */
    result_values[1] = PointerGetDatum(result_tgeompoints);

    /* Setting Boolean value */
    if (result_tnpoints == NULL && result_tgeompoints == NULL){
        result_is_null[0] = 1;
        result_is_null[1] = 1;
    }else if(result_tnpoints != NULL && result_tgeompoints == NULL){
        result_is_null[0] = 0;
        result_is_null[1] = 1;
    }else if(result_tnpoints == NULL && result_tgeompoints != NULL){
        result_is_null[0] = 1;
        result_is_null[1] = 0;
    }else {
        result_is_null[0] = 0;
        result_is_null[1] = 0;
    }

    /* Form tuple and return */
    resultTuple = heap_form_tuple(resultTupleDesc, result_values, result_is_null);
    // result_combined = DatumGetPointer(HeapTupleGetDatum(resultTuple));
    result_combined = HeapTupleGetDatum(resultTuple);
    //printf("Inside Heap Tuple Function Block!!\n");
    // fflush(stdout);
    return result_combined;
}

/**
 * @ingroup
 * @brief Cast a Tgeompoint instant into temporal network point complex tuple.
 */
Datum
tgeompointinst_tnpointComplexinst(FunctionCallInfo fcinfo, const TInstant *inst)
{
    GSERIALIZED *gs = DatumGetGserializedP(tinstant_value(inst));
    Npoint *np = geom_npoint(gs);
    Datum result_c;
    if (np == NULL){
      // result_c = PointerGetDatum(tgeompoint_tnpointHeapTupleDatum(fcinfo, NULL, (Temporal *) inst));
      result_c = tgeompoint_tnpointHeapTupleDatum(fcinfo, NULL, (Temporal *) inst);
      return result_c;
    }
    TInstant *result = tinstant_make(PointerGetDatum(np), T_TNPOINT, inst->t);
    pfree(np);
    // result_c = PointerGetDatum(tgeompoint_tnpointHeapTupleDatum(fcinfo, (Temporal *) result, (Temporal *) inst));
    result_c = tgeompoint_tnpointHeapTupleDatum(fcinfo, (Temporal *) result, (Temporal *) inst);
    return result_c;
}


/**
 * @ingroup
 * @brief Return a reduced copy of the temporal instant value.
 */
TInstant **
tinstant_reduce_and_copy(TInstant **inst, const int new_length)
{
    if(sizeof(TInstant *) * new_length >= VARSIZE(inst)){
        return (TInstant **)inst;
    } else {
        int new_size = sizeof(TInstant *) * new_length;
        TInstant **result = palloc0(new_size);
        memcpy(result, inst, new_size);
        return result;
    }
}

// /**
//  * @ingroup
//  * @brief Return a reduced copy of the temporal sequence value.
//  */
// TSequence *
// tsequence_reduce_and_copy(const TSequence *inst, const int new_length)
// {
//     if(sizeof(TSequence *) * new_length >= VARSIZE(inst)){
//         return inst;
//     } else {
//         int new_size = sizeof(TSequence *) * new_length;
//         TSequence *result = palloc0(new_size);
//         memcpy(result, inst, new_size);
//         return result;
//     }
// }

/**
 * @ingroup libmeos_temporal_constructor
 * @brief Construct a temporal instant set from an array of temporal instants
 * and free the array and the instants after the creation.
 *
 * @param[in] instants Array of instants
 * @param[in] count Number of elements in the array
 * @param[in] merge True when overlapping instants are allowed as required in
 * merge operations
 * @see tinstantset_make
 */
TInstantSet *
tinstantset_make_free_c(TInstant **instants, int count, bool merge)
{
  if (count == 0)
  {
    pfree(instants);
    return NULL;
  }
  TInstantSet *result = tinstantset_make((const TInstant **) instants,
    count, merge);
  // pfree_array((void **) instants, count);
  return result;
}

/**
 * @ingroup libmeos_temporal_constructor
 * @brief Construct a temporal sequence set from an array of temporal
 * sequences and free the array and the sequences after the creation.
 *
 * @param[in] sequences Array of sequences
 * @param[in] count Number of elements in the array
 * @param[in] normalize True when the resulting value should be normalized.
 * @see tsequenceset_make
 */
TSequenceSet *
tsequenceset_make_free_c(TSequence **sequences, int count, bool normalize)
{
  if (count == 0)
  {
    pfree(sequences);
    return NULL;
  }
  TSequenceSet *result = tsequenceset_make((const TSequence **) sequences,
    count, normalize);
  // pfree_array((void **) sequences, count);
  return result;
}

/**
 * @ingroup
 * @brief Cast a Tgeompoint instant into temporal network point complex tuple.
 */
Datum
tgeompointinstset_tnpointComplexinstset(FunctionCallInfo fcinfo,const TInstantSet *is)
{
    TInstant **instants = palloc0(sizeof(TInstant *) * is->count);
    TInstant **instants_tgeom = palloc0(sizeof(TInstant *) * is->count);
    int count1 = 0;
    int count2 = 0;
    for (int i = 0; i < is->count; i++)
    {
        const TInstant *inst = tinstantset_inst_n(is, i);
        TInstant *inst1 = tgeompointinst_tnpointinst(inst);
        if (inst1 == NULL)
        {
            instants_tgeom[count2] = (TInstant *)inst; //assigning TInstant cast as "assigning const to non-const variable throws warning"
            count2++;
        }else{
            instants[count1] = inst1;
            count1++;
        }
    }
    TInstant **instants_r = tinstant_reduce_and_copy(instants,count1);
    TInstant **instants_tgeom_r = tinstant_reduce_and_copy(instants_tgeom,count2);
    TInstantSet *result = tinstantset_make_free_c(instants_r, count1, MERGE_NO);
    TInstantSet *result_geom = tinstantset_make_free_c(instants_tgeom_r, count2, MERGE_NO);
    Datum result_combined = tgeompoint_tnpointHeapTupleDatum(fcinfo,(Temporal *)result,(Temporal *)result_geom);
    // count1=0;
    // count2=0;
    // pfree_array((void **) instants, is->count);
    // pfree_array((void **) instants_tgeom, is->count);
    return result_combined;
}

/**
 * @ingroup libmeos_temporal_constructor
 * @brief Construct a temporal sequence from an array of temporal instants.
 *
 * @param[in] instants Array of instants
 * @param[in] output_seq_count Number of elements in the output array
 * @param[in] input_seq_count Number of elements in the input array
 * @param[in] lower_inc,upper_inc True when the respective bound is inclusive
 * @param[in] linear True when the interpolation is linear
 * @param[in] normalize True when the resulting value should be normalized
 * @param[in] iterations is 0 and input_seq_count for the first and last iterations respectively
 * @sqlfunc tbool_seq(), tint_seq(), tfloat_seq(), ttext_seq(), etc.
 */
TSequence *
tsequence_make_c(const TInstant **instants, int output_seq_count, int input_seq_count, bool lower_inc,
  bool upper_inc, bool linear, bool normalize, int iterations)
{
  if(output_seq_count==1)
    //instantenous sequence must have inclusive bounds
    return tsequence_make(instants, output_seq_count,true,true,
      linear,normalize);
  else{
    if(iterations==0 && input_seq_count==1){
    //first case could be removed as it is included in above conditions, as output_seq_count will be 1 if input_seq_count is 1
    return tsequence_make(instants,output_seq_count, lower_inc, upper_inc,
      linear,normalize);
    }else if(iterations==0 && input_seq_count>1){
      return tsequence_make(instants,output_seq_count,lower_inc,true,
        linear,normalize); 
    }else if(iterations>0 && iterations==input_seq_count-1){
      return tsequence_make(instants,output_seq_count,true, upper_inc,
        linear,normalize); 
    }else{
      return tsequence_make(instants, output_seq_count,true,true,
        linear,normalize);
    }  
  }
}

/**
 * @ingroup libmeos_temporal_cast
 * @brief Cast a temporal geometric sequence as a temporal network complex type.
 */
Datum
tgeompointseq_tnpointComplexseq(FunctionCallInfo fcinfo, const TSequence *seq)
{
    TSequence **sequences_tn = palloc0(sizeof(TSequence *) * seq->count);
    TSequence **sequences_tg = palloc0(sizeof(TSequence *) * seq->count);
    TInstant **instants_temp = palloc0(sizeof(TInstant *) * seq->count);
    // TSequence *sequence_temp = palloc0(sizeof(TSequence *) * seq->count);
    TInstant *buffer_point = palloc0(sizeof(TInstant *));

    int instants_temp_count = 0;

    bool prev_is_tgeom;
    bool this_is_tgeom;
    int count_tn = 0;
    int count_tg = 0;

    for (int i = 0; i < seq->count; i++)
    {
        const TInstant *inst = tsequence_inst_n(seq, i);
        TInstant *inst1 = tgeompointinst_tnpointinst(inst);
        if(i==0){
              if(inst1==NULL){
                instants_temp[instants_temp_count]=(TInstant *)inst;
                // instants_temp[instants_temp_count]=tsequence_inst_n(seq, i);
                // memcpy(instants_temp, inst, sizeof(TInstant *));
                instants_temp_count++;
                //for next iteration
                prev_is_tgeom = true;
              }
              else{
                // instants_temp[instants_temp_count]=inst1;
                instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                // memcpy(instants_temp, inst1, sizeof(TInstant *));
                instants_temp_count++;
                // printf("Debug #first instants");
                // fflush(stdout);
                //for next iteration
                prev_is_tgeom = false;

                //storing this tgeompoint to be used if there is no tnpoint in next iteration
                buffer_point = (TInstant *)inst;
              }
        }else{
              if (inst1 == NULL){
                          this_is_tgeom = true;
                          if(prev_is_tgeom == true){
                              //pass this inst to same sequence of tgeompoints
                              instants_temp[instants_temp_count]= (TInstant *)inst;
                              instants_temp_count++;

                              //set prev_is_tgeom = true for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }else{
                              //changes
                              // printf("Debug #second instant: %d \n", instants_temp_count);
                              // fflush(stdout);
                              // if(i==0 && seq->count==1){
                              //   sequences_tn[count_tn] = tsequence_make((const TInstant **)instants_temp, 
                              //     instants_temp_count,
                              //     seq->period.lower_inc, seq->period.upper_inc,
                              //     MOBDB_FLAGS_GET_LINEAR(seq->flags),true);    
                              // }else if(i==0 && seq->count>1){
                              //   sequences_tn[count_tn] = tsequence_make((const TInstant **)instants_temp, 
                              //     instants_temp_count,
                              //     seq->period.lower_inc, true,
                              //     MOBDB_FLAGS_GET_LINEAR(seq->flags),true); 
                              // }else if(i>0 && i==seq->count-1){
                              //   sequences_tn[count_tn] = tsequence_make((const TInstant **)instants_temp, 
                              //     instants_temp_count,
                              //     true, seq->period.upper_inc,
                              //     MOBDB_FLAGS_GET_LINEAR(seq->flags),true); 
                              // }else{
                              //   sequences_tn[count_tn] = tsequence_make((const TInstant **)instants_temp, instants_temp_count,
                              //     true,true,
                              //     MOBDB_FLAGS_GET_LINEAR(seq->flags),true);
                              // }
                            sequences_tn[count_tn]=tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                              count_tn++;

                              //free the temporary array of instants after making the sequence
                              //pfree_array_c((void **) instants_temp, instants_temp_count);
                              //initialize the temporary array of instants counter to 0
                              instants_temp_count = 0;

                              //pass the last point of tnpoint sequence as the first point of next tgeompoint sequence
                              //this point will be in both tnpoint and tgeompoint sequence
                              //earlier point
                              instants_temp[instants_temp_count]=buffer_point;
                              instants_temp_count++;

                              //this point
                              instants_temp[instants_temp_count]=(TInstant *)inst;
                              instants_temp_count++;
                              //set prev_is_tgeom = true for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }
              }else{
                          //storing this tgeompoint to be used if there is no tnpoint in next iteration
                          buffer_point = (TInstant *)inst;

                          //seq_instant_tnpoint = tsequence_make((const TInstant **)temp_instant_tnpoint, 1,true,true,true,true);
                          this_is_tgeom = false;
                          if(prev_is_tgeom == false){
                              //if(this routeid == previous routeid) else different sequences
                              if(tnpointinst_route(inst1)== tnpointinst_route(instants_temp[instants_temp_count-1])){
                                //pass this inst1 to same sequence of tnpoints
                                instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                                instants_temp_count++;
                              }else{
                                sequences_tn[count_tn] = tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                                count_tn++;
                                instants_temp_count =0;
                                //
                                instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                                instants_temp_count++;
                              }                              

                              //set prev_is_tgeom = false for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }else{
                              //pass the first point of tnpoint sequence also as the last point of previous tgeompoint sequence
                              //this point will be in both tnpoint and tgeompoint sequence
                              instants_temp[instants_temp_count]=(TInstant *)inst;
                              //instants_temp_count++;

                              // printf("Debugging : %d \n", instants_temp_count);
                              // fflush(stdout);
                              // sequences_tg[count_tg] = tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                              sequences_tg[count_tg] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                              count_tg++;

                              //free the temporary array of instants after making the sequence

                              //pfree_array_c((void **) instants_temp, instants_temp_count);
                              //initialize the temporary array of instants counter to 0
                              instants_temp_count = 0;

                              //this point
                              instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                              instants_temp_count++;
                              //set prev_is_tgeom = false for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }
              }
        }
        
        if(i==seq->count-1 && instants_temp_count>0){
            if(inst1 == NULL){
                // printf("Debug #3: %d \n", instants_temp_count);
                // fflush(stdout);
                sequences_tg[count_tg] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                // = tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                count_tg++;

                //free the temporary array of instants after making the sequence
                // pfree_array((void **) instants_temp, instants_temp_count);
                //initialize the temporary array of instants counter to 0
                instants_temp_count = 0;
            }else{
                sequences_tn[count_tn] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                // = tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                count_tn++;

                //free the temporary array of instants after making the sequence
                // pfree_array((void **) instants_temp, instants_temp_count);
                //initialize the temporary array of instants counter to 0
                instants_temp_count = 0;
            }
        }
    }
    // printf("Two Counts: %d and %d \n", count_tn, count_tg);
    //TSequence **sequences_tn_r = tsequencearray_reduce_and_copy(sequences_tn,count_tn);
    //TSequence **sequences_tg_r = tsequencearray_reduce_and_copy(sequences_tg,count_tg);
    TSequenceSet *result_tnpoints = tsequenceset_make_free_c(sequences_tn, count_tn, MERGE_NO);
    TSequenceSet *result_tgeom = tsequenceset_make_free_c(sequences_tg, count_tg, MERGE_NO);

    Datum result_combined=tgeompoint_tnpointHeapTupleDatum(fcinfo,(Temporal *)result_tnpoints,(Temporal *)result_tgeom);

    //PG_FREE_IF_COPY(sequence_temp, 0);
    //pfree(instants_temp);
    // pfree(buffer_point);
    return result_combined;
}

/**
 * @brief Cast a temporal geometric point as a temporal network point.
 */
Datum
tgeompointseqset_tnpointComplexseqset(FunctionCallInfo fcinfo, const TSequenceSet *ss)
{
  //TSequence **sequences = palloc(sizeof(TSequence *) * ss->count);
  // TSequence **sequences_tn = palloc0(sizeof(TSequence *) * ss->count);
  // TSequence **sequences_tg = palloc0(sizeof(TSequence *) * ss->count);
  int total_instant_count = 0;
  for (int i = 0; i < ss->count; i++){
    const TSequence *seq_ = tsequenceset_seq_n(ss, i);
    total_instant_count+=seq_->count;
  }
  TSequence **sequences_tn = palloc0(sizeof(TSequence *) * total_instant_count);
  TSequence **sequences_tg = palloc0(sizeof(TSequence *) * total_instant_count);
  
  int count_tn = 0;
  int count_tg = 0;

  for (int i = 0; i < ss->count; i++)
  {
    const TSequence *seq = tsequenceset_seq_n(ss, i);
    int instants_temp_count = 0;
    bool prev_is_tgeom;
    bool this_is_tgeom;
    TInstant **instants_temp = palloc0(sizeof(TInstant *) * seq->count);
    TInstant *buffer_point = palloc0(sizeof(TInstant *));

    // HeapTuple read_heap_tuple = (HeapTuple) tgeompointseq_tnpointComplexseq(fcinfo, seq);
    // printf("HeapTupple Size: %d \n", read_heap_tuple->t_len);
    // //headtupple -> t_len
    // //headtupple -> t_data;
    // //headtupple -> t_self;
    // fflush(stdout);

    for (int i = 0; i < seq->count; i++)
    {
        const TInstant *inst = tsequence_inst_n(seq, i);
        TInstant *inst1 = tgeompointinst_tnpointinst(inst);
        if(i==0){
              if(inst1==NULL){
                instants_temp[instants_temp_count]=(TInstant *)inst;
                // instants_temp[instants_temp_count]=tsequence_inst_n(seq, i);
                // memcpy(instants_temp, inst, sizeof(TInstant *));
                instants_temp_count++;
                //for next iteration
                prev_is_tgeom = true;
              }
              else{
                // instants_temp[instants_temp_count]=inst1;
                instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                // memcpy(instants_temp, inst1, sizeof(TInstant *));
                instants_temp_count++;
                // printf("Debug #first instants");
                // fflush(stdout);
                //for next iteration
                prev_is_tgeom = false;

                //storing this tgeompoint to be used if there is no tnpoint in next iteration
                buffer_point = (TInstant *)inst;
              }
        }else{
              if (inst1 == NULL){
                          this_is_tgeom = true;
                          if(prev_is_tgeom == true){
                              //pass this inst to same sequence of tgeompoints
                              instants_temp[instants_temp_count]= (TInstant *)inst;
                              instants_temp_count++;

                              //set prev_is_tgeom = true for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }else{
                              //changes
                              // printf("Debug #second instant: %d \n", instants_temp_count);
                              // fflush(stdout);
                              sequences_tn[count_tn] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                              //tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                              count_tn++;

                              //free the temporary array of instants after making the sequence
                              //pfree_array_c((void **) instants_temp, instants_temp_count);
                              //initialize the temporary array of instants counter to 0
                              instants_temp_count = 0;

                              //pass the last point of tnpoint sequence as the first point of next tgeompoint sequence
                              //this point will be in both tnpoint and tgeompoint sequence
                              //earlier point
                              instants_temp[instants_temp_count]=buffer_point;
                              instants_temp_count++;

                              //this point
                              instants_temp[instants_temp_count]=(TInstant *)inst;
                              instants_temp_count++;
                              //set prev_is_tgeom = true for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }
              }else{
                          //storing this tgeompoint to be used if there is no tnpoint in next iteration
                          buffer_point = (TInstant *)inst;

                          //seq_instant_tnpoint = tsequence_make((const TInstant **)temp_instant_tnpoint, 1,true,true,true,true);
                          this_is_tgeom = false;
                          if(prev_is_tgeom == false){
                              //if(this routeid == previous routeid) else different sequences
                              if(tnpointinst_route(inst1)== tnpointinst_route(instants_temp[instants_temp_count-1])){
                                //pass this inst1 to same sequence of tnpoints
                                instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                                instants_temp_count++;
                              }else{
                                sequences_tn[count_tn] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                                //tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                                count_tn++;
                                instants_temp_count =0;
                                //
                                instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                                instants_temp_count++;
                              }                              

                              //set prev_is_tgeom = false for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }else{
                              //pass the first point of tnpoint sequence also as the last point of previous tgeompoint sequence
                              //this point will be in both tnpoint and tgeompoint sequence
                              instants_temp[instants_temp_count]=(TInstant *)inst;
                              //instants_temp_count++;

                              // printf("Debugging : %d \n", instants_temp_count);
                              // fflush(stdout);
                              sequences_tg[count_tg] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                              //tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                              count_tg++;

                              //free the temporary array of instants after making the sequence

                              //pfree_array_c((void **) instants_temp, instants_temp_count);
                              //initialize the temporary array of instants counter to 0
                              instants_temp_count = 0;

                              //this point
                              instants_temp[instants_temp_count]=tgeompointinst_tnpointinst(inst);
                              instants_temp_count++;
                              //set prev_is_tgeom = false for next iteration
                              prev_is_tgeom = this_is_tgeom;
                          }
              }
        }
        
        if(i==seq->count-1 && instants_temp_count>0){
            if(inst1 == NULL){
                // printf("Debug #3: %d \n", instants_temp_count);
                // fflush(stdout);
                sequences_tg[count_tg] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                //tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                count_tg++;

                //free the temporary array of instants after making the sequence
                // pfree_array((void **) instants_temp, instants_temp_count);
                //initialize the temporary array of instants counter to 0
                instants_temp_count = 0;
            }else{
                sequences_tn[count_tn] = tsequence_make_c((const TInstant **)instants_temp, 
                                  instants_temp_count, seq->count,
                                  seq->period.lower_inc, seq->period.upper_inc,
                                  MOBDB_FLAGS_GET_LINEAR(seq->flags),true, i);
                // tsequence_make((const TInstant **)instants_temp, instants_temp_count,true,true,true,true);
                count_tn++;

                //free the temporary array of instants after making the sequence
                // pfree_array((void **) instants_temp, instants_temp_count);
                //initialize the temporary array of instants counter to 0
                instants_temp_count = 0;
            }
        }
    }    
  }
  // Datum result = tsequenceset_make_free(sequences, ss->count, true);
  //Datum result_combined=tgeompoint_tnpointHeapTupleDatum(fcinfo,(Temporal *)result_tnpoints,(Temporal *)result_tgeom);
  TSequenceSet *result_tnpoints = tsequenceset_make_free_c(sequences_tn, count_tn, MERGE_NO);
  TSequenceSet *result_tgeom = tsequenceset_make_free_c(sequences_tg, count_tg, MERGE_NO);

  Datum result_combined=tgeompoint_tnpointHeapTupleDatum(fcinfo,(Temporal *)result_tnpoints,(Temporal *)result_tgeom);
  return result_combined;
}

Datum 
tgeompoint_tnpointComplex(FunctionCallInfo fcinfo, const Temporal *temp)
{
    int32_t srid_tpoint = tpoint_srid(temp);
    int32_t srid_ways = get_srid_ways();
    ensure_same_srid(srid_tpoint, srid_ways);
    Datum result;
    ensure_valid_tempsubtype(temp->subtype);
    if (temp->subtype == TINSTANT)
        result = tgeompointinst_tnpointComplexinst(fcinfo, (TInstant *) temp);
    else if (temp->subtype == TINSTANTSET)
        result = tgeompointinstset_tnpointComplexinstset(fcinfo, (TInstantSet *) temp);
    else if (temp->subtype == TSEQUENCE)
        result = tgeompointseq_tnpointComplexseq(fcinfo, (TSequence *) temp);
    else /* temp->subtype == TSEQUENCESET */
     result = tgeompointseqset_tnpointComplexseqset(fcinfo, (TSequenceSet *) temp);
    
    return result;
}

// Datum 
// tnpointComplex_tgeompoint(const Temporal *temp)
// {
//         return result;
// }

/**
 * @brief Cast a temporal network complex point as a temporal geometric point.
 */
Temporal *
tnpointComplex_tgeompoint(const Temporal *temp)
{
  Temporal *result;
  ensure_valid_tempsubtype(temp->subtype);
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