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
 * @brief Ever spatial relationships for temporal network points.
 *
 * These relationships compute the ever spatial relationship between the
 * arguments and return a Boolean. These functions may be used for filtering
 * purposes before applying the corresponding temporal spatial relationship.
 *
 * The following relationships are supported:
 * contains, disjoint, intersects, touches, and dwithin
 */

#include "npoint/tnpoint_spatialrels.h"

/* PostgreSQL */
#include <postgres.h>
#include <utils/palloc.h>
#include <fmgr.h>
/* MobilityDB */
#include "general/lifting.h"
#include "point/tpoint_spatialfuncs.h"
#include "point/tpoint_spatialrels.h"
#include "npoint/tnpoint_spatialfuncs.h"
/* Changes */
#include <meos_internal.h>
#include <utils/array.h>
#include "pg_general/temporal_catalog.h"
#include "pg_npoint/tnpoint_static.h"
#include "general/temporal_util.h"
#include <executor/spi.h>
#include "point/tpoint_out.h"
#include <funcapi.h>

/*****************************************************************************
 * Generic binary functions for tnpoint <rel> (geo | Npoint)
 *****************************************************************************/

/**
 * @brief Generic spatial relationships for a temporal network point and a geometry
 */
Datum
spatialrel_geo_tnpoint_ext(FunctionCallInfo fcinfo,
  Datum (*func)(Datum, Datum))
{
  Datum geom = PG_GETARG_DATUM(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  Datum result = spatialrel_tnpoint_geo(temp, geom, func, INVERT);
  PG_FREE_IF_COPY(temp, 1);
  PG_RETURN_DATUM(result);
}

/**
 * @brief Generic spatial relationships for a temporal network point and a geometry
 */
static Datum
spatialrel_tnpoint_geo_ext(FunctionCallInfo fcinfo,
  Datum (*func)(Datum, Datum))
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Datum geom = PG_GETARG_DATUM(1);
  Datum result = spatialrel_tnpoint_geo(temp, geom, func, INVERT_NO);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_DATUM(result);
}

/**
 * @brief Generic spatial relationships for a temporal network point and a network point
 */
static Datum
spatialrel_npoint_tnpoint_ext(FunctionCallInfo fcinfo,
  Datum (*func)(Datum, Datum))
{
  Npoint *np = PG_GETARG_NPOINT_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  Datum result = spatialrel_tnpoint_npoint(temp, np, func, INVERT);
  PG_FREE_IF_COPY(temp, 1);
  PG_RETURN_DATUM(result);
}

/**
 * @brief Generic spatial relationships for a temporal network point and a network point
 */
static Datum
spatialrel_tnpoint_npoint_ext(FunctionCallInfo fcinfo,
  Datum (*func)(Datum, Datum))
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np = PG_GETARG_NPOINT_P(1);
  Datum result = spatialrel_tnpoint_npoint(temp, np, func, INVERT_NO);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_DATUM(result);
}

/**
 * @brief Return true if the temporal network points ever satisfy the spatial
 * relationship
 *
 * @param[in] fcinfo Catalog information about the external function
 * @param[in] func Spatial relationship
 */
static Datum
spatialrel_tnpoint_tnpoint_ext(FunctionCallInfo fcinfo,
  Datum (*func)(Datum, Datum))
{
  Temporal *temp1 = PG_GETARG_TEMPORAL_P(0);
  Temporal *temp2 = PG_GETARG_TEMPORAL_P(1);
  int result = spatialrel_tnpoint_tnpoint(temp1, temp2, func);
  PG_FREE_IF_COPY(temp1, 0);
  PG_FREE_IF_COPY(temp2, 1);
  if (result < 0)
    PG_RETURN_NULL();
  PG_RETURN_BOOL(result == 1 ? true : false);
}

/*****************************************************************************
 * Ever contains
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Contains_geo_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the geometry contains the trajectory of the temporal network
 * point
 * @sqlfunc contains()
 */
PGDLLEXPORT Datum
Contains_geo_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_geo_tnpoint_ext(fcinfo, &geom_contains);
}

/*****************************************************************************
 * Ever disjoint
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Disjoint_geo_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the geometry and the trajectory of the temporal network
 * point are disjoint
 * @sqlfunc disjoint()
 */
PGDLLEXPORT Datum
Disjoint_geo_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_geo_tnpoint_ext(fcinfo, &geom_disjoint2d);
}

PG_FUNCTION_INFO_V1(Disjoint_npoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the network point and the trajectory of the temporal
 * network point are disjoint
 * @sqlfunc disjoint()
 */
PGDLLEXPORT Datum
Disjoint_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_npoint_tnpoint_ext(fcinfo, &geom_disjoint2d);
}

PG_FUNCTION_INFO_V1(Disjoint_tnpoint_geo);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * geometry are disjoint
 * @sqlfunc disjoint()
 */
PGDLLEXPORT Datum
Disjoint_tnpoint_geo(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_geo_ext(fcinfo, &geom_disjoint2d);
}

PG_FUNCTION_INFO_V1(Disjoint_tnpoint_npoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * network point are disjoint
 * @sqlfunc disjoint()
 */
PGDLLEXPORT Datum
Disjoint_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_npoint_ext(fcinfo, &geom_disjoint2d);
}

PG_FUNCTION_INFO_V1(Disjoint_tnpoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the temporal points are ever disjoint
 * @sqlfunc disjoint()
 */
PGDLLEXPORT Datum
Disjoint_tnpoint_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_tnpoint_ext(fcinfo, &datum2_point_ne);
}

/*****************************************************************************
 * Ever intersects
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Intersects_geo_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the geometry and the trajectory of the temporal network
 * point intersect
 * @sqlfunc intersects()
 */
PGDLLEXPORT Datum
Intersects_geo_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_geo_tnpoint_ext(fcinfo, &geom_intersects2d);
}

PG_FUNCTION_INFO_V1(Intersects_npoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the network point and the trajectory of the temporal network
 * point intersect
 * @sqlfunc intersects()
 */
PGDLLEXPORT Datum
Intersects_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_npoint_tnpoint_ext(fcinfo, &geom_intersects2d);
}

PG_FUNCTION_INFO_V1(Intersects_tnpoint_geo);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * geometry intersect
 * @sqlfunc intersects()
 */
PGDLLEXPORT Datum
Intersects_tnpoint_geo(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_geo_ext(fcinfo, &geom_intersects2d);
}

PG_FUNCTION_INFO_V1(Intersects_tnpoint_npoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the network
 * point intersect
 * @sqlfunc intersects()
 */
PGDLLEXPORT Datum
Intersects_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_npoint_ext(fcinfo, &geom_intersects2d);
}

PG_FUNCTION_INFO_V1(Intersects_tnpoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the temporal points are ever disjoint
 * @sqlfunc intersects()
 */
PGDLLEXPORT Datum
Intersects_tnpoint_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_tnpoint_ext(fcinfo, &datum2_point_eq);
}

/*****************************************************************************
 * Ever dwithin
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Dwithin_geo_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the geometry and the trajectory of the temporal network
 * point are within the given distance
 * @sqlfunc dwithin()
 */
PGDLLEXPORT Datum
Dwithin_geo_tnpoint(PG_FUNCTION_ARGS)
{
  Datum geom = PG_GETARG_DATUM(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  double dist = PG_GETARG_FLOAT8(2);
  Datum result = spatialrel3_tnpoint_geom(temp, geom, dist, &geom_dwithin2d,
    INVERT);
  PG_FREE_IF_COPY(temp, 1);
  PG_RETURN_DATUM(result);
}

PG_FUNCTION_INFO_V1(Dwithin_npoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the network point and the trajectory of the temporal network
 * point are within the given distance
 * @sqlfunc dwithin()
 */
PGDLLEXPORT Datum
Dwithin_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  Npoint *np = PG_GETARG_NPOINT_P(0);
  Temporal *temp = PG_GETARG_TEMPORAL_P(1);
  double dist = PG_GETARG_FLOAT8(2);
  Datum result = spatialrel3_tnpoint_npoint(temp, np, dist, &geom_dwithin2d,
    INVERT);
  PG_FREE_IF_COPY(temp, 1);
  PG_RETURN_DATUM(result);
}

PG_FUNCTION_INFO_V1(Dwithin_tnpoint_geo);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * geometry are within the given distance
 * @sqlfunc dwithin()
 */
PGDLLEXPORT Datum
Dwithin_tnpoint_geo(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Datum geom = PG_GETARG_DATUM(1);
  double dist = PG_GETARG_FLOAT8(2);
  Datum result = spatialrel3_tnpoint_geom(temp, geom, dist, &geom_dwithin2d,
    INVERT_NO);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_DATUM(result);
}

PG_FUNCTION_INFO_V1(Dwithin_tnpoint_npoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * network point are within the given distance
 * @sqlfunc dwithin()
 */
PGDLLEXPORT Datum
Dwithin_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np = PG_GETARG_NPOINT_P(1);
  double dist = PG_GETARG_FLOAT8(2);
  Datum result = spatialrel3_tnpoint_npoint(temp, np, dist, &geom_dwithin2d,
    INVERT_NO);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_DATUM(result);
}

PG_FUNCTION_INFO_V1(Dwithin_tnpoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectories of the temporal network points are within
 * the given distance
 * @sqlfunc dwithin()
 */
PGDLLEXPORT Datum
Dwithin_tnpoint_tnpoint(PG_FUNCTION_ARGS)
{
  Temporal *temp1 = PG_GETARG_TEMPORAL_P(0);
  Temporal *temp2 = PG_GETARG_TEMPORAL_P(1);
  double dist = PG_GETARG_FLOAT8(2);
  int result = dwithin_tnpoint_tnpoint(temp1, temp2, dist);
  PG_FREE_IF_COPY(temp1, 0);
  PG_FREE_IF_COPY(temp2, 1);
  if (result < 0)
    PG_RETURN_NULL();
  PG_RETURN_FLOAT8(result);
}

/*****************************************************************************
 * Ever touches
 *****************************************************************************/

PG_FUNCTION_INFO_V1(Touches_geo_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the geometry and the trajectory of the temporal network
 * point touch
 * @sqlfunc touches()
 */
PGDLLEXPORT Datum
Touches_geo_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_geo_tnpoint_ext(fcinfo, &geom_touches);
}

PG_FUNCTION_INFO_V1(Touches_npoint_tnpoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the network point and the trajectory of the temporal
 * network point touch
 * @sqlfunc touches()
 */
PGDLLEXPORT Datum
Touches_npoint_tnpoint(PG_FUNCTION_ARGS)
{
  return spatialrel_npoint_tnpoint_ext(fcinfo, &geom_touches);
}

PG_FUNCTION_INFO_V1(Touches_tnpoint_geo);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * geometry touch
 * @sqlfunc touches()
 */
PGDLLEXPORT Datum
Touches_tnpoint_geo(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_geo_ext(fcinfo, &geom_touches);
}

PG_FUNCTION_INFO_V1(Touches_tnpoint_npoint);
/**
 * @ingroup mobilitydb_temporal_spatial_rel
 * @brief Return true if the trajectory of the temporal network point and the
 * network point touch
 * @sqlfunc touches()
 */
PGDLLEXPORT Datum
Touches_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  return spatialrel_tnpoint_npoint_ext(fcinfo, &geom_touches);
}

/*****************************************************************************/
/*****************************************************************************/
/**
 * @ingroup libmeos_temporal_cast
 * @brief Cast a temporal network point sequence as a temporal network segment.
 */
Nsegment *
tnpointseq_nsegment(const TSequence *seq)
{
  Npoint *np = DatumGetNpointP(tinstant_value(tsequence_inst_n(seq, 0)));
  int64 rid = np->rid;
  double pos1=1;
  double pos2=0;
  //check whether the type is linear or stepwise interpolation?
  //check open and closed bracket
  for (int i = 0; i < seq->count; i++)
  {
    const TInstant *inst = tsequence_inst_n(seq, i);
    np = DatumGetNpointP(tinstant_value(inst));
    pos1 = pos1 < np->pos ? pos1: np->pos;
    pos2 = pos2 > np->pos ? pos2: np->pos;
  }
  Nsegment * result = nsegment_make(rid,pos1,pos2);
  //pfree??
  return result;
}

/**
 * @brief Generic spatial relationships for a temporal network point and a
 * network point
 */
Datum
spatialrel1_tnpoint_npoint(const Temporal *temp, const Npoint *np, const double precision)
{
  bool result;
  //double precision = precision || MOBDB_EPSILON;
  ensure_valid_tempsubtype(temp->subtype);
  if (temp->subtype == TINSTANT){
    Npoint *np2 = DatumGetNpointP(tinstant_value((TInstant *) temp));
    result =  np->rid == np2->rid && fabs(np->pos - np2->pos) < precision;
  }  
  else if (temp->subtype == TINSTANTSET){
    const TInstantSet *ti = (TInstantSet *) temp;
     for (int i = 0; i < ti->count; i++)
     {
       Npoint *np2 = DatumGetNpointP(tinstant_value(tinstantset_inst_n(ti, i)));
       result =  np->rid == np2->rid && fabs(np->pos - np2->pos) < precision;
       if(result)
        break;
     }
  }
  else if (temp->subtype == TSEQUENCE){
    //tnpointseq_step_npoints
    //pending for type=stepwise
    const TSequence *seq = (TSequence *) temp;
    Nsegment *ns = tnpointseq_nsegment(seq);
    //what if + - MOBDB_EPSILON is smaller or greater than 0, 1??
    result = np->rid == ns->rid && ((ns->pos1)-precision <= np->pos && np->pos <= (ns->pos2)+precision);
  }
  else /* temp->subtype == TSEQUENCESET */{
    const TSequenceSet *ts = (TSequenceSet *) temp;    
    //tnpointseqset_step_npoints
    for (int i = 0; i < ts->count; i++)
    {
       const TSequence *seq = tsequenceset_seq_n(ts, i);
       Nsegment *ns = tnpointseq_nsegment(seq);
       result = np->rid == ns->rid && ((ns->pos1)-precision <= np->pos && np->pos <= (ns->pos2)+precision);
       if(result)
        break;
    }
  }
  //PG_RETURN_BOOL(result);
  return result;
}

/**
 * @brief Generic spatial relationships for a temporal network point and a
 * network segment
 */
Datum
spatialrel1_tnpoint_nsegment(const Temporal *temp, const Nsegment *np, const double precision)
{
  bool result;
  //double precision = precision || MOBDB_EPSILON;
  ensure_valid_tempsubtype(temp->subtype);
  if (temp->subtype == TINSTANT){
    Npoint *np2 = DatumGetNpointP(tinstant_value((TInstant *) temp));
    result =  np->rid == np2->rid && (fabs(np->pos1 - np2->pos) < precision || fabs(np->pos2 - np2->pos) < precision);
  }
  else if (temp->subtype == TINSTANTSET){
    const TInstantSet *ti = (TInstantSet *) temp;
     for (int i = 0; i < ti->count; i++)
     {
       Npoint *np2 = DatumGetNpointP(tinstant_value(tinstantset_inst_n(ti, i)));
       result =  np->rid == np2->rid && (fabs(np->pos1 - np2->pos) < precision || fabs(np->pos2 - np2->pos) < precision);
       if(result)
        break;
     }
  }
  else if (temp->subtype == TSEQUENCE){
    //tnpointseq_step_npoints
    //pending for type=stepwise
    const TSequence *seq = (TSequence *) temp;
    Nsegment *ns = tnpointseq_nsegment(seq);
    //what if + - MOBDB_EPSILON is smaller or greater than 0, 1??
    //lower part of np is greater or equal to lower part of ns
    //upper part of np is smaller or equal to upper part of ns
    result = np->rid == ns->rid && ( ( (np->pos1-precision)<=ns->pos1 && ns->pos1<=(np->pos2+precision) ) || ( (np->pos1-precision)<=ns->pos2 && ns->pos2<=(np->pos2+precision) ) );
    // result = np->rid == ns->rid && ((ns->pos1)-precision <= np->pos && np->pos <= (ns->pos2)+precision);
  }
  else /* temp->subtype == TSEQUENCESET */{
    const TSequenceSet *ts = (TSequenceSet *) temp;    
    //tnpointseqset_step_npoints
    for (int i = 0; i < ts->count; i++)
    {
       const TSequence *seq = tsequenceset_seq_n(ts, i);
       Nsegment *ns = tnpointseq_nsegment(seq);
       // result = np->rid == ns->rid && ((ns->pos1)-precision <= np->pos && np->pos <= (ns->pos2)+precision);
      result = np->rid == ns->rid && ( ( (np->pos1-precision)<=ns->pos1 && ns->pos1<=(np->pos2+precision) ) || ( (np->pos1-precision)<=ns->pos2 && ns->pos2<=(np->pos2+precision) ) );
       if(result)
        break;
    }
  }
  //PG_RETURN_BOOL(result);
  return result;
}

/**
 * @brief Generic spatial relationships for a temporal network point and a
 * network point without translating them into geometry
 */
// spatialrel_tnpoint_npoint_ext(FunctionCallInfo fcinfo,
//   Datum (*func)(Datum, Datum))
Datum
spatialrel1_tnpoint_npoint_ext(FunctionCallInfo fcinfo)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np = PG_GETARG_NPOINT_P(1);
  bool result = spatialrel1_tnpoint_npoint(temp, np, MOBDB_EPSILON);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_BOOL(result);
}

/**
 * @brief Generic spatial relationships for a temporal network point and a
 * network point without translating them into geometry
 */
// spatialrel_tnpoint_npoint_ext(FunctionCallInfo fcinfo,
//   Datum (*func)(Datum, Datum))
Datum
spatialrel1_tnpoint_nsegment_ext(FunctionCallInfo fcinfo)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Nsegment *ns = PG_GETARG_NSEGMENT_P(1);
  bool result = spatialrel1_tnpoint_nsegment(temp, ns, MOBDB_EPSILON);
  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_BOOL(result);
}

/**
 * @brief Passes relationships for a temporal network point and a
 * network segment
 */
Nsegment **
passes1_tnpoint_nsegment_ext(const Temporal *temp, const Nsegment *ns1, const double precision, int *count)
{
  // Nsegment **result_ns = NULL;
  // if (temp->subtype == TINSTANT){
  //  **result_ns = palloc0(sizeof(Nsegment *) * 1);
  // }else{
  //  Nsegment **result_ns = palloc0(sizeof(Nsegment *) * temp->count);
  // }
  
  bool result;
  ensure_valid_tempsubtype(temp->subtype);
  if (temp->subtype == TINSTANT){
    Nsegment **result_ns = palloc0(sizeof(Nsegment *) * 1);
    Npoint *np2 = DatumGetNpointP(tinstant_value((TInstant *) temp));
      result_ns[0] = nsegment_make(np2->rid,np2->pos,np2->pos);
      *count = 1;
      return result_ns;
  }  
  else if (temp->subtype == TINSTANTSET){
    Nsegment **result_ns = palloc0(sizeof(Nsegment *) * ((TInstantSet *)temp)->count);
    const TInstantSet *ti = (TInstantSet *) temp;
    int count1 =0;
     for (int i = 0; i < ti->count; i++)
     {
       Npoint *np2 = DatumGetNpointP(tinstant_value(tinstantset_inst_n(ti, i)));
       result =  ns1->rid == np2->rid && (fabs(ns1->pos1 - np2->pos) < precision || fabs(ns1->pos2 - np2->pos) < precision);
       if(result){
         result_ns[count1++] = nsegment_make(np2->rid,np2->pos,np2->pos);
         *count = count1;
       }
     }
     return result_ns;
  }
  else if (temp->subtype == TSEQUENCE){
    //tnpointseq_step_npoints
    //pending for type=stepwise
    Nsegment **result_ns = palloc0(sizeof(Nsegment *) * ((TSequence *)temp)->count);
    const TSequence *seq = (TSequence *) temp;
    Nsegment *ns2 = tnpointseq_nsegment(seq);
    double lowPoint=1;
    double highPoint=0;
    if(ns1->rid==ns2->rid){
      lowPoint = ns1->pos1 > ns2->pos1 ? ns1->pos1:ns2->pos1;
      highPoint = ns1->pos2 < ns2->pos2 ? ns1->pos2:ns2->pos2;
      // printf("The lowpoint is: %f ; ns1->pos1 is: %f; ns2->pos1 is: %f \n", lowPoint, ns1->pos1, ns2->pos1);
      // printf("The highpoint is: %f ; ns1->pos2 is: %f; ns2->pos2 is: %f \n", highPoint, ns1->pos2, ns2->pos2);
      // fflush(stdout);
      result_ns[0] = nsegment_make(ns1->rid, lowPoint, highPoint);
      *count = 1;
    }
    return result_ns;
  }
  else /* temp->subtype == TSEQUENCESET */{
  Nsegment **result_ns = palloc0(sizeof(Nsegment *) * ((TSequenceSet *)temp)->count);
    const TSequenceSet *ts = (TSequenceSet *) temp;   
    double lowPoint=1;
    double highPoint=0;
    int count1=0;
    for (int i = 0; i < ts->count; i++)
    {
       const TSequence *seq = tsequenceset_seq_n(ts, i);
       Nsegment *ns2 = tnpointseq_nsegment(seq);
      result = ns1->rid == ns2->rid && 
      ( ( (ns1->pos1-precision)<=ns2->pos1 && ns2->pos1<=(ns1->pos2+precision) ) ||
        ( (ns1->pos1-precision)<=ns2->pos2 && ns2->pos2<=(ns1->pos2+precision) ) );
       if(result){
        lowPoint = ns1->pos1 > ns2->pos1 ? ns1->pos1:ns2->pos1;
        highPoint = ns1->pos2 < ns2->pos2 ? ns1->pos2:ns2->pos2;
        result_ns[count1++] = nsegment_make(ns1->rid, lowPoint, highPoint);
        *count = count1;
       }       
    }
    return result_ns;
  }
}

/*****************************************************************************/
PG_FUNCTION_INFO_V1(contains1_tnpoint_npoint);
/**
 * Return true if the trajectory of the temporal network point and the
 * network point are within the given distance
 */
PGDLLEXPORT Datum
contains1_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  return spatialrel1_tnpoint_npoint_ext(fcinfo);
}

/*****************************************************************************/
PG_FUNCTION_INFO_V1(contains1_tnpoint_nsegment);
/**
 * Return true if the trajectory of the temporal network point and the
 * network segment are within the given distance
 */
PGDLLEXPORT Datum
contains1_tnpoint_nsegment(PG_FUNCTION_ARGS)
{
  return spatialrel1_tnpoint_nsegment_ext(fcinfo);
}

/*****************************************************************************/
/**
 * @brief Convert a C array of network segment values into a PostgreSQL array
 */
ArrayType *
nsegmentarr1_array(Nsegment **nsegmentarr, int count)
{
  return construct_array((Datum *)nsegmentarr, count, type_oid(T_NSEGMENT),
    sizeof(Nsegment), false, 'd');
}
PG_FUNCTION_INFO_V1(passes1_tnpoint_nsegment);
/**
 * Return true if the trajectory of the temporal network point and the
 * network segment are within the given distance
 */
PGDLLEXPORT Datum
passes1_tnpoint_nsegment(PG_FUNCTION_ARGS)
{
  const Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Nsegment *ns = PG_GETARG_NSEGMENT_P(1);
  Nsegment **result;
  int count;
  if(spatialrel1_tnpoint_nsegment(temp, ns, MOBDB_EPSILON))
  {
    result = passes1_tnpoint_nsegment_ext(temp, ns, MOBDB_EPSILON, &count);
  }else{
    PG_RETURN_NULL();
  }
  ArrayType *result_final = nsegmentarr1_array(result, count);
  //pfree_array((void **) result_final, count);
  // PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_POINTER(result_final);
}

/*****************************************************************************/
PG_FUNCTION_INFO_V1(closest1_tnpoint_npoint);
/**
 * Return the closest point in the trajectory of the temporal network point and the
 * network point in same routeid
 */
PGDLLEXPORT Datum
closest1_tnpoint_npoint(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);
  Npoint *np = PG_GETARG_NPOINT_P(1);
  PG_FREE_IF_COPY(temp, 0);
  Npoint *result;

  ensure_valid_tempsubtype(temp->subtype);
  if(spatialrel1_tnpoint_npoint(temp,np,0))
    result = np;
  else {
    if (temp->subtype == TINSTANT){
      Npoint *np2 = DatumGetNpointP(tinstant_value((TInstant *) temp));
      if(np->rid == np2->rid)
        result = np2;
      else
        PG_RETURN_NULL();
    }  
    else if (temp->subtype == TINSTANTSET){
      const TInstantSet *ti = (TInstantSet *) temp;
      double difference=1.01;//max distance between two point in a route is 1 or positive_infinity
       for (int i = 0; i < ti->count; i++)
       {
         Npoint *np2 = DatumGetNpointP(tinstant_value(tinstantset_inst_n(ti, i)));
         double temp_ = fabs(np2->pos-np->pos);
         if(np->rid == np2->rid && temp_ < difference){
            difference = temp_;
            result = np2; //check
         }else
          PG_RETURN_NULL();
       }
    }
    else if (temp->subtype == TSEQUENCE){
      const TSequence *seq = (TSequence *) temp;
      Nsegment *ns = tnpointseq_nsegment(seq);
      if(np->rid==ns->rid){
        double pos1 = nsegment_start_position(ns);
        double pos2 = nsegment_end_position(ns);
        if (fabs(pos1-np->pos) < fabs(pos2-np->pos))
          result=npoint_make(np->rid, pos1);
        else
          result=npoint_make(np->rid, pos2);
      }else
        PG_RETURN_NULL();
    }
    else /* temp->subtype == TSEQUENCESET */{
      const TSequenceSet *ts = (TSequenceSet *) temp;
      double difference=1.01; //considering a number greater than 1; as the maximum difference is 1 in tnpoints' start and end position for a route
      bool matchedflag = false;
      //pending for type=stepwise
      for (int i = 0; i < ts->count; i++)
      {
         const TSequence *seq = tsequenceset_seq_n(ts, i);
         Nsegment *ns = tnpointseq_nsegment(seq);
         
         if(np->rid==ns->rid){
          matchedflag = true;
          double pos1 = nsegment_start_position(ns);
          double pos2 = nsegment_end_position(ns);
          double temp1 = fabs(pos1-np->pos);
          double temp2 = fabs(pos2-np->pos);
          if(temp1<temp2){
            if(temp1<difference){
             result=npoint_make(np->rid, pos1);
             difference=temp1;
            }
          }
          else{
            if(temp2<difference){
             result=npoint_make(np->rid, pos2);
             difference=temp2;
            }
          }        
        }
      }
      if(!matchedflag)
          PG_RETURN_NULL();
    }
  }
  PG_RETURN_NPOINT_P(result);
}

/*****************************************************************************/
/**
 * @brief Find the Closest Point in different routes near a point.
 */
// Npoint *
// find_closest_point(const GSERIALIZED *gs)
// {
// //   /* Ensure validity of operation */
// //   ensure_non_empty(gs);
// //   ensure_point_type(gs);
// //   int32_t srid_geom = gserialized_get_srid(gs);
// //   int32_t srid_ways = get_srid_ways();
// //   ensure_same_srid(srid_geom, srid_ways);

// //   char *geomstr = ewkt_out(0, PointerGetDatum(gs), OUT_DEFAULT_DECIMAL_DIGITS);
// //   char sql[512];
// //   // sprintf(sql, "SELECT npoint(gid, ST_LineLocatePoint(the_geom, '%s')) "
// //   //   "FROM public.ways WHERE ST_DWithin(the_geom, '%s', %lf) "
// //   //   "ORDER BY ST_Distance(the_geom, '%s') LIMIT 1", geomstr, geomstr,
// //   //   DIST_EPSILON, geomstr);
// //   sprintf(sql, "SELECT npoint(gid, ST_LineLocatePoint(the_geom, '%s')) "
// //   "FROM public.ways WHERE ST_DWithin(the_geom, '%s', %lf) "
// //   "ORDER BY ST_Distance(the_geom, '%s') LIMIT 1", geomstr, geomstr,
// //   DIST_EPSILON, geomstr);

// // //   SELECT gid, 
// // // ST_ClosestPoint(the_geom,'SRID=4326;Point(4.3041517 50.7926352)'),
// // // npoint(gid, ST_LineLocatePoint(the_geom, 'SRID=4326;Point(4.3041517 50.7926352)')),
// // // ST_Distance(the_geom,'SRID=4326;Point(4.3041517 50.7926352)') AS distanceToPoint 
// // // FROM public.ways 
// // // WHERE ST_Intersects(st_buffer('SRID=4326;Point(4.3041517 50.7926352)'::geography, 500)::geometry, the_geom) 
// // // ORDER BY ST_Distance(the_geom,'SRID=4326;Point(4.3041517 50.7926352)') 
// // // LIMIT 100;

// //   pfree(geomstr);
// //   Npoint *result = palloc(sizeof(Npoint));
// //   bool isNull = true;
// //   SPI_connect();
// //   int ret = SPI_execute(sql, true, 1);
// //   uint64 proc = SPI_processed;
// //   if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
// //   {
// //     SPITupleTable *tuptable = SPI_tuptable;
// //     Datum value = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isNull);
// //     if (!isNull)
// //     {
// //       /* Must allocate this in upper executor context to keep it alive after SPI_finish() */
// //       Npoint *np = DatumGetNpointP(value);
// //       memcpy(result, np, sizeof(Npoint));
// //     }
// //   }
// //   SPI_finish();
// //   if (isNull)
// //   {
// //     pfree(result);
// //     return NULL;
// //   }
// //   return result;
//   return NULL;
// }
/**
 * @brief Find the length of the npoint.
 */
double
find_ways_length(const int64 rid, bool inMeters)
{
  char sql[512];
  double result;
    if(inMeters)
      sprintf(sql, "SELECT length_m FROM public.ways WHERE gid = %ld", rid);
    else
      sprintf(sql, "SELECT length FROM public.ways WHERE gid = %ld", rid);

    bool isNull = true;
    SPI_connect();

    int ret = SPI_execute(sql, true, 1);
    uint64 proc = SPI_processed;

    if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
    {
      SPITupleTable *tuptable = SPI_tuptable;
      Datum value = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isNull);
      if (!isNull)
      {
        /* Must allocate this in upper executor context to keep it alive after SPI_finish() */
        // Datum *result_0 = (Datum *)DatumGetPointer(value);//(* ((int64 *) DatumGetPointer(X)))
        result = DatumGetFloat8(value);//(* ((int64 *) DatumGetPointer(X)))
        //Datum *result_0 = value;
        //memcpy(result, result_0, sizeof(double));
      }
    }

    SPI_finish();
    if (isNull)
    {
      // pfree(result);
      // return NULL;
      elog(ERROR, "Cannot get the length for route %ld", rid);

    }
    return result;
}

/*****************************************************************************/
PG_FUNCTION_INFO_V1(find_npoint_distance);
/**
 * Return the distances and related attributes between tempral network point trajectory and the
 * geometry
 */
PGDLLEXPORT Datum
find_npoint_distance(PG_FUNCTION_ARGS)
{
  Npoint *np1 = PG_GETARG_NPOINT_P(0);
  Npoint *np2 = PG_GETARG_NPOINT_P(1);
  bool inMetersFlag = PG_GETARG_BOOL(2);

  if((int64)np1->rid==(int64)np2->rid){
    // PG_RETURN_DATUM(PointerGetDatum(find_ways_length(np1->rid,inMetersFlag)));
    // PG_RETURN_INT64( route_length(np1->rid) );
    // PG_RETURN_DATUM( route_length(np1->rid) );
    double fraction = fabs(np2->pos - np1->pos);
    PG_RETURN_DATUM( fraction * find_ways_length(np1->rid, inMetersFlag) );
  }else
  // PG_FREE_IF_COPY(temp, 0);
    PG_RETURN_NULL();
}

/*****************************************************************************/
#define OUT_DEFAULT_DECIMAL_DIGITS 15
/**
 * @brief Get the nearby roads (road_id).
 */
int64*
get_nearby_roads_ext(const GSERIALIZED *gs, const int limit)
{
  /* Ensure validity of operation */
  ensure_non_empty(gs);
  ensure_point_type(gs);
  int32_t srid_geom = gserialized_get_srid(gs);
  int32_t srid_ways = get_srid_ways();
  ensure_same_srid(srid_geom, srid_ways);

  char *geomstr = ewkt_out(0, PointerGetDatum(gs), OUT_DEFAULT_DECIMAL_DIGITS);
  char sql[512];
  sprintf(sql, "SELECT gid FROM public.ways ORDER BY ST_Distance(the_geom,'%s') LIMIT %d", geomstr, limit);
  
  pfree(geomstr);
  bool isNull = true;
  SPI_connect();

  int ret = SPI_execute(sql, true, limit);
  uint64 proc = SPI_processed;

  int64 *result;
  if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
  {
    SPITupleTable *tuptable = SPI_tuptable;
    for(uint64 i=0;i<tuptable->numvals;i++){
      Datum value = SPI_getbinval(tuptable->vals[i], tuptable->tupdesc, 1, &isNull);
      if (!isNull)
      {
        result[i] = DatumGetInt64(value);
      }
    }  
  }

  SPI_finish();
  if (isNull)
  {
    elog(ERROR, "Cannot get the nearest routeid");
  }
  return result;
  // return 108;
}
char*
get_nearby_roads_ext1(const GSERIALIZED *gs, const int limit)
{
  /* Ensure validity of operation */
  ensure_non_empty(gs);
  ensure_point_type(gs);
  int32_t srid_geom = gserialized_get_srid(gs);
  int32_t srid_ways = get_srid_ways();
  ensure_same_srid(srid_geom, srid_ways);

  char *geomstr = ewkt_out(0, PointerGetDatum(gs), OUT_DEFAULT_DECIMAL_DIGITS);
  char sql[512];
  sprintf(sql, "SELECT gid FROM public.ways ORDER BY ST_Distance(the_geom,'%s') LIMIT %d", geomstr, limit);
  
  pfree(geomstr);
  bool isNull = true;
  SPI_connect();

  int ret = SPI_execute(sql, true, 1);
  uint64 proc = SPI_processed;

  // int64 *result;
  // if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
  // {
  //   SPITupleTable *tuptable = SPI_tuptable;
  //   Datum value = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isNull);
  //   if (!isNull)
  //   {
  //     result = DatumGetInt64(value);
  //   }
  // }

  char *result;
  if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
  {
    SPITupleTable *tuptable = SPI_tuptable;
    result = SPI_getvalue(tuptable->vals[0], tuptable->tupdesc, 1);
    // if (!isNull)
    // {
    //   result = DatumGetCString(value);
    // }
  }

  SPI_finish();
  // if (isNull)
  // {
  //   elog(ERROR, "Cannot get the nearest routeid");
  // }
  return result;


  // return 108;
}

Datum
get_nearby_roads_ext_2(FunctionCallInfo fcinfo)
// , const GSERIALIZED *gs, const int limit)
{
  GSERIALIZED *gs =((GSERIALIZED *)PG_DETOAST_DATUM(PG_GETARG_DATUM(0)));
  int limit = PG_GETARG_INT32(1);
  FuncCallContext *funcctx;
  int call_cntr;
  int max_calls;
  int64 *result_in=palloc(sizeof(int64)*limit);
  int64 result;
  // bool isNull[1]={0};

  // Datum tupple_array[1]; Used to Construct the Composite Return Value
  // HeapTuple tuple;
  // Datum result; /*Returned Composite Value*/

    /* Ensure validity of operation */
    ensure_non_empty(gs);
    ensure_point_type(gs);
    int32_t srid_geom = gserialized_get_srid(gs);
    int32_t srid_ways = get_srid_ways();
    ensure_same_srid(srid_geom, srid_ways);

    char *geomstr = ewkt_out(0, PointerGetDatum(gs), OUT_DEFAULT_DECIMAL_DIGITS);
    char sql[512];
    sprintf(sql, "SELECT gid FROM public.ways ORDER BY ST_Distance(the_geom,'%s') LIMIT %d", geomstr, limit);
    
    pfree(geomstr);
    bool isNull = true;
    SPI_connect();

    int ret = SPI_execute(sql, true, limit);
    uint64 proc = SPI_processed;

    // int64 result;
    if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
    {

      SPITupleTable *tuptable = SPI_tuptable;
      for(uint64 i=0;i<tuptable->numvals;i++){
        Datum value = SPI_getbinval(tuptable->vals[i], tuptable->tupdesc, 1, &isNull);
        if (!isNull)
        {
          result_in[i] = DatumGetInt64(value);
        }
      }

    }

    SPI_finish();
    if (isNull)
    {
      elog(ERROR, "Cannot get the nearest routeid");
    }
    
    /*for first function call*/
    if(SRF_IS_FIRSTCALL()){

      /*Initialize the FuncCallContext */
      funcctx = SRF_FIRSTCALL_INIT();

      /*Switch to memory context appropriate for multiple function calls */
      MemoryContext oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

      /* total number of tuples to be returned */
      funcctx->max_calls = limit;

      MemoryContextSwitchTo(oldcontext);
    }

     /* stuff done on every call of the function */
    funcctx = SRF_PERCALL_SETUP();

    call_cntr = funcctx->call_cntr;
    max_calls = funcctx->max_calls;

    if (call_cntr < max_calls)    /* do when there is more left to send */
    {
      result = result_in[call_cntr];
      SRF_RETURN_NEXT(funcctx, result);
    }else{
      SRF_RETURN_DONE(funcctx);
    }
    
}

// CREATE TYPE nearbyTuple AS (
//   road_id bigint
//   shortest_distance float,
//   closest_point geometry,
//   closest_npoint npoint
// );

Datum
get_nearby_roads_ext_3(FunctionCallInfo fcinfo)
// , const GSERIALIZED *gs, const int limit)
{
  GSERIALIZED *gs =((GSERIALIZED *)PG_DETOAST_DATUM(PG_GETARG_DATUM(0)));
  int limit = PG_GETARG_INT32(1);
  FuncCallContext *funcctx;
  int call_cntr;
  int max_calls;
  int64 *result_in_gid=palloc(sizeof(int64)*limit);
  float *result_in_distance=palloc(sizeof(float)*limit);
  Npoint **result_in_npoint=palloc(sizeof(Npoint)*limit);
  bool *result_in_bool=palloc(sizeof(bool)*limit);
  // int64 result;

  bool result_is_null[3]={0,0,0};
  Datum result_values[3]; //Used to Construct the Composite Return Value
  HeapTuple tuple;
  Datum result_combined; /*Returned Composite Value*/

    /* Ensure validity of operation */
    ensure_non_empty(gs);
    ensure_point_type(gs);
    int32_t srid_geom = gserialized_get_srid(gs);
    int32_t srid_ways = get_srid_ways();
    ensure_same_srid(srid_geom, srid_ways);

    char *geomstr = ewkt_out(0, PointerGetDatum(gs), OUT_DEFAULT_DECIMAL_DIGITS);
    char sql[512];
    sprintf(sql, "SELECT gid, ST_Distance(the_geom,'%s'), npoint(gid, ST_LineLocatePoint(the_geom, ST_ClosestPoint(the_geom,'%s'))) FROM public.ways ORDER BY ST_Distance(the_geom,'%s') LIMIT %d", geomstr, geomstr, geomstr, limit);
    
    pfree(geomstr);
    bool isNull = true;
    SPI_connect();

    int ret = SPI_execute(sql, true, limit);
    uint64 proc = SPI_processed;

    // int64 result;
    if (ret > 0 && proc > 0 && SPI_tuptable != NULL)
    {

      SPITupleTable *tuptable = SPI_tuptable;
      for(uint64 i=0;i<tuptable->numvals;i++){
        Datum value = SPI_getbinval(tuptable->vals[i], tuptable->tupdesc, 1, &isNull);
        result_in_bool[i] = isNull;
        if (!isNull)
        {
          result_in_gid[i] = DatumGetInt64(value);
          result_in_distance[i] = DatumGetFloat8(SPI_getbinval(tuptable->vals[i], tuptable->tupdesc, 2, &isNull));
          result_in_npoint[i] = DatumGetNpointP(SPI_getbinval(tuptable->vals[i], tuptable->tupdesc, 3, &isNull));
        }
      }

    }

    SPI_finish();
    if (isNull)
    {
      elog(ERROR, "Cannot get the nearest routeid");
    }
    
    /*for first function call*/
    if(SRF_IS_FIRSTCALL()){

      /*Initialize the FuncCallContext */
      funcctx = SRF_FIRSTCALL_INIT();

      /*Switch to memory context appropriate for multiple function calls */
      MemoryContext oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

      /* total number of tuples to be returned */
      funcctx->max_calls = limit;

      MemoryContextSwitchTo(oldcontext);
    }

     /* stuff done on every call of the function */
    funcctx = SRF_PERCALL_SETUP();

    call_cntr = funcctx->call_cntr;
    max_calls = funcctx->max_calls;

    if (call_cntr < max_calls)    /* do when there is more left to send */
    {
      result_values[0]=result_in_gid[call_cntr];
      result_values[1]=result_in_distance[call_cntr];
      result_values[2]=PointerGetDatum(result_in_npoint[call_cntr]);
      if(!result_in_bool[call_cntr]){
        // result_is_null[3]= {1,1,1};
      }
      /* Form tuple and return */
      tuple = heap_form_tuple(funcctx->tuple_desc, result_values, result_is_null);
      result_combined = HeapTupleGetDatum(tuple);
      // result = result_in_gid[call_cntr];
      SRF_RETURN_NEXT(funcctx, result_combined);
    }else{
      SRF_RETURN_DONE(funcctx);
    }    
}

PG_FUNCTION_INFO_V1(get_nearby_roads);
/**
 * Return the distances and related attributes between tempral network point trajectory and the
 * geometry
 */
PGDLLEXPORT Datum
get_nearby_roads(PG_FUNCTION_ARGS)
{
  // GSERIALIZED *gs = PG_GETARG_GSERIALIZED_P(0);
  // GSERIALIZED *gs =((GSERIALIZED *)PG_DETOAST_DATUM(PG_GETARG_DATUM(0)));
  // int limit = PG_GETARG_INT32(1);

  // int64 *result;
  // result = get_nearby_roads_ext(gs,limit);
  // //PG_RETURN_NULL();
  // PG_RETURN_INT64(result);
  // PG_RETURN_DATUM(result);

  // char *result;
  // result = get_nearby_roads_ext(gs,limit);
  // //PG_RETURN_NULL();
  // PG_RETURN_CSTRING(result);
  return get_nearby_roads_ext_2(fcinfo);
}


/**
 * @brief Find the distance between two temporal network points.
 */
double *
find_npoint_distance_1(const Npoint *pointa, const Npoint *pointb, bool inMeters)
{
  double result=0;
  if(pointa->rid==pointb->rid){
    // double *length = find_ways_length(pointa->rid,inMeters);
    // if(length){
    if(1){
      //return DatumGetPointer( ((double)PointerGetDatum(length)) * fabs(pointa->pos-pointb->pos) );
      // return find_ways_length(pointa->rid,inMeters);
      return NULL;
    }else
      return NULL;
  }else
    return NULL;
}

/*****************************************************************************/
PG_FUNCTION_INFO_V1(distances1_tnpoint_geometry);
/**
 * Return the distances and related attributes between tempral network point trajectory and the
 * geometry
 */
PGDLLEXPORT Datum
distances1_tnpoint_geometry(PG_FUNCTION_ARGS)
{
  Temporal *temp = PG_GETARG_TEMPORAL_P(0);  
  // Point     *center = PG_GETARG_POINT_P(1);

  // GSERIALIZED *gs = PG_GETARG_GSERIALIZED_P(1);


  // double    radius = PG_GETARG_FLOAT8(1);

  //get the nearest rid

  //query the nearest rid's to get the closest point

  //check if the closest point lies in the route id

  //return closest point

  //calculate the geometric distance between the given point and the closest point

  //get the distance between the closest point and closest point in the trajectory

  PG_FREE_IF_COPY(temp, 0);
  PG_RETURN_NULL();
}