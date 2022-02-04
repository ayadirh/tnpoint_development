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
 * @file time_gist.c
 * R-tree GiST index for time types.
 */

#ifndef __TIME_GIST_H__
#define __TIME_GIST_H__

#include <postgres.h>
#include <catalog/pg_type.h>
#include "timetypes.h"

/*****************************************************************************/

extern Datum period_gist_union(PG_FUNCTION_ARGS);
extern Datum timestampset_gist_compress(PG_FUNCTION_ARGS);
extern Datum period_gist_compress(PG_FUNCTION_ARGS);
extern Datum periodset_gist_compress(PG_FUNCTION_ARGS);
extern Datum period_gist_penalty(PG_FUNCTION_ARGS);
extern Datum period_gist_picksplit(PG_FUNCTION_ARGS);
extern Datum period_gist_same(PG_FUNCTION_ARGS);
extern Datum period_gist_fetch(PG_FUNCTION_ARGS);

extern int common_entry_cmp(const void *i1, const void *i2);

extern bool period_index_consistent_leaf(const Period *key, const Period *query,
  StrategyNumber strategy);
extern bool period_gist_consistent_internal(const Period *key, const Period *query,
  StrategyNumber strategy);
extern bool period_index_recheck(StrategyNumber strategy);

#endif

/*****************************************************************************/
