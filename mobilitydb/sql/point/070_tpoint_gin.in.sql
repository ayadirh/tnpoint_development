/*
 * tpoint_gin.c
 * GIN index for temporal points.
 */
 
  -- CREATE FUNCTION comparePartial
 CREATE FUNCTION tnpoint_gin_compare(internal, internal)
  RETURNS int
  AS 'MODULE_PATHNAME', 'tnpoint_gin_compare'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

 CREATE FUNCTION tnpoint_gin_extract_value(internal,internal,internal)
  RETURNS internal
  AS 'MODULE_PATHNAME', 'tnpoint_gin_extract_value'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

 CREATE FUNCTION tnpoint_gin_extract_query(tnpoint, internal, smallint, internal, internal, internal, internal)
  RETURNS internal
  AS 'MODULE_PATHNAME', 'tnpoint_gin_extract_query'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

 -- CREATE FUNCTION triConsistent
 CREATE FUNCTION tnpoint_gin_consistent(internal, smallint, tnpoint, int32, internal, internal, internal, internal)
  RETURNS bool
  AS 'MODULE_PATHNAME', 'tnpoint_gin_consistent'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

 -- VOID OPTIONS
 CREATE OPERATOR CLASS tnpoint_gin_ops
  DEFAULT FOR TYPE tnpoint USING gin AS
  -- STORAGE stbox,
  -- -- overlaps
  OPERATOR  1    && (tnpoint, tnpoint),
  OPERATOR  1    && (nsegment, nsegment),
  OPERATOR  1    && (npoint, npoint),
  -- -- contains
  -- OPERATOR  2    @> (tnpoint, tnpoint),
  -- OPERATOR  2    @> (tnpoint, nsegment),
  -- OPERATOR  2    @> (tnpoint, npoint),  
  -- OPERATOR  2    @> (nsegment, nsegment),
  -- OPERATOR  2    @> (nsegment, npoint),
  -- OPERATOR  2    @> (npoint, npoint),
  OPERATOR  2    @> (tnpoint, int), 
  OPERATOR  2    @> (nsegment, int),
  OPERATOR  2    @> (npoint, int),
  -- -- is contained by
  -- OPERATOR  3    <@ (tnpoint, tnpoint),
  -- OPERATOR  3    <@ (nsegment, tnpoint),
  -- OPERATOR  3    <@ (npoint, tnpoint),
  -- OPERATOR  3    <@ (nsegment, nsegment),
  -- OPERATOR  3    <@ (npoint, nsegment),
  -- OPERATOR  3    <@ (npoint, npoint),
  --
  OPERATOR  3    <@ (int, tnpoint),
  OPERATOR  3    <@ (int, nsegment),
  OPERATOR  3    <@ (int, npoint),
  -- -- equal
  OPERATOR  4    = (npoint, npoint),
  OPERATOR  4    = (nsegment, nsegment),
  OPERATOR  4    = (tnpoint, tnpoint),
  -- function
  FUNCTION  1  tnpoint_gin_compare(internal,internal),
  FUNCTION  2  tnpoint_gin_extract_value(internal,internal,internal),
  FUNCTION  3  tnpoint_gin_extract_query(tnpoint, internal, smallint, internal, internal, internal, internal),
  FUNCTION  4  tnpoint_gin_consistent(internal, smallint, tnpoint, int32, internal, internal, internal, internal)
  -- FUNCTION  5  tnpoint_gin_compare_partial(internal,internal),
  -- FUNCTION  6  tnpoint_gin_tri_consistent(),
  -- FUNCTION  7  tnpoint_gin_options()
  ;
  -- INBUILT GIN OPERATOR CLASSES
    -- && (anyarray,anyarray)   //does the two values overlap
    -- @> (anyarray,anyarray)   //does the first value contain the second one
    -- <@ (anyarray,anyarray) //is the first value contained by the second
    -- = (anyarray,anyarray)    //are the two values equal

    --GIN ARRAY STRATEGIES
    -- OPERATION | STRATEGY NUMBER
    -- OVERLAP | 1
    -- CONTAINS | 2
    -- IS CONTAINED BY | 3
    -- EQUAL | 4
