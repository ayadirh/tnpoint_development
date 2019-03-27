﻿-------------------------------------------------------------------------------
-- Tests of operators that do not involved indexes for time types.
-- File TimeOps.c
-------------------------------------------------------------------------------

SELECT timestamptz '2000-01-01' -|- period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' -|- period '(2000-01-01, 2000-01-03]';

SELECT timestamptz '2000-01-01' -|- periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-01' -|- periodset '{(2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' -|- period '[2000-01-01, 2000-01-03]';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' -|- periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT period '[2000-01-01, 2000-01-03]' -|- timestamptz '2000-01-01';
SELECT period '[2000-01-01, 2000-01-03]' -|- timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03]' -|- period '[2000-01-01, 2000-01-03]';
SELECT period '[2000-01-01, 2000-01-03]' -|- periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' -|- timestamptz '2000-01-01';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' -|- timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' -|- period '[2000-01-01, 2000-01-03]';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' -|- periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

-------------------------------------------------------------------------------

SELECT timestamptz '2000-01-01' + timestamptz '2000-01-01';
SELECT timestamptz '2000-01-01' + timestamptz '2000-01-02';
SELECT timestamptz '2000-01-01' + timestampset '{2000-01-02, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-01' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-05' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-06' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-01' + period '[2000-01-02, 2000-01-03]';
SELECT timestamptz '2000-01-01' + period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' + period '(2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-02' + period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-03' + period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-03' + period '[2000-01-01, 2000-01-03)';
SELECT timestamptz '2000-01-05' + period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' + periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-01' + periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-03' + periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-04' + periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-04' + periodset '{[2000-01-02, 2000-01-03],[2000-01-05, 2000-01-05]}';
SELECT timestamptz '2000-01-05' + periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-06' + periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' + timestamptz '2000-01-01';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' + timestampset '{2000-01-03, 2000-01-05, 2000-01-07}';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' + period '[2000-01-01, 2000-01-03]';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' + periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT period '[2000-01-01, 2000-01-03]' + timestamptz '2000-01-01';
SELECT period '[2000-01-01, 2000-01-03]' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '(2000-01-01, 2000-01-03]' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03)' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03]' + period '[2000-01-01, 2000-01-03]';
SELECT period '[2000-01-01, 2000-01-03]' + period '(2000-01-03, 2000-01-05]';
SELECT period '[2000-01-01, 2000-01-03]' + periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' + timestamptz '2000-01-01';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' + timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' + period '[2000-01-01, 2000-01-03]';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' + periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

-------------------------------------------------------------------------------

SELECT timestamptz '2000-01-01' - timestamp '2000-01-01';
SELECT timestamptz '2000-01-01' - timestamp '2000-01-02';
SELECT timestamptz '2000-01-01' - timestampset '{2000-01-02, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-01' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-05' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-06' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-01' - period '[2000-01-02, 2000-01-03]';
SELECT timestamptz '2000-01-01' - period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' - period '(2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-02' - period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-03' - period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-03' - period '[2000-01-01, 2000-01-03)';
SELECT timestamptz '2000-01-05' - period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' - periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-01' - periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-03' - periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-04' - periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-04' - periodset '{[2000-01-02, 2000-01-03],[2000-01-05, 2000-01-05]}';
SELECT timestamptz '2000-01-05' - periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-06' - periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' - timestamptz '2000-01-01';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' - timestampset '{2000-01-03, 2000-01-05, 2000-01-07}';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' - period '[2000-01-01, 2000-01-03]';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' - periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT period '[2000-01-01, 2000-01-03]' - timestamptz '2000-01-01';
SELECT period '[2000-01-01, 2000-01-03]' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '(2000-01-01, 2000-01-03]' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03)' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03]' - period '[2000-01-01, 2000-01-03]';
SELECT period '[2000-01-01, 2000-01-03]' - period '(2000-01-03, 2000-01-05]';
SELECT period '[2000-01-01, 2000-01-03]' - periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' - timestamptz '2000-01-01';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' - timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' - period '[2000-01-01, 2000-01-03]';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' - periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

-------------------------------------------------------------------------------

SELECT timestamptz '2000-01-01' * timestamp '2000-01-01';
SELECT timestamptz '2000-01-01' * timestamp '2000-01-02';
SELECT timestamptz '2000-01-01' * timestampset '{2000-01-02, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-01' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-05' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-06' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestamptz '2000-01-01' * period '[2000-01-02, 2000-01-03]';
SELECT timestamptz '2000-01-01' * period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' * period '(2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-02' * period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-03' * period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-03' * period '[2000-01-01, 2000-01-03)';
SELECT timestamptz '2000-01-05' * period '[2000-01-01, 2000-01-03]';
SELECT timestamptz '2000-01-01' * periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-01' * periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-03' * periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-04' * periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-04' * periodset '{[2000-01-02, 2000-01-03],[2000-01-05, 2000-01-05]}';
SELECT timestamptz '2000-01-05' * periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';
SELECT timestamptz '2000-01-06' * periodset '{[2000-01-02, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' * timestamptz '2000-01-01';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' * timestampset '{2000-01-03, 2000-01-05, 2000-01-07}';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' * period '[2000-01-01, 2000-01-03]';
SELECT timestampset '{2000-01-01, 2000-01-03, 2000-01-05}' * periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT period '[2000-01-01, 2000-01-03]' * timestamptz '2000-01-01';
SELECT period '[2000-01-01, 2000-01-03]' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '(2000-01-01, 2000-01-03]' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03)' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT period '[2000-01-01, 2000-01-03]' * period '[2000-01-01, 2000-01-03]';
SELECT period '[2000-01-01, 2000-01-03]' * period '(2000-01-03, 2000-01-05]';
SELECT period '[2000-01-01, 2000-01-03]' * periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' * timestamptz '2000-01-01';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' * timestampset '{2000-01-01, 2000-01-03, 2000-01-05}';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' * period '[2000-01-01, 2000-01-03]';
SELECT periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}' * periodset '{[2000-01-01, 2000-01-03],[2000-01-04, 2000-01-05]}';

-------------------------------------------------------------------------------