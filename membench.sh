#!/bin/bash

TIME=/usr/bin/time

bench_log=/tmp/bench.tmp
csv=/tmp/csv1.tmp
csv_tmp=/tmp/csv2.tmp
rm -f $bench_log $csv $csv_tmp

for qftn in `seq 2 16`; do
    echo -n "$qftn " >> $bench_log;
    $TIME -f "%M" ./qvm qft/qft$qftn.mc > /dev/null 2>> $bench_log;
    #echo "" >> $bench_log;
  done

mv -f $bench_log benchmark_data.csv
