#!/bin/bash

qpts=(10000 25000 50000 75000 100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000 25000000 50000000 75000000 100000000)

for q in "${qpts[@]}"; do
    ./random-data.py \
        -q "$q" \
        -o "datafiles/data-$q.csv"
done
