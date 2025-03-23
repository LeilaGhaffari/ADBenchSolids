#!/bin/bash

qpts=(1000 10000 25000 50000 75000 100000 250000 500000 750000 1000000 2500000 5000000 7500000 10000000 25000000 50000000 75000000 100000000)

for arg in "$@"; do
    if [[ $arg == --models=* ]]; then
        models="${arg#*=}"
        shift
    elif [[ $arg == --output_file=* ]]; then
        file_name="${arg#*=}"
        shift
    else
        echo "Unknown argument: $arg"
        echo "Usage: $0 --models=model1,model2 --output_file=filename"
        exit 1
    fi
done

if [ -z "$models" ] || [ -z "$file_name" ]; then
    echo "Usage: $0 --models=model1,model2 --output_file=filename"
    exit 1
fi

echo "AD Benchmark: $models" > "$file_name"
echo -e "\n" >> "$file_name"

for iter in {1..3}; do
    echo "Iteration $iter:" >> "$file_name"
    for q in "${qpts[@]}"; do
        /home/leila/ADBenchSolids/build/elasticity-exec \
            -models="$models" \
            -data="/home/leila/ADBenchSolids/runs/benchmark/datafiles/data-$q.csv" >> "$file_name"
    done
    echo -e "\n" >> "$file_name"
done
