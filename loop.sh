#!/bin/bash

for i in $(eval echo {$1..$2})
do
   echo "measuring from realisation $i"
   python measure_from_catalogue.py $3 $4 $i 0 3 2 all
done
