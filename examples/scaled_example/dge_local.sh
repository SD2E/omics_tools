#!/bin/bash

chmod -R 777 scripts/

for i in {0..12}
do
    echo i=$i
    date
    nohup ./scripts/dge_$i.sh &
    date
done