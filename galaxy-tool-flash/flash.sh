#!/bin/bash

outlocation=$(mktemp -d /home/galaxy/galaxy/database/XXXXXX)
SCRIPTDIR=$(dirname "$(readlink -f "$0")")

python $SCRIPTDIR"/flash_wrapper.py" -i $1 -of $outlocation -t $4 -f $5 -m $6 -x $7 -M $8
mv $outlocation"/log.log" $3
mv $outlocation"/merged.zip" $2
rm -rf $outlocation 
