#!/bin/bash

for dld in $(seq 110 202); do
	.build/ptssdse ${dld} >> dump.txt
done
