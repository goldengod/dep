#!/bin/bash

dirs="../tipo00/* ../tipo01/* ../tipo10/* ../tipo11/* ../RND/* ../LOLIB/*"
files=$(ls -l $dirs | sed -r 's/\s+/ /g' | cut -d' ' -f9 | grep -v "greedy$")

for inp in $files; do
	out="${inp}.greedy"
	echo "--- INPUT: ${inp} --- OUTPUT: ${out} ---"
	./greedy $inp > $out
	echo "---"
done


