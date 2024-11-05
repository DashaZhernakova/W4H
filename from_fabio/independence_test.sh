#!/bin/bash
for correction in T F
do
	echo "Doing with correction = ${correction}"
	Rscript "independence_tests.R" $correction
done

