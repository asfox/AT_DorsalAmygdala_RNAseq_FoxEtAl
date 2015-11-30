#!/bin/bash

fslroi $1 $2 0 109 0 148 0 109 
fslswapdim $2 -x y z $2 
gzip $2
