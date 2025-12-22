#!/usr/bin/bash 
apptainer run --nv -B /ddn/:/ddn docker://genomicpariscentre/guppy-gpu:6.3.8 guppy_basecaller -i $IN -s $OUT --flowcell FLO-MIN114 --kit SQK-LSK114 -x cuda:0 --fast5_out