#!/bin/bash
module load python/3.6.3
module load python/3.6.3/matplotlib/1.5.1
module load python/3.6.3/numpy/1.14.1
module load bedtools/2.28.0

python3 /scratch/Shares/rinn/Michael/firre_timecourse/proseq/bin/mumerge/mumerge/mumerge.py \
-i sample_info.txt \
-o split_bidir_predictions_mumerge

