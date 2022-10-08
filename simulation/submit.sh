#!/bin/bash

for av in 0 2; do
	for J in 1 6; do
		dir=short
		for L in 4 6 8 10 12; do
			name="L_${L}_dJ_${J}_${dir}_av_${av}_llr"
			sed "s/NAME/${name}/g" job_script.sh > job_script_${name}.sh
			sed "s/LENGTH/${L}/g; s/CHAOS/${J}/g; s/DIR/${dir}/g; s/OFFSET/0/g; s/AVOID/${av}/g" llr_vs_reweighting.R > llr_vs_reweighting_${name}.R
			sbatch job_script_${name}.sh llr_vs_reweighting_${name}.R
		done

		dir=long
		L=4
		name="L_${L}_dJ_${J}_${dir}_av_${av}_llr"
		sed "s/NAME/${name}/g" job_script.sh > job_script_${name}.sh
		sed "s/LENGTH/${L}/g; s/CHAOS/${J}/g; s/DIR/${dir}/g; s/OFFSET/16/g; s/AVOID/${av}/g" llr_vs_reweighting.R > llr_vs_reweighting_${name}.R
		sbatch job_script_${name}.sh llr_vs_reweighting_${name}.R
	done
done
