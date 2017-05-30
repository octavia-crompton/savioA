#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Account:
#SBATCH --account=fc_firewater
#
# Partition:
#SBATCH --partition=savio
#
# Wall clock limit:
#SBATCH --time=00:30:00
#
## Command(s) to run:
python build_model.py B2i_plane_randv_Sf 0



