#!/bin/bash
#----------------------------------------------------------
#SBATCH -p defq
###SBATCH -p compute-c6400
#SBATCH -J pandurata_ULX
#SBATCH -o pandurata_ULX.o%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
###SBATCH -t 24:00:00
#SBATCH --mail-user=fragilep@cofc.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#----------------------------------------------------------

./pandurata
