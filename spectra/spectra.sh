#!/bin/tcsh
#
#  regular haswell node: 128GB
#  shared haswell node: max 118GB
#  bigmem: 750GB max
#
#SBATCH --job-name ncremap
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH -C haswell
#XXSBATCH -n 1 --mem 118GB -q shared
#SBATCH -q regular
#XXSBATCH --clusters=escori
#XXSBATCH --qos=bigmem
#XXSBATCH --mem=128GB          

date
#cd ~/scratch2/dyamond/dy2
#cd ~/chrys/F2010y1/run

#unbuffer ncl ~/codes/nclscript/spectra/ke.ncl   't1=22.000' 't2=24.000' 'tinc=0.125'  
ncl ~/codes/nclscript/spectra/ke.ncl  'plvl=500'  't1=1890.' 't2=1920.' 'tinc=0.5'  


date
