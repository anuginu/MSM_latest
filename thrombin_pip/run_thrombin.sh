#!/bin/bash
#PBS -l select=32:ncpus=272 -lplace=excl
cd /homes/anuginueni/thrombin_pip/
python ga_main_diheds_thrombin.py
