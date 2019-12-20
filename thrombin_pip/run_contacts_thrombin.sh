#!/bin/bash
#PBS -l select=1:ncpus=100 -lplace=excl
cd /homes/anuginueni/thrombin_pip/
python ga_main_contacts_thrombin.py
