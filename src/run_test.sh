#!/bin/bash

BASE='/biggin/b240/ball5406/isocitrate_lyase/ICL2/ICL2_apo/fullprotein'

python3 -t "$BASE/Run1/md_prot.gro" -x $BASE/Run1/md_prot_cent_ts1000.xtc -ref $BASE/Run1/6edw_fix.pdb
