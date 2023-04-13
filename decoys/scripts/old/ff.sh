#!/bin/bash


for i in {0..9}
do
  echo "nohup python ~/xray/md_xray/scripts/ff.py $i &"
  nohup python ~/xray/md_xray/scripts/ff.py $i &
done