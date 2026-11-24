#!/bin/bash
p=$1; if [ $((p*p)) -le 1084 ]; then sq=-square; else sq=; fi
/usr/bin/time -f "p=$p wall=%e" ../../omf5/src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -p=$p $sq > q$p.tmp 2> q$p.err && mv q$p.tmp q$p.txt
