#!/bin/bash
p=$1
/usr/bin/time -f "p=$p wall=%e s" ../omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -p=$p -square > s$p.tmp 2> s$p.err && mv s$p.tmp s$p.txt
echo "$p rc=$? $(tail -1 s$p.err)" >> done.txt
