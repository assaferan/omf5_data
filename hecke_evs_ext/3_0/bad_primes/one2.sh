#!/bin/sh
# one.sh d p: T_p and T_{p^2} at a prime p dividing d, omf5 main f3e0c64
cd /scratch/home/assaferan/g20/badp
/usr/bin/time -f "%e s %M KB" -o out2/$1.$2.time ./omf5_f3e0c64 -genus=/scratch/home/assaferan/g20/qf5db.sage -isom -disc=$1 -p=$2 -square > out2/$1.$2.dat 2> out2/$1.$2.log
echo "$1 $2 rc=$? $(tail -1 out2/$1.$2.time)" >> done2.txt
