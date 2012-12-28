#! /bin/bash
for f in `ls -1 *.pdb.txt`
do
r=`basename $f .pdb.txt`
mv $f $r.dat
done
