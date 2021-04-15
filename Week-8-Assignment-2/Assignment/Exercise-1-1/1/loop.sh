#! /bin/bash
BIN='/home/JSWL/bin/vasp_ncl'
rm WAVECAR SUMMARY.stanene
for i in 4.64 4.65 4.66 4.67 4.68 4.69 4.70 4.71 4.72 ; do
for j in 0.5000 0.5050 0.5080 0.5085 0.5090 0.5095 0.5100 ; do
cat >POSCAR  <<!
stanene
$i
        1.0000000000         0.0000000000         0.0000000000
       -0.5000000000         0.8660254037         0.0000000000
        0.0000000000         0.0000000000        20.0000000000/$i
   Sn
    2
Direct
     0.333333349         0.666666698         0.500000000
     0.666666675         0.333333349         $j
!
echo "a= $i  z= $j" ; mpirun -np 16 $BIN
E=`awk '/F=/ {print $0}' OSZICAR` ; echo $i $j $E >>SUMMARY.stanene
done
echo "" >>SUMMARY.stanene
done
cat SUMMARY.stanene
