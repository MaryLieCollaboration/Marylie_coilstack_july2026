#! /bin/sh
cp $1.mip fort.11
if [ -e $1.aim ]; then
cp $1.aim fort.7
fi
if [ -e $1.var ]; then
cp $1.var fort.8
fi
if [ -e $1.sca ]; then
cp $1.sca fort.3
fi
if [ -e $1.ssv ]; then
cp $1.ssv fort.4
fi
if [ -e $1.sq ]; then
cp $1.sq fort.9
fi
if [ -e $1.sv ]; then
cp $1.sv fort.10
fi
if [ -e $1.inparts ]; then
cp $1.inparts fort.34
fi
if [ -e $1.map ]; then
cp $1.map fort.15
fi
if [ -e $1.cor ]; then
cp $1.cor fort.37
fi
if [ -e $1.map2 ]; then
cp $1.map2 fort.38
fi
ml3_183610.x
if [ -e fort.12 ]; then
cp fort.12 $1.out
fi
if [ -e fort.14 ]; then
mv fort.14 $1.dis
fi
if [ -e fort.16 ]; then
mv fort.16 $1.outmap
fi
if [ -e fort.17 ]; then
mv fort.17 next.mip
fi
if [ -e fort.18 ]; then
mv fort.18 $1.uuu
fi
if [ -e fort.20 ]; then
mv fort.20 $1.wcl
fi
if [ -e fort.22 ]; then
mv fort.22 $1.dmp
fi
if [ -e fort.24 ]; then
mv fort.24 $1.env
fi
if [ -e fort.25 ]; then
mv fort.25 $1.xray
fi
if [ -e fort.26 ]; then
mv fort.26 $1.olp
fi
if [ -e fort.99 ]; then
mv fort.99 $1.tam
fi
if [ -e fort.34 ]; then
mv fort.34 $1.msc
fi
