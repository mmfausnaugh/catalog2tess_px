#! /bin/csh -f

rm fpg_pars.txt
ln -s fpg_pars2_m1.txt fpg_pars.txt
~aml/asmxn/tess2/fpg1/starspx2 1 194.562 11.858 0.0 < starsg1.txt > starsg1_spx2.txt

~aml/asmxn/tess2/fpg1/getpix1 < starsg1_spx2.txt > pixels1.txt
(~aml/asmxn/tess2/fpg1/radec1 1 194.562 11.858 0.0 < pixels1.txt > radec1_1.txt) >& ra.serr
~aml/asmxn/tess2/fpg1/fpg_a1 1 194.562 11.858 0.0 < starsg1.txt > starsg1_fpga1_1.txt
echo starsg1_spx2.txt
head starsg1_spx2.txt
echo radec1_1.txt
head radec1_1.txt
echo ra.serr
head -50 ra.serr


