#! /bin/csh -f

rm fpg_pars.txt
ln -s fpg_pars2_m1.txt fpg_pars.txt
~aml/asmxn/tess2/fpg1/chkc1 1 194.562 11.858 0.0 < starsg1.txt > chkc.txt
head chkc.txt


