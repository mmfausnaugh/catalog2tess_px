#! /bin/csh -f

rm fpg_pars.txt
ln -s fpg_pars1_m1.txt fpg_pars.txt
#~aml/asmxn/tess2/fpg1/starspx1 1 194.562 11.858 0.0 < stars1.txt > stars1_spx1.txt
#~aml/asmxn/tess2/fpg1/starsig1 0.5 < stars1_spx1.txt > stars1_sig1.txt

~aml/asmxn/tess2/fpg1/starspx1 1 194.562 11.858 0.0 < starsg1.txt > starsg1_spx1.txt
# ~aml/asmxn/tess2/fpg1/starsig1 0.5 < starsg1_spx1.txt > starsg1_sig1.txt
~aml/asmxn/tess2/fpg1/starsig2 0.7 0.07 13.0 10 < starsg1_spx1.txt > starsg2_sig1.txt

rm fpg_pars.txt
ln -s fpg_pars1.txt fpg_pars.txt
(~aml/asmxn/tess2/fpg1/runfpg1 1 194.562 11.858 0.0 0.0 0.0 0.0 0.03 0.2 1 50.0 < starsg2_sig1.txt > runfpg2.txt) >& rf2.serr
