#! /bin/csh -f

#rm fpg_pars.txt
#ln -s fpg_pars2_m1.txt fpg_pars.txt
#~aml/asmxn/tess2/fpg1/starspx1 1 194.562 11.858 0.0 < stars1.txt > stars1_spx1.txt
#~aml/asmxn/tess2/fpg1/starsig1 0.5 < stars1_spx1.txt > stars1_sig1.txt

#~aml/asmxn/tess2/fpg1/starspx2 1 194.562 11.858 0.0 < starsg1.txt > starsg2_spx1.txt
#~aml/asmxn/tess2/fpg1/starsig2 0.7 0.07 13.0 10 < starsg2_spx1.txt > starsg2_sig1.txt

# rm fpg_pars.txt
# ln -s fpg_pars2.txt fpg_pars.txt
# (~aml/asmxn/tess2/fpg1/runfpg2 1 194.562 11.858 0.0 0.0 0.0 0.0 0.03 0.1 1 50.0 < starsg2_sig1.txt > runfpg2.txt) >& rf2.serr

rm fpg_pars.txt
ln -s fpg_pars2_m1.txt fpg_pars.txt
~aml/asmxn/tess2/fpg1/starspx3 1 147.418 -46.028 -214.759 < star_table_sector1_camera1ed1.txt > stars_tsc1_spx2.txt
~aml/asmxn/tess2/fpg1/starsig2 0.7 0.07 13.0 10 < stars_tsc1_spx2.txt > stars_tsc1_sig2.txt

rm fpg_pars.txt
ln -s fpg_pars2_eul0.txt fpg_pars.txt
(~aml/asmxn/tess2/fpg1/runfpg3 1 147.418 -46.028 -214.759 0.0 0.0 0.0 0.03 0.1 1 50.0 < stars_tsc1_sig2.txt > runfpg3a.txt) >& rf3a.serr

head -3 runfpg3a.txt > fpg_pars3a.txt
tail -n +4 fpg_pars2.txt >> fpg_pars3a.txt

rm fpg_pars.txt
ln -s fpg_pars3a.txt fpg_pars.txt
(~aml/asmxn/tess2/fpg1/runfpg3 1 147.418 -46.028 -214.759 0.0 0.0 0.0 0.03 0.1 1 50.0 < stars_tsc1_sig2.txt > runfpg3b.txt) >& rf3b.serr
