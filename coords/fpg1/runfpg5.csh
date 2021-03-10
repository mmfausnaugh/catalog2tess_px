#! /bin/csh -f


set codedir=~aml/asmxn/tess2/fpg1

echo $codedir""/program

set icam=$argv[1]

echo camera = $icam

$codedir""/stars_aber1 1 0.0 0.0 0.0 < star_table_sector1_camera1ed1.txt > stars_sec1_cam1ed1_aber.txt

rm fpg_pars.txt
ln -s fpg_pars2_m1.txt fpg_pars.txt
$codedir""/starspx4 $icam 147.418 -46.028 -214.759 < stars_sec1_cam1ed1_aber.txt > stars_tsc1_spx4.txt 
$codedir""/starsig3 0.7 0.07 13.0 10 < stars_tsc1_spx4.txt > stars_tsc1_sig3.txt

rm fpg_pars.txt
ln -s fpg_pars2_eul0.txt fpg_pars.txt
($codedir""/runfpg4 $icam 147.418 -46.028 -214.759 0.03 0.1 1 50.0 < stars_tsc1_sig3.txt > runfpg4a.txt) >& rf4a.serr

mv runfpg_stars.txt runfpg_stars4a.txt
head -3 runfpg4a.txt > fpg_pars4a.txt
tail -n +4 fpg_pars2.txt >> fpg_pars4a.txt

rm fpg_pars.txt
ln -s fpg_pars4a.txt fpg_pars.txt
($codedir""/runfpg4 $icam 147.418 -46.028 -214.759 0.03 0.1 1 50.0 < stars_tsc1_sig3.txt > runfpg4b.txt) >& rf4b.serr
