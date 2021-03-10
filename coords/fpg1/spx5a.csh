#! /bin/csh -f


set codedir=~aml/asmxn/tess2/fpg1

echo $codedir""/program

$codedir""/stars_aber1 1 0.0 0.0 0.0 < star_table_sector1_camera1ed1.txt > stars_sec1_cam1ed1_aber.txt

rm fpg_pars.txt
ln -s fpg_pars2.txt fpg_pars.txt
$codedir""/starspx4 1 147.418 -46.028 -214.759 < stars_sec1_cam1ed1_aber.txt > stars_tsc1_spx4a.txt 
$codedir""/stars_px5 1 147.418 -46.028 -214.759 < stars_sec1_cam1ed1_aber.txt > stars_tsc1_spx5a.txt 
