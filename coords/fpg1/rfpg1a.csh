#! /bin/csh -f


set codedir=~aml/asmxn/tess2/fpg1
echo $codedir""/program

set datadir=/nfs/tess/data/roland/FP/

set icam=$argv[1]

echo camera = $icam

$codedir""/stars_aber1 1 23.30373901 -15.76175465 -7.02924008 < $datadir""/gsc_cam1_radec.ntxt > $datadir""/gsc_cam1_aber.txt

rm fpg_pars.txt
ln -s fpg_pars2.txt fpg_pars.txt
$codedir""/starspx4 $icam 250.546599604 32.4748335111 -13.0205693574 < $datadir""/gsc_cam1_aber.txt > $datadir""/gsc_cam1_px4.txt
