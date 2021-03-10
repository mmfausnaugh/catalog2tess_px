#! /bin/csh -f


set codedir=~aml/asmxn/tess2/fpg1

echo $codedir""/program

set icam=$argv[1]

echo camera = $icam

$codedir""/stars_aber2 1 1 1 0.0 0.0 0.0 < gs_table_sector1_camera$icam""_calibration_nohdr.txt > taba1.txt
$codedir""/stars_aber2 6 1 2 0.0 0.0 0.0 < taba1.txt > taba2.txt

$codedir""/stars_aber2 1 1 1 30.0 0.0 0.0 < gs_table_sector1_camera$icam""_calibration_nohdr.txt > tabb1.txt
$codedir""/stars_aber2 6 1 2 30.0 0.0 0.0 < tabb1.txt > tabb2.txt

$codedir""/stars_aber2 1 1 1 0.0 30.0 0.0 < gs_table_sector1_camera$icam""_calibration_nohdr.txt > tabc1.txt
$codedir""/stars_aber2 6 1 2 0.0 30.0 0.0 < tabc1.txt > tabc2.txt

$codedir""/stars_aber2 1 1 1 0.0 0.0 100.0 < gs_table_sector1_camera$icam""_calibration_nohdr.txt > tabd1.txt
$codedir""/stars_aber2 6 1 2 0.0 0.0 100.0 < tabd1.txt > tabd2.txt
