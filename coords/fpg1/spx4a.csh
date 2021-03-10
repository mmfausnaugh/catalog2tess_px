#! /bin/csh -f

set codedir=~aml/asmxn/tess2/fpg1

echo $codedir""/program

# v_x, v_y, v_z should be the components of the spacecraft velocity in km/s relative
# to the Solar System barycenter in the J2000 equatorial frame of reference.
# For input code=1, star_table.txt should be an ASCII text file with three columns
# RA(J2000), Dec(J2000), and TMag, one line per star.
# Other input and output formats may be selected - see the code for details.
# $codedir""/stars_aber2 1 1 1 <v_x> <v_y> <v_z> < star_table.txt > stars_aber.txt
$codedir""/stars_aber2 1 1 1 0.0 0.0 0.0 < star_table.txt > stars_aber.txt

# fpg_pars.txt must be a text file with the values of the geometric parameters that accurately
#describe the chosen camera.
# icam = camera number (choices are 1-4)
# ra_sc, dec_sc, roll_sc are parameters (angle in degrees - J2000) describing the spacecraft
# relative to the J2000 equatorial reference frame.  The equivalent 3-2-3 Euler angles
# (rotations around z, y, and z axes) are ra_sc, 90 - dec_sc, and roll_sc + 180.
rm fpg_pars.txt
ln -s fpg_pars2.txt fpg_pars.txt
# $codedir""/starspx4 icam ra_sc dec_sc roll_sc < stars_aber.txt > stars_spx4a.txt
$codedir""/starspx4 1 147.418 -46.028 -214.759 < stars_aber.txt > stars_spx4a.txt

