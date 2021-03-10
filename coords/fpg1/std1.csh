#! /bin/csh -f

set n=1
while($n<10)
    ~aml/asmxn/tess2/fpg1/stdif1 starsg1_fpg0.txt starsg1_fpg$n"".txt 0.015 > std_g1_$n"".txt
    mv std_xy.txt stdg1_xy$n"".txt
    mv std_ccd.txt stdg1_ccd$n"".txt
    mv std_fits.txt stdg1_fits$n"".txt
    @ n++
end




