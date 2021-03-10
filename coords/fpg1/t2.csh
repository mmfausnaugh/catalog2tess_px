#! /bin/csh -f

rm fpg_pars.txt
ln -s fpg_pars1.txt fpg_pars.txt
~aml/asmxn/tess2/fpg1/fpg_a1 1 194.562 11.858 0.0 < starsg1.txt > starsg1_fpg0.txt

set n=1
while ($n <10)
    echo $n
    rm fpg_pars.txt
    ln -s fpg_pars1_m$n"".txt fpg_pars.txt
    ~aml/asmxn/tess2/fpg1/fpg_a1 1 194.562 11.858 0.0 < starsg1.txt > starsg1_fpg$n"".txt
    @ n++
end
