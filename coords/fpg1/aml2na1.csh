#! /bin/csh -f

set n=1
while($n<5)
    rm fpg_pars.txt
    ln -s fpg_pars2_$n"".txt fpg_pars.txt
    echo Camera $n
    # cat fpg_pars.txt
    # cp fpg_pars2_$n"".txt fl145
    # (~aml/asmxn/tess2/fpg1/aml2na1 $n > fl145/a2$n"".txt) >& fl145/a2$n"".err
    # cp fpg_pars2_$n"".txt fl145
    # (~aml/asmxn/tess2/fpg1/aml2na2 $n > fl145/a3$n"".txt) >& fl145/a3$n"".err
    (~aml/asmxn/tess2/fpg1/aml2na3 $n > fl145/a4$n"".txt) >& fl145/a4$n"".err
    @ n++
end
exit
    set n=1
    rm fpg_pars.txt
    ln -s fpg_pars2_$n""a.txt fpg_pars.txt
    echo Camera $n
    # cat fpg_pars.txt
    cp fpg_pars2_$n""a.txt fl146
    (~aml/asmxn/tess2/fpg1/aml2na1 $n > fl146/a2$n""a.txt) >& fl146/a2$n""a.err



