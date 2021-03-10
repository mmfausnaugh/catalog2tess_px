AVGFPG1_OBJS = avgfpg1.o

STATFPG1_OBJS = statfpg1.o

STARFPG1_OBJS = starfpg1.o

STARFPGDIF1_OBJS = starfpgdif1.o

MODIFY_PARS1_OBJS = modify_pars1.o twister2.o seed_tw_ran.o

MODIFY_PARS2_OBJS = modify_pars2.o twister2.o seed_tw_ran.o

FPG_A1_OBJS = fpg_a1.o mat_ra3.o vec.o

FPG_A2_OBJS = fpg_a2.o mat_ra3.o vec.o

RADEC1_OBJS = radec1.o mat_ra3.o vec.o

RADEC2_OBJS = radec2.o mat_ra3.o vec.o

AML2NA1_OBJS = aml2na1.o mat_ra3.o vec.o

AML2NA2_OBJS = aml2na2.o mat_ra3.o vec.o

AML2NA3_OBJS = aml2na3.o mat_ra3.o vec.o

CHKC1_OBJS = chkc1.o mat_ra3.o vec.o

RUNFPG1_OBJS = runfpg1.o fpg1.o mrq_fpg1.o mat_ra3.o vec.o

RUNFPG2_OBJS = runfpg2.o fpg2.o mrq_fpg2.o mat_ra3.o vec.o

RUNFPG3_OBJS = runfpg3.o fpg3.o mrq_fpg3.o mat_ra3.o vec.o

RUNFPG4_OBJS = runfpg4.o fpg4.o mrq_fpg4.o mat_ra3.o vec.o

RUNFPG5_OBJS = runfpg5.o fpg4.o mrq_fpg5.o mat_ra3.o vec.o

FPGTR1_OBJS = fpgtr1.o fpg4.o mrq_fpg5.o mat_ra3.o vec.o

TESTTWISTER2_OBJS = testtwister2.o twister2.o seed_tw_ran.o

ABER_STARS1_OBJS = aber_stars1.o vec.o

STARS_ABER1_OBJS = stars_aber1.o vec.o

STARS_ABER2_OBJS = stars_aber2.o vec.o

STARSPX1_OBJS = starspx1.o mat_ra3.o vec.o

STARSPX2_OBJS = starspx2.o mat_ra3.o vec.o

STARSPX3_OBJS = starspx3.o mat_ra3.o vec.o

STARSPX4_OBJS = starspx4.o mat_ra3.o vec.o

STARSPX5_OBJS = starspx5.o mat_ra3.o vec.o

STARSPX6_OBJS = starspx6.o mat_ra3.o vec.o

STARSPX7_OBJS = starspx7.o mat_ra3.o vec.o

STARSPX8_OBJS = starspx8.o mat_ra3.o vec.o

STARSPX9_OBJS = starspx9.o mat_ra3.o vec.o

STARS_PX5_OBJS = stars_px5.o mat_ra3.o vec.o

STARSIG1_OBJS = starsig1.o twister2.o seed_tw_ran.o

STARSIG2_OBJS = starsig2.o twister2.o seed_tw_ran.o

STARSIG3_OBJS = starsig3.o twister2.o seed_tw_ran.o

STGRID1_OBJS = stgrid1.o mat_ra3.o vec.o

STDIF1_OBJS = stdif1.o

ANRES1_OBJS = anres1.o

ANRES2_OBJS = anres2.o

ANRES3_OBJS = anres3.o

RMDUP1_OBJS = rmdup1.o

CROSSCHK1_OBJS = crosschk1.o

RUNSUNPOS1_OBJS = runsunpos1.o astron_rb.o mat_ra3.o vec.o

COMBEUL1_OBJS = combeul1.o mat_ra3.o vec.o

ORIENT1_OBJS = orient1.o mat_ra3.o vec.o

GETPIX1_OBJS = getpix1.o

GETRDRDM1_OBJS = getrdrdm1.o

GETRDM1_OBJS = getrdm1.o

RADLDST1_OBJS = radldst1.o

RHIST1_OBJS = rhist1.o

ROFT1_OBJS = roft1.o

GET3COLRDM1_OBJS = get3colrdm1.o

all: fpg1

vec.o:			vec.c vec.h
			gcc -c -o vec.o vec.c

mat_ra3.o:		mat_ra3.c mat_ra3.h
			gcc -c -o mat_ra3.o mat_ra3.c

astron_rb.o:		astron_rb.c
			gcc -c -o astron_rb.o astron_rb.c

fpg_a1.o:		fpg_a1.c
			gcc -c -o fpg_a1.o fpg_a1.c

fpg_a1:			${FPG_A1_OBJS}
			gcc -o fpg_a1 ${FPG_A1_OBJS} -lm

fpg_a2.o:		fpg_a2.c
			gcc -c -o fpg_a2.o fpg_a2.c

fpg_a2:			${FPG_A2_OBJS}
			gcc -o fpg_a2 ${FPG_A2_OBJS} -lm

aber_stars1.o:		aber_stars1.c
			gcc -c -o aber_stars1.o aber_stars1.c

aber_stars1:		${ABER_STARS1_OBJS}
			gcc -o aber_stars1 ${ABER_STARS1_OBJS} -lm

stars_aber1.o:		stars_aber1.c
			gcc -c -o stars_aber1.o stars_aber1.c

stars_aber1:		${STARS_ABER1_OBJS}
			gcc -o stars_aber1 ${STARS_ABER1_OBJS} -lm

stars_aber2.o:		stars_aber2.c
			gcc -c -o stars_aber2.o stars_aber2.c

stars_aber2:		${STARS_ABER2_OBJS}
			gcc -o stars_aber2 ${STARS_ABER2_OBJS} -lm

radec1.o:		radec1.c
			gcc -c -o radec1.o radec1.c

radec1:			${RADEC1_OBJS}
			gcc -o radec1 ${RADEC1_OBJS} -lm

radec2.o:		radec2.c
			gcc -c -o radec2.o radec2.c

radec2:			${RADEC2_OBJS}
			gcc -o radec2 ${RADEC2_OBJS} -lm

aml2na1.o:		aml2na1.c
			gcc -c -o aml2na1.o aml2na1.c

aml2na1:		${AML2NA1_OBJS}
			gcc -o aml2na1 ${AML2NA1_OBJS} -lm

aml2na2.o:		aml2na2.c
			gcc -c -o aml2na2.o aml2na2.c

aml2na2:		${AML2NA2_OBJS}
			gcc -o aml2na2 ${AML2NA2_OBJS} -lm

aml2na3.o:		aml2na3.c
			gcc -c -o aml2na3.o aml2na3.c

aml2na3:		${AML2NA3_OBJS}
			gcc -o aml2na3 ${AML2NA3_OBJS} -lm

chkc1.o:		chkc1.c
			gcc -c -o chkc1.o chkc1.c

chkc1:			${CHKC1_OBJS}
			gcc -o chkc1 ${CHKC1_OBJS} -lm

mrq_fpg1.o:		mrq_fpg1.c
			gcc -c -o mrq_fpg1.o mrq_fpg1.c

mrq_fpg2.o:		mrq_fpg2.c
			gcc -c -o mrq_fpg2.o mrq_fpg2.c

mrq_fpg3.o:		mrq_fpg3.c
			gcc -c -o mrq_fpg3.o mrq_fpg3.c

mrq_fpg4.o:		mrq_fpg4.c
			gcc -c -o mrq_fpg4.o mrq_fpg4.c

mrq_fpg5.o:		mrq_fpg5.c
			gcc -c -o mrq_fpg5.o mrq_fpg5.c

fpg1.o:			fpg1.c
			gcc -c -o fpg1.o fpg1.c

fpg2.o:			fpg2.c
			gcc -c -o fpg2.o fpg2.c

fpg3.o:			fpg3.c
			gcc -c -o fpg3.o fpg3.c

fpg4.o:			fpg4.c
			gcc -c -o fpg4.o fpg4.c

runfpg1.o:		runfpg1.c
			gcc -c -o runfpg1.o runfpg1.c

runfpg1:		${RUNFPG1_OBJS}
			gcc -o runfpg1 ${RUNFPG1_OBJS} -lm

runfpg2.o:		runfpg2.c
			gcc -c -o runfpg2.o runfpg2.c

runfpg2:		${RUNFPG2_OBJS}
			gcc -o runfpg2 ${RUNFPG2_OBJS} -lm

runfpg3.o:		runfpg3.c
			gcc -c -o runfpg3.o runfpg3.c

runfpg3:		${RUNFPG3_OBJS}
			gcc -o runfpg3 ${RUNFPG3_OBJS} -lm

runfpg4.o:		runfpg4.c
			gcc -c -o runfpg4.o runfpg4.c

runfpg4:		${RUNFPG4_OBJS}
			gcc -o runfpg4 ${RUNFPG4_OBJS} -lm

runfpg5.o:		runfpg5.c
			gcc -c -o runfpg5.o runfpg5.c

runfpg5:		${RUNFPG5_OBJS}
			gcc -o runfpg5 ${RUNFPG5_OBJS} -lm

fpgtr1.o:		fpgtr1.c
			gcc -c -o fpgtr1.o fpgtr1.c

fpgtr1:			${FPGTR1_OBJS}
			gcc -o fpgtr1 ${FPGTR1_OBJS} -lm

twister2.o:		twister2.c
			gcc -c -o twister2.o twister2.c

seed_tw_ran.o:		seed_tw_ran.c
			gcc -c -o seed_tw_ran.o seed_tw_ran.c 

testtwister2.o:		testtwister2.c
			gcc -c -o testtwister2.o testtwister2.c

testtwister2:		${TESTTWISTER2_OBJS}
			gcc -o testtwister2 ${TESTTWISTER2_OBJS} -lm

modify_pars1.o:		modify_pars1.c
			gcc -c -o modify_pars1.o modify_pars1.c

modify_pars1:		${MODIFY_PARS1_OBJS}
			gcc -o modify_pars1 ${MODIFY_PARS1_OBJS} -lm

modify_pars2.o:		modify_pars2.c
			gcc -c -o modify_pars2.o modify_pars2.c

modify_pars2:		${MODIFY_PARS2_OBJS}
			gcc -o modify_pars2 ${MODIFY_PARS2_OBJS} -lm

starspx1.o:		starspx1.c
			gcc -c -o starspx1.o starspx1.c

starspx1:		${STARSPX1_OBJS}
			gcc -o starspx1 ${STARSPX1_OBJS} -lm

starspx2.o:		starspx2.c
			gcc -c -o starspx2.o starspx2.c

starspx2:		${STARSPX2_OBJS}
			gcc -o starspx2 ${STARSPX2_OBJS} -lm

starspx3.o:		starspx3.c
			gcc -c -o starspx3.o starspx3.c

starspx3:		${STARSPX3_OBJS}
			gcc -o starspx3 ${STARSPX3_OBJS} -lm

starspx4.o:		starspx4.c
			gcc -c -o starspx4.o starspx4.c

starspx4:		${STARSPX4_OBJS}
			gcc -o starspx4 ${STARSPX4_OBJS} -lm

starspx5.o:		starspx5.c
			gcc -c -o starspx5.o starspx5.c

starspx5:		${STARSPX5_OBJS}
			gcc -o starspx5 ${STARSPX5_OBJS} -lm

starspx6.o:		starspx6.c
			gcc -c -o starspx6.o starspx6.c

starspx6:		${STARSPX6_OBJS}
			gcc -o starspx6 ${STARSPX6_OBJS} -lm

starspx7.o:		starspx7.c
			gcc -c -o starspx7.o starspx7.c

starspx7:		${STARSPX7_OBJS}
			gcc -o starspx7 ${STARSPX7_OBJS} -lm

starspx8.o:		starspx8.c
			gcc -c -o starspx8.o starspx8.c

starspx8:		${STARSPX8_OBJS}
			gcc -o starspx8 ${STARSPX8_OBJS} -lm

starspx9.o:		starspx9.c
			gcc -c -o starspx9.o starspx9.c

starspx9:		${STARSPX9_OBJS}
			gcc -o starspx9 ${STARSPX9_OBJS} -lm

stars_px5.o:		stars_px5.c
			gcc -c -o stars_px5.o stars_px5.c

stars_px5:		${STARS_PX5_OBJS}
			gcc -o stars_px5 ${STARS_PX5_OBJS} -lm

starsig1.o:		starsig1.c
			gcc -c -o starsig1.o starsig1.c

starsig1:		${STARSIG1_OBJS}
			gcc -o starsig1 ${STARSIG1_OBJS} -lm

starsig2.o:		starsig2.c
			gcc -c -o starsig2.o starsig2.c

starsig2:		${STARSIG2_OBJS}
			gcc -o starsig2 ${STARSIG2_OBJS} -lm

starsig3.o:		starsig3.c
			gcc -c -o starsig3.o starsig3.c

starsig3:		${STARSIG3_OBJS}
			gcc -o starsig3 ${STARSIG3_OBJS} -lm

stgrid1.o:		stgrid1.c
			gcc -c -o stgrid1.o stgrid1.c

stgrid1:		${STGRID1_OBJS}
			gcc -o stgrid1 ${STGRID1_OBJS} -lm

stdif1.o:		stdif1.c
			gcc -c -o stdif1.o stdif1.c

stdif1:			${STDIF1_OBJS}
			gcc -o stdif1 ${STDIF1_OBJS} -lm

anres1.o:		anres1.c
			gcc -c -o anres1.o anres1.c

anres1:			${ANRES1_OBJS}
			gcc -o anres1 ${ANRES1_OBJS} -lm

anres2.o:		anres2.c
			gcc -c -o anres2.o anres2.c

anres2:			${ANRES2_OBJS}
			gcc -o anres2 ${ANRES2_OBJS} -lm

anres3.o:		anres3.c
			gcc -c -o anres3.o anres3.c

anres3:			${ANRES3_OBJS}
			gcc -o anres3 ${ANRES3_OBJS} -lm

rmdup1.o:		rmdup1.c
			gcc -c -o rmdup1.o rmdup1.c

rmdup1:			${RMDUP1_OBJS}
			gcc -o rmdup1 ${RMDUP1_OBJS} -lm

crosschk1.o:		crosschk1.c
			gcc -c -o crosschk1.o crosschk1.c

crosschk1:		${CROSSCHK1_OBJS}
			gcc -o crosschk1 ${CROSSCHK1_OBJS} -lm

runsunpos1.o:		runsunpos1.c
			gcc -c -o runsunpos1.o runsunpos1.c

runsunpos1:		${RUNSUNPOS1_OBJS}
			gcc -o runsunpos1 ${RUNSUNPOS1_OBJS} -lm

combeul1.o:		combeul1.c
			gcc -c -o combeul1.o combeul1.c

combeul1:		${COMBEUL1_OBJS}
			gcc -o combeul1 ${COMBEUL1_OBJS} -lm

avgfpg1.o:		avgfpg1.c
			gcc -c -o avgfpg1.o avgfpg1.c

avgfpg1:		${AVGFPG1_OBJS}
			gcc -o avgfpg1 ${AVGFPG1_OBJS} -lm

statfpg1.o:		statfpg1.c
			gcc -c -o statfpg1.o statfpg1.c

statfpg1:		${STATFPG1_OBJS}
			gcc -o statfpg1 ${STATFPG1_OBJS} -lm

starfpg1.o:		starfpg1.c
			gcc -c -o starfpg1.o starfpg1.c

starfpg1:		${STARFPG1_OBJS}
			gcc -o starfpg1 ${STARFPG1_OBJS} -lm

starfpgdif1.o:		starfpgdif1.c
			gcc -c -o starfpgdif1.o starfpgdif1.c

starfpgdif1:		${STARFPGDIF1_OBJS}
			gcc -o starfpgdif1 ${STARFPGDIF1_OBJS} -lm

orient1.o:		orient1.c
			gcc -c -o orient1.o orient1.c

orient1:		${ORIENT1_OBJS}
			gcc -o orient1 ${ORIENT1_OBJS} -lm

getpix1.o:		getpix1.c
			gcc -c -o getpix1.o getpix1.c

getpix1:		${GETPIX1_OBJS}
			gcc -o getpix1 ${GETPIX1_OBJS} -lm

getrdrdm1.o:		getrdrdm1.c
			gcc -c -o getrdrdm1.o getrdrdm1.c

getrdrdm1:		${GETRDRDM1_OBJS}
			gcc -o getrdrdm1 ${GETRDRDM1_OBJS} -lm

getrdm1.o:		getrdm1.c
			gcc -c -o getrdm1.o getrdm1.c

getrdm1:		${GETRDM1_OBJS}
			gcc -o getrdm1 ${GETRDM1_OBJS} -lm

radldst1.o:		radldst1.c
			gcc -c -o radldst1.o radldst1.c

radldst1:		${RADLDST1_OBJS}
			gcc -o radldst1 ${RADLDST1_OBJS} -lm

rhist1.o:		rhist1.c
			gcc -c -o rhist1.o rhist1.c

rhist1:			${RHIST1_OBJS}
			gcc -o rhist1 ${RHIST1_OBJS} -lm

roft1.o:		roft1.c
			gcc -c -o roft1.o roft1.c

roft1:			${ROFT1_OBJS}
			gcc -o roft1 ${ROFT1_OBJS} -lm

get3colrdm1.o:		get3colrdm1.c
			gcc -c -o get3colrdm1.o get3colrdm1.c

get3colrdm1:		${GET3COLRDM1_OBJS}
			gcc -o get3colrdm1 ${GET3COLRDM1_OBJS} -lm

