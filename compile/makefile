#=====================================================#
#                                                     #
#     This work is distributed under CC_BY_NC_SA      #
#                                                     #
# Created by Ruben M. Cabezon and Domingo Garcia-Senz #
#               ruben.cabezon@unibas.ch               #
#               domingo.garcia@upc.edu                #
#                                                     #
#                        (2017)                       #
#                                                     #
#=====================================================#
#                                                     #
#                   SPHYNX: makefile                  #
#                                                     #
#=====================================================#

#intel, gnu, ompi
compiler = intel

path = ../sphynx/

#path to simulation directory
#in this way we can have different parameters for different simulations
#and initial datas.
pathsim = ../run_NR/
executable = run

comp = $(strip $(compiler))

ifeq ($(comp),intel)
FC = mpiifort -c -fpp #compile
LK = mpiifort -o $(pathsim)$(executable)  #link
SF = -r8 -shared-intel -qopenmp -xSSE2 -fPIC #-qopt-report
SF2 = $(SF) -mcmodel=large -O3 -DEQMASS #-traceback -CB #-g -check uninit -ftrapuv
else ifeq ($(comp),ompi)
FC = mpif90 -c  #compile
LK = mpif90 -o $(pathsim)$(executable)  #link
SF = -r8 -mcmodel=large -shared-intel -openmp #-qopt-report
SF2 = $(SF) -O3 #-traceback
else ifeq ($(comp),gnu)
FC = mpif90 -c  #compile
LK = mpif90 -o $(pathsim)$(executable)  #link
SF = -fdefault-real-8 -fdefault-double-8 -mcmodel=large -fopenmp
SF2 = $(SF) -O3 #-fbacktrace
else
$(error COMPILER not defined)
endif

CF = $(SF2)  #compiler flags
CF_special = $(SF) -O3 -traceback -qopt-report

OBJS =\
	apply_PBC.o\
	buildtreemod_grav.o\
	calculate_ro0.o\
	calculate_density.o\
	calculate_divv.o\
	calculate_hNR.o\
	calculate_hpowers.o\
	calculate_IAD.o\
	calculate_norm.o\
	calculate_omega.o\
	calculate_switches.o\
	conservemod.o\
	cubicspline.o\
	eosid_t.o\
	eosid_u.o\
	eosmod.o\
	helmholtz_EOS.o\
	const_eos_mod.o\
	helm_table_storage_mod.o\
	read_helm_table.o\
	outputmod.o\
	estabilmod.o\
	findneighbors.o\
	first_findneighbors.o\
	indexx.o\
	indexxi.o\
	init.o\
	init_scenario.o\
	sphynx_hybrid.o\
	masscenter.o\
	momeqnmod.o\
	parameters.o\
	printreport.o\
	profile.o\
	readdata.o\
	testparameters.o\
	timectrlmod.o\
	treewalkmod_grav_mefec.o\
	update.o\
	wkernel.o\
	wkernel_noderiv.o
#...... SPH ......

sphynxSPH:\
  $(OBJS); $(LK) $(CF) $(OBJS)

apply_PBC.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

buildtreemod_grav.o:\
  parameters.o\
  profile.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_ro0.o:\
  parameters.o\
  apply_PBC.o\
  wkernel.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_density.o:\
  parameters.o\
  apply_PBC.o\
  wkernel.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_divv.o:\
  parameters.o\
  wkernel_noderiv.o\
  profile.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_hNR.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_hpowers.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_IAD.o:\
  parameters.o\
  wkernel_noderiv.o\
  profile.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_norm.o:\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_omega.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

calculate_switches.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

conservemod.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

const_eos_mod.o:\
  $(path)EOS/$(@:.o=.f90); $(FC) $(CF) $(path)EOS/$(@:.o=.f90)

cubicspline.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

eosmod.o:\
  parameters.o\
  profile.o\
  eosid_t.o\
  eosid_u.o\
  helmholtz_EOS.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

eosid_t.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

eosid_u.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

outputmod.o:\
  parameters.o\
  indexx.o\
  $(pathsim)$(@:.o=.f90); $(FC) $(CF) $(pathsim)$(@:.o=.f90)

estabilmod.o:\
  parameters.o\
  indexx.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

findneighbors.o:\
  parameters.o\
  apply_PBC.o\
  calculate_hpowers.o\
  profile.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

first_findneighbors.o:\
  parameters.o\
  apply_PBC.o\
  calculate_hpowers.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

helmholtz_EOS.o:\
  parameters.o\
  const_eos_mod.o\
  helm_table_storage_mod.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

helm_table_storage_mod.o:\
  $(path)EOS/$(@:.o=.f90); $(FC) $(CF) $(path)EOS/$(@:.o=.f90)

indexx.o:\
  $(path)$(@:.o=.f); $(FC) $(CF) $(path)$(@:.o=.f)

indexxi.o:\
  $(path)$(@:.o=.f); $(FC) $(CF) $(path)$(@:.o=.f)

init.o:\
  parameters.o\
  cubicspline.o\
  calculate_norm.o\
  read_helm_table.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

init_scenario.o:\
  parameters.o\
  indexx.o\
  $(pathsim)$(@:.o=.f90); $(FC) $(CF) $(pathsim)$(@:.o=.f90)

sphynx_hybrid.o:\
  parameters.o\
  profile.o\
  testparameters.o\
  printreport.o\
  init.o\
  readdata.o\
  init_scenario.o\
  buildtreemod_grav.o\
  treewalkmod_grav_mefec.o\
  first_findneighbors.o\
  findneighbors.o\
  calculate_ro0.o\
  calculate_density.o\
  calculate_hNR.o\
  calculate_omega.o\
  calculate_hpowers.o\
  calculate_divv.o\
  calculate_IAD.o\
  calculate_switches.o\
  eosmod.o\
  momeqnmod.o\
  update.o\
  timectrlmod.o\
  conservemod.o\
  outputmod.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

masscenter.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

momeqnmod.o:\
  parameters.o\
  profile.o\
  wkernel_noderiv.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

parameters.o:\
  $(pathsim)$(@:.o=.f90); $(FC) $(CF) $(pathsim)$(@:.o=.f90)

printreport.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

profile.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

readdata.o:\
  parameters.o\
  calculate_hpowers.o\
  masscenter.o\
  $(pathsim)$(@:.o=.f90); $(FC) $(CF) $(pathsim)$(@:.o=.f90)

read_helm_table.o:\
  helm_table_storage_mod.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

testparameters.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

timectrlmod.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

treewalkmod_grav_mefec.o:\
  parameters.o\
  profile.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

update.o:\
  parameters.o\
  profile.o\
  masscenter.o\
  estabilmod.o\
  calculate_hpowers.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

wkernel.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)

wkernel_noderiv.o:\
  parameters.o\
  $(path)$(@:.o=.f90); $(FC) $(CF) $(path)$(@:.o=.f90)


#cleaning
clean:
	rm *.o *.mod

clena:
	rm *.o *.mod
