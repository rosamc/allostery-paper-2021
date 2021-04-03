import sys,os
sys.path.append("./bin")
import utils_matchHill_repo as utils_matchHill
import CG_c4_N6_samesitesFalse
import itertools
import numpy as np
psfunc4=CG_c4_N6_samesitesFalse.interfaceps_a_GRF_CG_c4_N6_samesitesFalse_x
GRFfunc4=CG_c4_N6_samesitesFalse.interface_GRF_CG_c4_N6_samesitesFalse_x


targets=[4,5,6]
betas=[0.25,0.5,0.75]
##initially I had tested more combinations of parameter sets and done like this
#plims_c=[6]
#plims=[4,6]
#conformations=[2,4]

#combis=list(itertools.product(targets,plims,plims_c,conformations,betas)) 
#combi=combis[jid%len(combis)]
#n_target=combi[0]
#c=combi[3]
#beta=combi[4]

#To reproduce the search for the seeds corresponding to plim=4 and 4 conformations:
init_seeds = [0,1,2,6,7,8,12,13,14]
jid=int(sys.argv[1])-1 #range(1,290*9+1)
init_seed=init_seeds[jid%9]
if init_seed<3:
    n_target=4
elif init_seed<10:
    n_target=5
else:
    n_target=6

seeds=np.arange(init_seed,18*290+1,18)
seed=seeds[jid//9]
beta=betas[init_seed%3]

c=4

if c==2:
	#in case of exploring 2 conformations
    npars=13
    idx_c=12
    psfunc=psfunc2
    GRFfunc=GRFfunc2
else:
    psfunc=psfunc4
    GRFfunc=GRFfunc4
    npars=27
    idx_c=24

plim=4
plimc=6
#plimc=8
min_p=1*10**(-plim)
max_p=1*10**(plim)
min_pc=1*10**(-plimc)
max_pc=1*10**(plimc)
print("#n_target:", n_target)
print("#min_p:%g, max_p:%g,min_pc:%g, max_pc=%g"%(min_p,max_p,min_pc,max_pc))
print("seed:%d"%seed)
print("c:%d"%c)

nmax=5000000

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path2=dir_path.replace("/home","/n/scratch3/users/r")
name_save="nt%d_plim%d_plimc%d_c%d_beta%g"%(n_target,plim,plimc,c,beta)
outi=os.path.join(dir_path2,name_save+'_out_%d.out'%jid)
outf=os.path.join("./final","final_%s_%d.out"%(name_save,jid))
bestpars=utils_matchHill.find_parset(n_target,nmax=nmax,seed=seed,min_pos=0.5,max_pos=1.2,min_stp=0.5,max_stp=1.3,npars=npars,min_pc=min_pc, max_pc=max_pc,idx_c=idx_c,min_p=min_p,max_p=max_p,beta=beta,
                                    psfunc=psfunc,GRFfunc=GRFfunc,fitto="curve",printminerror=15,minerror=2,norm=True,fout_intermediate=outi,fout_final=outf)
