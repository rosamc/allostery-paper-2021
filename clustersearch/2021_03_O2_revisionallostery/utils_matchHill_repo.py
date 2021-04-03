import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append("/home/rm335/repos/sharedposstpNov19/GeneRegulatoryFunctions/utilsGRF") #path to utilsGRF in cluster
sys.path.append("/Users/rosamartinezcorral/Dropbox (HMS)/work2/shared/utilsGRF") #path to utilsGRF locally
import BoundaryFinder as BF

def hill(n,xar):
    return [x**n/(1+x**n) for x in xar]

def MWC(pars,x,n=4):
    K_1,K_2,l2=pars
    alpha=x*K_1
    c_alpha=x*K_2
    L=l2
    f=(alpha*(1+alpha)**(n-1)+c_alpha*L*(1+c_alpha)**(n-1))/((1+alpha)**n+L*(1+c_alpha)**n)
    return f

def error_parset_hillfit(parset,n_target=1,psfunc=None,GRFfunc=None,npoints=1000,norm=False,plot=False):
    x05=1
    xvals_h=np.logspace(np.log10(0.005),np.log10(5),npoints)
    #xvals_h=np.arange(0.005,7,0.001)linspace(0.005,7,npoints)
    yvals_hill=hill(n_target,xvals_h)
    yvals_hill=np.array(yvals_hill)

    #
    xvals=xvals_h
    if norm:
        x05=psfunc(parset)[2]
        #xvals=np.logspace(np.log10(x05)-2,np.log10(x05)+2,npoints)
        values=np.array([GRFfunc(parset,x*x05) for x in xvals_h])
        if plot:
            fig,ax=plt.subplots(1,1)
            plt.plot(xvals_h,yvals_hill,color="gray",linestyle=":",label="Hill")
            plt.plot(xvals_h,values,label="f",color="k")
            plt.xlabel("x")
            plt.ylabel("fractional saturation")
            plt.legend()
            plt.show()
    else:
    
        values=np.array([GRFfunc(parset,x) for x in xvals]) #fit directly, no normalisation
    #return np.sum((yvals_hill-values)**2)
    return np.sum(np.abs(yvals_hill-values))

def error_parset_psfit(parset,psfunc=None,ptarget=1,starget=1):

    ps,stp,x05=psfunc(parset)
    if ps<0:
        return 100
    else:
    
        return np.sum(np.abs(ptarget-ps)+np.abs(starget-stp))

    

def find_parset(n_target,getinit=True,initpars=None,printminerror=0.03,nmax=10000,minerror=0.001,seed=1,psfunc=None,GRFfunc=None,min_pc=1e-3,max_pc=1e3,idx_c=None,
    min_pos=0.5,max_pos=1.2,min_stp=0.5,max_stp=1.3,min_p=1e-3,max_p=1e3,npars=27, fitto="ps",plot=False,norm=True,fout_intermediate=None,fout_final=None,beta=0.25):
    #norm=True if normalise when fitting to curve
    
    if fout_intermediate is None:
        fhi=sys.stdout
        closei=False
    else:
        fhi=open(fout_intermediate,"w")
        closei=True

    if fout_final is None:
        fhf=sys.stdout
        closef=False
    else:
        fhf=open(fout_final,"w")
        closef=True

    if fitto=="ps":
        ptarget,starget=BF.position_steepness_hill(n_target)
        print(ptarget,starget,file=fhf)
    #else fit to hill curve

    #Now find parameter set that matches this ptarget, starget


    #found_parset=False
    bestpars_list=[]
    lastprinted=2*printminerror 
    
    np.random.seed(seed)
    print(seed,file=fhf)
    if getinit is True:
        found=False
    else:
        found=True
        parset=initpars
        p,s,x05=psfunc(parset)
    if idx_c is not None:
        minplogc=np.log10(min_pc)
        maxplogc=np.log10(max_pc)
        nparsnoc=idx_c
        nparsc=npars-idx_c
    else:
        nparsnoc=npars
        nparsc=0
    minplog=np.log10(min_p)
    maxplog=np.log10(max_p)
    minpars_list=[min_p]*nparsnoc+[min_pc]*nparsc
    #    print(minlogpars_list)
    maxpars_list=[max_p]*nparsnoc+[max_pc]*nparsc
    while not found:
        parset=10**np.random.uniform(minplog,maxplog,size=nparsnoc)
        if nparsc>0:
            parset=np.concatenate((parset,10**np.random.uniform(minplogc,maxplogc,size=nparsc)))
            #print("initial parset", parset)

        p,s,x05=psfunc(parset)
        if p>0:
            if p>min_pos and p<max_pos and s>min_stp and s<max_stp:
                found=True

    if fitto=="ps":
        besterror=np.abs(p-ptarget)+np.abs(s-starget)
    else:
        besterror=error_parset_hillfit(parset,n_target=n_target,psfunc=psfunc,GRFfunc=GRFfunc,plot=plot,norm=norm)
    bestparset=parset
    error=besterror
    print("starting", besterror,file=fhf)
    n=0
    if fitto=="ps":
        error_limits=[0.25,0.01,0.005]
    else:
        error_limits=[40,20,5]

    while error>minerror and n<nmax:
        #newparset=parset.copy()
        newparset=bestparset.copy()
        parset=bestparset
        
        if error>error_limits[0]:
            minperturb=0.1
            maxperturb=10
            pmut=0.5
        #else:
        elif error>error_limits[1]:
            minperturb=0.5
            maxperturb=2
            pmut=0.5
        elif error>error_limits[2]:
            minperturb=0.8
            maxperturb=1.2
            pmut=0.25
        else:
            minperturb=0.9
            maxperturb=1.1
            pmut=0.25

        for p in range(len(parset)):
            if np.random.uniform()<pmut:
                newpar=10**np.random.uniform(np.log10(parset[p]*minperturb),np.log10(parset[p]*maxperturb))

                if newpar<minpars_list[p]:
                    newpar=minpars_list[p]
                if newpar>maxpars_list[p]:
                    newpar=maxpars_list[p]
                newparset[p]=newpar
            

        parset=newparset
        n+=1
        if fitto=="ps":
            p,s,x05=psfunc(parset)
            
            if p>0:
                error=np.abs(p-ptarget)+np.abs(s-starget)
                keep=False
        else:
            error=error_parset_hillfit(parset,n_target=n_target,psfunc=psfunc,GRFfunc=GRFfunc,plot=plot,norm=norm)
            if plot:
                print(error)
                sys.stdout.flush()

        keep=False
        if error>besterror:
            if np.random.uniform()<beta: #np.exp((besterror-error)):
                keep=True
        else:
            keep=True
        if keep:
            bestparset=parset
            besterror=error
            print(besterror, end=",",file=fhi)
            if error<printminerror:
                if error<0.9*lastprinted:
                    print("\n",file=fhi)
                    #print("#newerror:",besterror)
                    print("#newresult:",",".join(map(str,bestparset)),file=fhi)
                    lastprinted=error
            #sys.stdout.flush()
            #slightly mutate parset
    print("",file=fhi)
    #if besterror<0.1:
    #    bestpars_list.append([bestparset,besterror])
    #if besterror<0.001:
    #    found_parset=True

    if n==nmax:
        print("reached nmax",file=fhf)
    else:
        print("converged",file=fhf)

    print("#besterror:",besterror,file=fhf)
    print("#result:",",".join(map(str,bestparset)),file=fhf)

    if closef:
        fhf.close()
    if closei:
        fhi.close()
        
    return bestparset
