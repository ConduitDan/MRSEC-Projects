import glob
import matplotlib.pyplot as plt
import os,time,sys
import numpy as np
import statistics as st
from scipy.stats import norm
from subprocess import Popen,call,PIPE,STDOUT
import subprocess
import pandas as pd
from Read_AA_Data import extractdata
from scipy.stats import ks_2samp


"""
This code creates a folder and runs the MC simulations and minimizes the error function:
For efficiency, we pass log of the parameters to run simulations.

parameters :  Epsilon: elastic modulus) Kappa_theta : bending modulus) Kappa_phi : bending modulus)

"""

def initialize_folder(ROOT,METHOD,DIRNAME,seed):
    #make the folder, write the sdee file and change the permissions
    os.chdir('%s'%ROOT)
    

    #call('mkdir %s'%DIRNAME , shell=True)
    os.chdir("./%s"%DIRNAME)
    with open('seed','w') as seedfile:
        seedfile.write('%s\n'%seed)
    cwd=os.getcwd()
    print(cwd)

    with open('copysource.sh',"w") as fcopy:
        fcopy.write( "mkdir source\nscp  ../../HEsourceML/* ./source/ ;\ncd source ;\nmake ;\nsleep 5\ncd .. ;")
        #fcopy.write( "mkdir source;\nscp  ../../HEsourceML/assemble ./source/; \n")
    call("chmod +x copysource.sh", shell=True)
    call("./copysource.sh", shell=True)





def extractdata_df(datafile):
    
    df=pd.read_csv(datafile,header=0,names=['l0', 'l1', 't0', 't1', 'p33', 'p12', 'p01', 'p20'])
    df2=df.astype(np.float64)
    
    return(df2)

#### function to run simualtions
#x is epsilon0 kappa0 kappaPhi0 

def runSim(x,df_hist_AA,dirdir,y): 
    
    print("x",x)
    x00=10**x[0] #kappa_l
    x10=10**x[1] #kappa_{dihedral} (theta in code - phi in plot)
    x20=10**x[2] #kappa_{bind} (phi in code - theta in plot)

    X=[x00,x10,x20]

    #### Set up and run simulations

    with open("minimize-record.csv","a+") as fout:
        fout.write('%d,%d,%d'%(x00,x10,x20))
    files=glob.glob('data*dat')
    #dirdir=os.getcwd()
    print(dirdir)
    if len(files)>0:
        timestr=time.strftime("%Y%m%d-%H%M%S")
        
        #call("mv %s/data.dat %s/data_%s.csv"%(dirdir,dirdir,timestr), shell=True)
        #call("mv %s/energy.dat %s/energy_%s.csv"%(dirdir,dirdir,timestr), shell=True) 
        #call("rm T4* trajlammps* ", shell=True)
    #sys.exit()
    with open('seed','r') as seedfile:
        seed=int(seedfile.readline().split()[0])
    with open('seed','w') as seedfile:
        seedfile.write(str(seed+1))
           
    runstring= "./source/assemble %s . 1"%seed
    
    runparamEn="%s %s %s"%(str(x00),str(x10),str(x20))

    runparamG="%s %s %s %s %s %s %s %s"%(str(y[0]),str(y[1]),str(y[2]),str(y[3]),str(y[4]),str(y[5]),str(y[6]),str(y[7]))
    
    frames=500
    #print('%s %s %s %s'%(runstring, runparamEn, runparamG ,str(frames) ))
    with open('runCPP.sh',"w") as frun:
        frun.write('#!/bin/bash\n\n%s %s %s %s'%(runstring, runparamEn, runparamG ,str(frames)))
    call("scp ../../HE-AA-lammps1.dat ./Initial_frame.dat", shell=True)    
    call("chmod +x runCPP.sh", shell=True)
    begintime=time.strftime("%H:%M:%S")
    print("start simulation at", begintime)
    subprocess.run(["./runCPP.sh"],stdout=PIPE , stderr=STDOUT)
    endtime=time.strftime("%H:%M:%S")
    print("simulation done at",endtime )
    
   
    #### Calculate the histogram
    print(os.getcwd())
    df_CG=extractdata_df('./data.dat')
    df_CG_combined=combine(df_CG)
    df_hist_CG=gethist(df_CG_combined,df_hist_AA['bins'])
    #return()
    return(errorSigma_dfs(df_hist_CG,df_hist_AA,X))


##now run minimization
from scipy.optimize import minimize
def Run_minimize(x0,df_AA,DIRNAME,METHOD,y):
    
    print('minimization using %s method'%METHOD)
    if METHOD=='L-BFGS-B':
        res = minimize(runSim, x0,args=(df_AA,DIRNAME,y), method='L-BFGS-B',bounds=((0, 10), (0, 10), (0, 10)),options={'disp': True,'maxiter':1000,'full_output':True})
        
        print(res.x)
    elif METHOD=='nelder-mead':
        res = minimize(runSim, x0,args=(df_AA,DIRNAME,y), method='nelder-mead',options={'disp': True,'maxiter':100})
        print('error %s'%res.x)




    

def readAAdata():
    df_AA=pd.read_csv('./AAdata.csv')
    print('######################in readAAdata ')
    print(df_AA.head())
    return(df_AA)
#

####
# 1 read tha AA data
####

def choose_name_method(METHODS,i):
    #path_here=os.getcwd()
    DIRDIR="/home/farri/Data/Haganlab/testOpt/HyperoptMinimizeThree/Minimization_Runs/"
    os.chdir("%s"%DIRDIR)
    ROOT=os.getcwd()
    print(ROOT)
    #X=10**np.array(x)
    timestr=time.strftime("%Y%m%d-%H%M%S")

    #print(timestr)
    seed=timestr.split('-')[1]
    #print(seed)
    METHOD=METHODS[i]
    DIRNAME='%s_%s'%(timestr,METHOD)
    
    #DIRNAME='Minimization-%s-%s'%(METHOD,timestr)
    #DIRNAME='Minimization'
    return(ROOT,METHOD,DIRNAME,seed)

######
def kl_divergence(p, q):
    from scipy.special import kl_div
    #(x, y, out=None)
    #return np.sum(np.where(q != 0, p * np.log(p / q), 0))
    return np.sum(np.where(q != 0, kl_div(p,q), 0))


def combine(df):
    #Index(['l0', 'l1', 't0', 't1', 'p33', 'p12', 'p01', 'p20']

    Phis=np.sort(np.concatenate((df['p01'],df['p20'],df['p33'],df['p12']),axis=0))
    Ls=np.sort(np.concatenate((df['l0'],df['l1']),axis=0))
    Thetas=np.sort(np.concatenate((df['t0'],df['t1']),axis=0))
    #print(Phis[:10])
    #print(Ls[:10])
    #print(Thetas[:10])

    df_comb=pd.DataFrame(columns=['l','theta','phi'],index=[x for x in range(len(Phis))])

    df_comb['phi']=Phis
    df_comb['l']=pd.Series(Ls)
    df_comb['theta']=pd.Series(Thetas)
    


    #print('######################in combine ')
    #print(df_comb.head())
    #print(df_comb.tail())
    #sys.exit()
    return(df_comb)
    

def gethist(df_params,bins_list): #get dataframe, return histogram and bins
    print('######################in gethist ')
    #print(df_params.head())
    #print(df_params.columns) 
    df_hist=pd.DataFrame(columns=['param','bins','hist'],index=[x for x in range(len(df_params.columns))])
    #print(df_hist.columns)
    for i,param_i in enumerate(df_params.columns):
        print(i)
        print(param_i)
        df_p=df_params[param_i].dropna()
        #print(bins_list[i])
        #print(df_p.to_list())
        hist,hist_bins=np.histogram(df_p.to_list(),bins=bins_list[i])
        hist=hist/np.sum(hist)
        df_hist.iloc[i]['param']=param_i
        df_hist.iloc[i]['bins']=hist_bins
        df_hist.iloc[i]['hist']=hist
    
    df_hist.set_index('param')
    
    print(df_hist.head())
    return(df_hist)


def calc_plot_err(df_hist_CG,df_hist_AA,X): ## change theta to Phi in plot
    
    param_names=[r'$err_{l}$',r'$err_{\phi}$',r'$err_{\theta}$']
    #df has ['param', 'bins', 'hist']
    #param is index : ['l', 'theta', 'phi']

    print(df_hist_CG.columns)
    print(df_hist_AA.columns)


    plt.figure(figsize=[12,4])
    
    tstr=time.strftime("%H%M%S")
    
    err=[]
    for i,param_i in enumerate(df_hist_CG['param']): #['l', 'theta', 'phi']
        
        plt.subplot(1,3,i+1)
        df_AA=df_hist_AA[df_hist_AA['param']==param_i] 
        df_CG=df_hist_CG[df_hist_CG['param']==param_i] 

        #df_CG = df_CG[df_CG.hist != 0]
        
        #print(len(df_CG['hist'].iloc[0]))
        #print(len(df_AA['hist'].iloc[0]))

        #print(df_CG['hist'].iloc[0][:10]) 
        #print(df_CG['hist'].iloc[0][-10:])

        #print(df_AA['hist'].iloc[0][:10]) 
        #print(df_AA['hist'].iloc[0][-10:])

        assert not np.any(np.isnan(df_AA['hist'].iloc[0]))
        assert not np.any(np.isnan(df_CG['hist'].iloc[0]))

        
        assert df_CG['bins'].iloc[0].all()==df_CG['bins'].iloc[0].all()
        oldbins=df_CG['bins'].iloc[0]
        newbins=[(oldbins[i]+oldbins[i+1])/2 for i in range(len(oldbins)-1)]
        
        err.append(100*kl_divergence(df_AA['hist'].iloc[0],df_CG['hist'].iloc[0]))
        

        #plt.figure()
        
               
        plt.plot(newbins,df_AA['hist'].iloc[0])
        plt.plot(newbins,df_CG['hist'].iloc[0])
                
        #dbin=.1*(max(df_AA['bins'].iloc[0]-min(df_AA['bins'].iloc[0])))
        #plt.xlim([min(df_AA['bins'].iloc[0])-dbin,max(df_AA['bins'].iloc[0])+dbin])
        plt.ylim([0,max(df_AA['hist'].iloc[0])*1.1])
        plt.yticks([])
        plt.title('%s=%.3f'%(param_names[i],err[i]))
        
        
        #statis,pvalue=ks_2samp(AA_plot_norm, CG_plot_norm)
        
        #pval[i]=pvalue
    plt.subplots_adjust(top=0.8)
    error=np.sum(err)
    plt.suptitle(r'$\kappa_{l} = %d,\kappa_{\phi} = %d , \kappa_{\theta} = %d , error=%.3f$'%(X[0],X[1],X[2],error),y=.98)  
    plt.tight_layout()  
    plt.savefig('%s.png'%(tstr))
    #plt.show()
    
     
    plt.close('all')
    return(err,error)    
        
    


### function to calculate error of distributions gets two data frame returns error
#
def errorSigma_dfs(df_hist_CG,df_hist_AA,X):
    #from scipy.stats import ks_2samp
    #from statistics import mean
    #### update

    err,error= calc_plot_err(df_hist_CG,df_hist_AA,X)

    
    #sigma=[1,1,1,1,1,1,1,1,1]
    #w_err=[sigma[i]*err_std[i] for i in range(len(err_std))] #+ [err_mu[i] for i in range(len(err_mu))]
    #error=1-mean(pval)
    #print(err)
    #error=sum(err)*100
    with open("minimize-record.csv","a+") as fout:
        fout.write(',%.5f\n'%(error))
    print(error)
    print("##############")
    return(error)

import numpy as np
from functools import partial
from hyperopt import fmin, tpe, hp , Trials

def f(space,df_hist_AA,dir,y):
    x1 = space['x1']
    x2 = space['x2']
    x3 = space ['x3']
    return runSim([x1,x2,x3],df_hist_AA,dir,y)




def main():

    # get the AA data with the columns: ['l0', 'l1', 't0', 't1', 'p33', 'p12', 'p01', 'p20']
    df_AA=readAAdata() 

    # For 2 type optimization get average of df_AA parameters as target values:
    Y= list(np.mean(df_AA))
    #print(list(Y))
    #sys.exit()

    # combine AA data new coliumn are: ['l', 'theta', 'phi']
    df_AA_combined=combine(df_AA)

    # get histogram of combined AA data , in the format: ['param', 'bins', 'hist']
    df_hist_AA=gethist(df_AA_combined,['auto' for x in df_AA_combined.columns])
   
    #check if things are ok
    print(df_hist_AA.head())
    #for row in df_hist_AA.iterrows():
    #    plt.figure()
    #    print(row[1]['param'])
    #    oldbins=row[1]['bins']
    #    newbins=[oldbins[i]+oldbins[i+1]/2 for i in range(len(oldbins)-1)]
    #    plt.plot(newbins,row[1]['hist'])
    #plt.show()

    #sys.exit()
    # set guess value and method
    #if len(sys.argv)==1:
    #    x0=np.log10(np.array([4000.0,80.0,600])) 
    #    m=2
    #elif len(sys.argv)>3:
    #    x0=np.log10(np.array([float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])]))
    #    m=2#int(sys.argv[4])
    #else:
    #    sys.exit("usage : guess (epsilon kappa_angle kappa_bind) , method_index ['L-BFGS-B','nelder-mead'] ")
    
    m=2
    METHODS=['L-BFGS-B','nelder-mead','hyperopt']    
    #print(x0)
    print(METHODS[m])

    
    #set the directory to run
    ROOT,METHOD,DIRNAME,seed=choose_name_method(METHODS,m)

    #### these to run minimize

    os.mkdir(DIRNAME)

    #copy this so we know which version we ran
    #call("scp ../AA-CG-Optimization-l-theta-phi-NoSeparation.py %s"%DIRNAME, shell=True)

    #initialize the folder
    initialize_folder(ROOT,METHOD,DIRNAME,seed) #copies the source code and compiles it in the directory

    #run Minimization
    #Run_minimize(x0,df_hist_AA,DIRNAME,METHOD,Y)

    space = {
    'x1': hp.uniform('x1', 2,3.5),
    'x2': hp.uniform('x2',  1, 3),
    'x3': hp.uniform('x3',  1, 3.2),
    }

    tpe_trials=Trials();

    fmin_objective = partial(f, 
                            df_hist_AA=df_hist_AA,
                            dir=DIRNAME,
                            y=Y
                            )

    #bestParams = fmin(fn = fmin_objective ,space = space)

    best = fmin(
        fn=fmin_objective,# f,
        space=space,
        algo=tpe.suggest,
        max_evals=100,
        trials=tpe_trials   
    )

    print("Found minimum after 1000 trials:")
    print(best)



#TEST runsim
def testRunSim(df_AA,DIRNAME):
    
    cwd=os.getcwd()
    print(cwd)
    x0=[3.700,8.0,1.000]
    begintime=time.strftime("%H:%M:%S")
    print(begintime)
    time.sleep(2)
    print(runSim(x0,df_AA,DIRNAME))
    endtime=time.strftime("%H:%M:%S")
    #tdelta = etime.strptime(s1, FMT)
    print(endtime )

## test errorfunction:

def c_main():
    
    

    X=np.array([0,0,0]) 
    # get the AA data with the columns: ['l0', 'l1', 't0', 't1', 'p33', 'p12', 'p01', 'p20']
    df_AA=readAAdata() 

    # For 2 type optimization get average of df_AA parameters as target values:
    Y= list(np.mean(df_AA))
    #print(list(Y))
    #sys.exit()

    # combine AA data new coliumn are: ['l', 'theta', 'phi']
    df_AA_combined=combine(df_AA)

    # get histogram of combined AA data , in the format: ['param', 'bins', 'hist']
    df_hist_AA=gethist(df_AA_combined,['auto' for x in df_AA_combined.columns])
    #print(os.getcwd())
    #sys.exit()
    y=[1,1,1]

    DIRNAME="/home/farri/Data/Haganlab/CPROJECTS/HEFOURTYPE/NEW_MLRUN/Final_Minimization/HyperoptMinimizeThree/Minimization_Runs/20210929-105255_hyperopt"

    runSim(X,df_hist_AA,DIRNAME,y)


    df_CG=extractdata_df("/home/farri/Data/Haganlab/CPROJECTS/HEFOURTYPE/NEW_MLRUN/Final_Minimization/HyperoptMinimizeThree/Minimization_Runs/20210929-105255_hyperopt/data.dat")
    print(df_CG.head())
    #sys.exit()
    df_CG_combined=combine(df_CG)
    print(df_CG_combined.head())
    print(df_hist_AA['bins'].to_list())
    #sys.exit()
    df_hist_CG=gethist(df_CG_combined,df_hist_AA['bins'].to_list())
    print(df_hist_CG.head())

    print(errorSigma_dfs(df_hist_CG,df_hist_AA,X))
"""
###testplot
#print(df_AA.head())
AA_norm= df_AA.apply(norm.fit,axis=0).to_numpy().T
#print(AA_norm)

#print(df_AA.head())
#print(AA_norm[0])


#sys.exit()
########## test sigma dist
for i,param in enumerate(df_AA.columns):
    
    #print(df_AA[param])
    param_AA,bins_AA=np.histogram(df_AA[param],bins='auto')
    
    #print(AA_norm[i])
    AA_plot_norm=norm.pdf(bins_AA, AA_norm[i][0], AA_norm[i][1])
    AA_plot_norm_test= norm.pdf(bins_AA, AA_norm[i][0], .015)  
    plt.figure()
               
        
    #plt.hist(df_AA[param],bins=bins_AA, edgecolor='gray',color='c', alpha=0.65,density=1)
    plt.plot(bins_AA,AA_plot_norm,label='AA_plot_norm')
    plt.plot(bins_AA,AA_plot_norm_test,label='AA_plot_norm_test')
    plt.legend()
    plt.show()
############
sys.exit()  
"""    
    



#### these to run test
'''
os.mkdir(DIRNAME)
initialize_folder(ROOT,METHOD,DIRNAME,seed)
with open("./minimize-record.csv","w+") as fout:
    testRunSim(df_AA,DIRNAME)

'''
#### these to test already available folder
'''
DIRNAME='Minimization-nelder-mead-20210125-172803'
with open("./test-minimize-record.csv","w+") as fout:
    print(DIRNAME)
    df_CG=extractdata('%s'%DIRNAME)
    print(errorSigma_dfs(df_CG,df_AA))

'''


if __name__ == "__main__":
    # execute only if run as a script
    main()


Activity

    2 collaborators uploaded v1 - 2

Write a comment(optional)Use the @ symbol to mention users and use the up and down arrow keys to scroll through autocomplete suggestions.
Write a comment

@mention users to notify them.
Drop files on this page to upload them into this folder.
