#!/usr/bin/env python
#Update 21/10/19 added variance explained C.I.
#Order of heatmap
from sys import exit
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import argparse
    import seaborn as sns
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    print('Python module does not seem to be loaded')
    exit()

####################Arg parsing #########################
parser = argparse.ArgumentParser(description="input positional parameters: PRS file, phenotypes file (includes covars and qcovars) and variable description file. Check --help for required IDs")
parser.add_argument('jobName',help="Base name of your prs (before R2 and Pval)")
parser.add_argument('-plt','--plotResults',nargs='?', const='heat',default=None,help="Whether to plot the results as a heatmap (default) or if -plt bars is used it prints a barplots of variance explained for each phenotype")
parser.add_argument('-nvar','--numberOfVariables',default=1,help="Number of real non correlated variables to correct for mulitple testing")
parser.add_argument('-fontsize','--fontsize',default=12,help="X and Y axis fontsize")
parser.add_argument('-transpose','--TransposePlot',help="Whether to transpose",action='store_true')
parser.add_argument('-R2confint','--R2CONFINT',help="Whether to plot R2confidence intervals",action='store_true')
parser.add_argument('-heatCol','--heatMapColor',help="Heatmap color from white to this (in hex) e.g. \'#a92e4a\' (use quotes)",default='#a92e4a')
parser.add_argument('-title','--AddTitle',help="Which title to add to the heatMap")
parser.add_argument('-kwrgs','--kwargs',help="string with a dictionary of kwargs for the plot (read python and seaborn docs on kwargs)",default="{}")
parser.add_argument('-noCbarTxt','--NoColorBarText',help="Use this option to remove the colorbar text (e.g. if several hetmaps will be used).",action='store_true')
parser.add_argument('-sigoffset','--SignificanceOffset',help="Use this option to move the significance mark (*) default = 2",default=2)
parser.add_argument('-vmin','--minValue',help="Use this option to move the significance mark (*) default = 2",default=None)
parser.add_argument('-vmax','--maxValue',help="Use this option to move the significance mark (*) default = 2",default=None)
parser.add_argument('-order','--colORDER',nargs='+',help='Use this flag to set the order of the columns for the joint barplot or heatmap. Just write the headers separated by a space. It is best to use this flag at the end as it receives an undefinite number of arguments',default=None)
parser.add_argument('-figextension','--FigureExtension',help="Use this option to save the figure as pdf or png defualt pdf, do not add a '.' ",default='pdf')




parser._actions[0].help='Print this help message and exit'
args = parser.parse_args()
kwargs=eval(args.kwargs)
################### Opening Files #######################
jobName=args.jobName
PRS_R2=pd.read_csv(jobName+'R2.tab',sep='\t',header=0,index_col=0)
PRS_PVAL=pd.read_csv(jobName+'Pval.tab',sep='\t',header=0,index_col=0)
PRS_PVAL=PRS_PVAL.apply(pd.to_numeric)
PRS_R2=PRS_R2.apply(pd.to_numeric)

PRS_R2=PRS_R2.sort_index()
if args.colORDER:
    PRS_R2=PRS_R2.loc[:,args.colORDER]

print ("Read both files")
############ Make plots #################
threshold=0.05/int(args.numberOfVariables)
print ("Using %s variables to correct"%(args.numberOfVariables))

#Generating plots
plotType=args.plotResults

if plotType=='bars':
    print ("Plotting bars")
    for column in PRS_R2:
        ########### Undocumment the lines below and change fz accordingly to make your plot look how you want ############ 
        #fz=args.fontsize
        #rc={'font.size': fz, 'axes.labelsize': fz, 'legend.fontsize': fz, 
        #'axes.titlesize': fz, 'xtick.labelsize': fz, 'ytick.labelsize': fz} #change the fz to make your plot nicer if it is not looking properly
        #sns.set(rc)
        ###################################################################################################################
        fig,ax=plt.subplots(1)
        sns.set_style('white')
        sns.despine()
        bp=sns.barplot(x=PRS_R2.index,y=PRS_R2.loc[:,column]*100,ci=None,color='lightsteelblue',ax=ax)
        ax.set_title(column)
        xticks=ax.get_xticks()
        for ix,val in enumerate(PRS_PVAL.index):
            if PRS_PVAL.loc[val,column]<threshold:
                ax.text(xticks[ix]-0.35,PRS_R2.loc[val,column]*100,"%.1e"%(PRS_PVAL.loc[val,column]),fontsize=10,fontweight='bold')
            else:
                ax.text(xticks[ix]-0.35,PRS_R2.loc[val,column]*100,"%.1e"%(PRS_PVAL.loc[val,column]),fontsize=10)     
        ax.set_ylabel('Var explained (%)')
        bot,top=ax.get_ylim()
        ax.set_ylim(bot,top+0.05)
        ax.set(**kwargs)
        fig.tight_layout()
        fig.savefig('barplot%s.%s'%(jobName+column,args.FigureExtension),bbox_inches='tight',dpi=300)
elif plotType=='barsjoint':
    if args.R2CONFINT:
        PRS_R2Low=pd.read_csv(jobName+'R2Low.tab',sep='\t',header=0,index_col=0)
        PRS_R2Upp=pd.read_csv(jobName+'R2Upp.tab',sep='\t',header=0,index_col=0)
        d={'R2':PRS_R2,'R2Low':PRS_R2Low,'R2Upp':PRS_R2Upp,'Pval':PRS_PVAL}
    else:
        d={'R2':PRS_R2,'Pval':PRS_PVAL}
    DF=pd.concat(d.values(),axis=1,keys=d.keys())
    DF['Index']=DF.index
    r2melted=pd.melt(DF,id_vars='Index')
    sns.set_style("whitegrid")
    cmapa=LinearSegmentedColormap.from_list('QIMR',colors=['#f4f3f3',args.heatMapColor],N=100)
    fig,ax=plt.subplots(1)
    r2Only=r2melted.loc[r2melted.variable_0=='R2']
    sns.barplot(x=r2Only.variable_1,y=r2Only.value*100,ci=None,palette="Blues",hue=r2Only.Index,ax=ax,order=None)
    ax.get_legend().remove()
    ax.set_ylabel("Variance explained (%)", fontsize=15)
    xticks=ax.get_xticks()
    xticklabels=[i.get_text() for i in ax.get_xticklabels()]
    pvalsLog=-np.log10(r2melted.loc[r2melted.variable_0=='Pval','value'])
    norm = plt.Normalize(pvalsLog.min(), pvalsLog.max())
    sm = plt.cm.ScalarMappable(cmap=cmapa, norm=norm)
    sm.set_array([])
    for ix,patch in enumerate(ax.patches):
        currR2=patch.get_height()
        if np.isnan(currR2):
            patch.set_facecolor(sm.to_rgba(np.nan))
        else:
            currDFpos=DF.loc[:,'R2']*100==currR2
            currPval=-np.log10([i for i in DF.loc[:,'Pval'][currDFpos].values.flatten() if not np.isnan(i)][0])
            currxpos=patch.get_x()
            if args.R2CONFINT:
                currLow=[i for i in DF.loc[:,'R2Low'][currDFpos].values.flatten() if not np.isnan(i)][0]
                currHigh=[i for i in DF.loc[:,'R2Upp'][currDFpos].values.flatten() if not np.isnan(i)][0]
                ax.vlines(currxpos+patch.get_width()/2,currLow*100,currHigh*100)
            patch.set_facecolor(sm.to_rgba(currPval))
            if currPval>-np.log10(threshold):
                x=patch.get_x()
                y=currR2
                ax.text(x+patch.get_width()/5,y+0.0001,'*',fontsize=15,fontweight='bold')
    ax.set_ylabel('Variance explained\n$R^2$',fontsize=16)
    ax.set_xlabel('')
    cb=ax.figure.colorbar(sm)
    cb.set_label('-log10(pvalue)',fontsize=14)
    for label in ax.get_xticklabels():
        label.set_fontsize(12)
        label.set_weight('bold')
    for label in ax.get_yticklabels():
        label.set_fontsize(12)
        label.set_weight('bold')

    ax.set_xlabel("")
    ax.set(**kwargs)
    fig.tight_layout()
    fig.savefig('barplotJoint%s.%s'%(jobName,args.FigureExtension),dpi=300)
else:
    print ("Plotting heatmap")
    fz=args.fontsize
    rc={'font.size': fz, 'axes.labelsize': fz, 'legend.fontsize': fz, 
    'axes.titlesize': fz, 'xtick.labelsize': fz, 'ytick.labelsize': fz} #change the fz to make your plot nicer if it is not looking properly
    sns.set(rc)
    tmp=PRS_PVAL<threshold
    tmp=tmp.replace(True,'*').replace(False,'')
    tmp2=PRS_PVAL<0.05
    tmp3=PRS_PVAL<threshold
    tmp[tmp2]='*'
    tmp[tmp3]='**'
 
    fig,ax = plt.subplots(1)
    cmapa=LinearSegmentedColormap.from_list('QIMR',colors=['#f4f3f3',args.heatMapColor],N=100)
    if args.NoColorBarText:
        cbarkws=None
    else:
        cbarkws={'label':'Var explained (%)\n* p<0.05'}
    if args.TransposePlot:
        sns.heatmap(PRS_R2.transpose()*100,vmin=args.minValue,vmax=args.maxValue,cmap=cmapa,ax=ax,cbar_kws=cbarkws,
            annot=tmp.transpose(),annot_kws={'fontsize':fz,'weight':'bold'},fmt = '',xticklabels=True, yticklabels=True)
    else:
        sns.heatmap(PRS_R2*100,cmap=cmapa,ax=ax,cbar_kws=cbarkws,
            annot=tmp,annot_kws={'fontsize':fz,'weight':'bold'},fmt = '',vmin=args.minValue,vmax=args.maxValue,xticklabels=True, yticklabels=True)
    fig.tight_layout()
    #fig.set_size_inches(20,20) #uncomment this line and change the size if your heatmap is looking ugly
    if args.AddTitle:
        ax.set_title(args.AddTitle)
    ax.set(**kwargs)
    fig.savefig('HeatMap%s.%s'%(jobName,args.FigureExtension),bbox_inches='tight',dpi=300)
        
# End of script
        
        
        
        