'''
Make control plots of data and MC templates per-year for data and MC 
'''

import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
from math import pi
from collections import OrderedDict
import array
import subprocess
from TIMBER.Analyzer import HistGroup
import sys
sys.path.insert(0,'./modules/')
from TwoDHistStyle import setTDRStyle

# Options for plotting
stack_style = {
    'edgecolor': (0, 0, 0, 0.5),
}
errorbar_style = {
    'linestyle': 'none',
    'marker': '.',      # display a dot for the datapoint
    'elinewidth': 2,    # width of the errorbar line
    'markersize': 10,   # size of the error marker
    'capsize': 0,       # size of the caps on the errorbar (0: no cap fr)
    'color': 'k',       # black 
}

# Function stolen from https://root-forum.cern.ch/t/trying-to-convert-rdf-generated-histogram-into-numpy-array/53428/3
def hist2array(hist, include_overflow=False, return_errors=False):
    '''Create a numpy array from a ROOT histogram without external tools like root_numpy.

    Args:
        hist (TH1): Input ROOT histogram
        include_overflow (bool, optional): Whether or not to include the under/overflow bins. Defaults to False. 
        return_errs (bool, optional): Whether or not to return an array containing the sum of the weights squared. Defaults to False.

    Returns:
        arr (np.ndarray): Array representing the ROOT histogram
        errors (np.ndarray): Array containing the sqrt of the sum of weights squared
    '''
    hist.BufferEmpty()
    root_arr = hist.GetArray()
    if isinstance(hist, ROOT.TH1):
        shape = (hist.GetNbinsX() + 2,)
    elif isinstance(hist, ROOT.TH2):
        shape = (hist.GetNbinsY() + 2, hist.GetNbinsX() + 2)
    elif isinstance(hist, ROOT.TH3):
        shape = (hist.GetNbinsZ() + 2, hist.GetNbinsY() + 2, hist.GetNbinsX() + 2)
    else:
        raise TypeError(f'hist must be an instance of ROOT.TH1, ROOT.TH2, or ROOT.TH3')

    # Get the array and, optionally, errors
    arr = np.ndarray(shape, dtype=np.float64, buffer=root_arr, order='C')
    if return_errors:
        errors = np.sqrt(np.ndarray(shape, dtype='f8', buffer=hist.GetSumw2().GetArray()))

    if not include_overflow:
        arr = arr[tuple([slice(1, -1) for idim in range(arr.ndim)])]
        if return_errors:
            errors = errors[tuple([slice(1, -1) for idim in range(errors.ndim)])]

    if return_errors:
        return arr, errors
    else:
        return arr


def poisson_conf_interval(k):
    """
    Calculate Poisson (Garwood) confidence intervals using ROOT's TH1 with kPoisson error option.
    
    Parameters:
    k (array): The number of counts (events) per bin.

    Returns:
    lower (array): Bin count - lower error.
    upper (array): Bin count + upper error.
    """
    lower = np.zeros_like(k, dtype=float)
    upper = np.zeros_like(k, dtype=float)
    #Temp hist to exploit ROOT's built-in CI calculating function
    hist = ROOT.TH1F("hist_delete", "", 1, 0, 1)
    hist.SetBinErrorOption(ROOT.TH1.kPoisson)
    hist.Sumw2()

    for i, count in enumerate(k):
        hist.SetBinContent(1, count)
        
        lower[i] = hist.GetBinContent(1) - hist.GetBinErrorLow(1)
        upper[i] = hist.GetBinContent(1) + hist.GetBinErrorUp(1)
        
    hist.Delete()
    
    return lower, upper    

def plot_stack(
    outname,
    data = None, # numpy array
    bkgs = {},  # {latex name : (numpy array, color)} ordered by yield - use OrderedDict
    sigs = {},  # {latex name : (numpy array, color)}
    edges = None,
    title = '',
    xtitle = '',
    subtitle = '',
    totalBkg = None,
    logyFlag = False,
    plot_ratio = False,
    lumiText = r'$138 fb^{-1} (13 TeV)$',
    extraText = 'Preliminary',
    units='GeV'):

    if not plot_ratio:
        plt.style.use([hep.style.CMS])
        fig, ax = plt.subplots()
    else:
        plt.style.use([hep.style.CMS])
        fig, (ax, rax) = plt.subplots(
            2, 1, dpi=200, figsize=(12,13), gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05}, sharex=True
        )

    # obtain background and signals
    bkg_stack = np.vstack([val[0] for key, val in bkgs.items()])
    bkg_stack = np.hstack([bkg_stack, bkg_stack[:,-1:]])
    bkg_stack = np.hstack([bkg_stack])
    bkg_colors = [val[1] for key, val in bkgs.items()]
    bkg_labels = [key for key, val in bkgs.items()]

    sig_stack = np.vstack([val[0] for key, val in sigs.items()])
    sig_stack = np.hstack([sig_stack, sig_stack[:,-1:]])
    sig_stack = np.hstack([sig_stack])
    sig_colors = [val[1] for key, val in sigs.items()]
    sig_labels = [key for key, val in sigs.items()]

    ax.stackplot(edges, bkg_stack, labels=bkg_labels, colors=bkg_colors, step='post', **stack_style)
    width = edges[1]-edges[0]
    if 'vsQCD' in outname: units = ''
    ax.set_ylabel(f'Events')
    if plot_ratio:
        rax.set_xlabel(xtitle)
    else:
        ax.set_xlabel(xtitle)
   
    # only plot data for preselection 
    #if ('preselection' in outname): 
    lower_errors, upper_errors = poisson_conf_interval(data)
    yerr = [data - lower_errors, upper_errors - data]
    bin_centers = (edges[:-1] + edges[1:])/2
    ax.errorbar(x=bin_centers, y=data, yerr=yerr, xerr=None, label='Data', **errorbar_style)

    # plot signals
    for key,val in sigs.items():
        sigarr = val[0]
        scaling = (totalBkg.max()/sigarr.max())
        ax.step(x=edges, y=np.hstack([sigarr,sigarr[-1]])*scaling, where='post', color=val[1], label=r'%s $\times$ %s'%(key,round(scaling,1)))
    
    if logyFlag:
        if ('preselection' in outname):
            if totalBkg.max() >= data.max():
                ax.set_ylim(0.01, totalBkg.max()*1e10)
            else:
                ax.set_ylim(0.01, data.max()*1e10)
        else:
            ax.set_ylim(0.01, totalBkg.max()*1e10)
        ax.set_yscale('log')
    else:
        if ('preselection' in outname):
            if totalBkg.max() >= data.max():
                ax.set_ylim(0, totalBkg.max()*1.5)
            else:
                ax.set_ylim(0, data.max()*1.5)
        else:
            ax.set_ylim(0, totalBkg.max()*1.5)

    ax.legend(frameon=True,fontsize=14)
    if ('preselection' not in outname): ax.margins(x=0)
    ax.margins(y=0)
    hep.cms.label(loc=0, ax=ax, label=extraText, rlabel='', data=True if 'preselection' in outname else False, fontsize=20)
    hep.cms.lumitext(lumiText,ax=ax, fontsize=20)

    if plot_ratio:
        yerr = np.nan_to_num(
            np.abs(
                poisson_conf_interval(data)
                - data
            )
            / (totalBkg + 1e-5)
        )
        rax.errorbar(x=bin_centers, y=data/(totalBkg+1e-5), yerr=yerr, xerr=None, **errorbar_style)
        '''
        hep.histplot(
            data / (totalBkg + 1e-5),
            yerr=yerr,
            ax=rax,
            histtype="errorbar",
            color="black",
            capsize=4,
        )
        ''' 
        rax.set_ylabel("Data/MC")
        rax.set_ylim([0,2])
        rax.grid()

    if 'MET' in outname:
        plt.savefig(outname)
    else:
        plt.savefig(outname,bbox_inches='tight')
    plt.close('all') # close all figures so that they don't consume extra memory


def rebin(inHist, newBins):
    '''
    inHist  : TH1
    newBins : array with new bins (must be a subset of old bins)

    returns : new rebinned histogram
    '''
    print(f'\tRebinning {inHist.GetName()} to {len(newBins)-1} bins')
    nbins = len(newBins)
    hOut = inHist.Rebin(nbins-1, inHist.GetName(), newBins)
    hOut.SetDirectory(0)
    return hOut


redir = 'root://cmsxrootd.fnal.gov/'
#fname = '{redir}/store/user/mhesford/XHYbbWW_semileptonic/plots/distributions/{proc}_{year}_distributions.root'
fname = 'plots/distributions/{proc}_{year}_distributions.root'
fTest = ROOT.TFile.Open(fname.format(redir=redir,proc='ttbar-allhad',year='18'),'READ')
histNames = [i.GetName() for i in fTest.GetListOfKeys()]
histTitles = [i.GetTitle() for i in fTest.GetListOfKeys()]


def fileExists(proc,year):
    try:
        f = subprocess.check_output(f'ls plots/distributions/ | grep {proc}_{year}_distributions.root',shell=True,text=True)
        return 1
    except:
        return 0


# TStyle for plotting histograms using ROOT
style = setTDRStyle()
style.cd()
# Text for HEM plots
cms_text = ROOT.TLatex(-3,2,'#splitline{#bf{CMS}}{#it{Work in Progress}}')
cms_text.SetTextSize(0.04)
cms_text.SetTextFont(42)
cms_text.SetTextColor(ROOT.kBlack)
lumi_text = ROOT.TLatex(3.14,2.5,'60 fb^{-1}, 2018 (13 TeV)')
lumi_text.SetTextSize(0.03)
lumi_text.SetTextFont(42)
lumi_text.SetTextColor(ROOT.kBlack)
lumi_text.SetTextAlign(31)

# Loop over years
for year in ['16','16APV','17','18']:
    for Rebin in [False]: #[False,True]
        for i,histName in enumerate(histNames):
            #if not 'HEM' in histName: continue
            if 'HEM' in histName:
                if year != '18': continue
                c = ROOT.TCanvas('c','c',2000,2000)
                obj = 'H_{bb}' if 'hbb' in histName else 'W_{qq}'
                procs = {
                    'Diboson' : [['WW','WZ','ZZ'],'Diboson'],
                    'ttbar' : [['ttbar-allhad','ttbar-semilep', 'ttbar-alllep'],'tt'],
                    'ST' : [['ST-t','STbar-t','ST-tW','STbar-tW','ST-slep'],'Single Top'],
                    'QCD' : [['QCDHT200','QCDHT300','QCDHT500','QCDHT700','QCDHT1000','QCDHT1500','QCDHT2000'],'QCD'],
                    'VJets' : [['WJetsHT200','WJetsHT400','WJetsHT600','WJetsHT800','WJetsLepHT200','WJetsLepHT400','WJetsLepHT600',
                                'WJetsLepHT800','WJetsLepHT1200','WJetsLepHT2500','ZJetsHT400','ZJetsHT600','ZJetsHT800'],'V+Jets'],
                    'DYJets' : [['DYJetsM4HT100','DYJetsM4HT200','DYJetsM4HT400','DYJetsM4HT600','DYJetsM50HT100','DYJetsM50HT200',
                                 'DYJetsM50HT400','DYJetsM50HT600','DYJetsM50HT800','DYJetsM50HT1200','DYJetsM50HT2500'],'Drell-Yan'],
                    'Data' : [['SingleMuonDataA','SingleMuonDataB','SingleMuonDataC','SingleMuonDataD','EGammaDataA','EGammaDataB','EGammaDataC','EGammaDataD'],'Data (2018)']
                }
                for proc,info in procs.items():
                    hem_group = HistGroup(f'{proc}')
                    for p in info[0]:
                        f = ROOT.TFile.Open(fname.format(proc=p,year=year),'READ')
                        hem = f.Get(histName)
                        hem.SetDirectory(0)
                        hem_group[p] = hem
                    hem_merged = hem_group.Merge()
                    hem_merged.GetXaxis().SetTitle(obj+' #phi')
                    hem_merged.GetYaxis().SetTitle(obj+' #eta')
                    hem_merged.SetTitle(info[1])
                    hem_merged.GetXaxis().SetRangeUser(-np.pi,np.pi)
                    hem_merged.GetYaxis().SetRangeUser(-2.4,2.4)
                    hem_merged.Draw('COLZ')
                    cms_text.Draw()
                    lumi_text.Draw()
                    c.Update()
                    c.Print('plots/distributions/png/{}.png'.format(histName+f'_{proc}'))
                    c.Clear()
                continue
                
            #if ('Wqq' not in histName and 'Higgs' not in histName) or year != '18': continue
            # plot results after kinematic preselection and after WW identification
            #if ('preselection' not in histName): continue# and ('stage1' not in histName): continue
            #if ('particleNet_mass' not in histName) and ('softdrop' not in histName) and ('vsQCD' not in histName): continue

            # plot results after H,W,W selection 
            #if ('stage2' not in histName): continue
            #if ('vsQCD' not in histName) and ('softdrop' not in histName) and ('particleNet_mass' not in histName): continue

            print(f'Plotting {histName}')

            if 'Higgs' in histName:
                obj = r'$H_{bb}$'
            elif 'Wqq' in histName:
                obj = r'$W_{qq}$'
            elif 'MET' in histName:
                obj = r'$p_{T}^{miss}$'
            elif 'Neutrino' in histName:
                obj = r'$\nu$'
            elif 'Electron' in histName:
                obj = r'e'
            elif 'Muon' in histName:
                obj = r'$\mu$'

            # Get binnings
            if ('particleNet_mass' in histName) or ('softdrop' in histName):
                edges = np.linspace(50,250,51)
                new_edges = array.array('d',np.linspace(50,250,26))
                new_edges_np = np.linspace(50,250,26)
                if 'particleNet_mass' in histName:
                    xtitle = r'%s $m_{reg}$ [GeV]'%(obj)
                else: 
                    xtitle = r'%s $m_{SD}$ [GeV]'%(obj)

            elif 'pt' in histName:
                if 'on_pt' in histName: #lepton pt:
                    edges = np.linspace(35,1200,24)
                    new_edges = array.array('d',np.linspace(35,1200,24))
                    new_edges_np = np.linspace(25,1200,13)
                elif 'MET' in histName: #MET pt
                    edges = np.linspace(0,1400,29)
                    new_edges = array.array('d',np.linspace(0,1400,15))
                    new_edges_np = np.linspace(0,1400,15)
                else: #jet pt
                    edges = np.linspace(300,2500,45)
                    new_edges = array.array('d',np.linspace(300,2500,23))
                    new_edges_np = np.linspace(300,2500,23)
                xtitle = r'%s $p_T$ [GeV]'%(obj)

            elif 'eta' in histName:
                edges = np.linspace(-3,3,41)
                new_edges = np.linspace(-3,3,21)
                new_edges_np = np.linspace(-3,3,21)
                xtitle = r'%s $\eta$'%(obj)

            elif 'phi' in histName:
                edges = np.linspace(-pi,pi,43)
                new_edges = np.linspace(-pi,pi,21)
                new_edges_np = np.linspace(-pi,pi,21)
                xtitle = r'%s $\phi$'%(obj)

            elif 'RelIso' in histName:
                edges = np.linspace(0,1,101)
                new_edges = np.linspace(0,1,51)
                new_edges_np = np.linspace(0,1,51)
                xtitle = r'%s $iso_{mini}$'%(obj)

            elif 'vsQCD' in histName:
                edges = np.linspace(0,1,51)
                new_edges = np.linspace(0,1,26)
                new_edges_np = np.linspace(0,1,26)
                if 'HbbvsQCD' in histName:
                    xtitle = r'%s $T_{Hbb}^{MD}$'%(obj)
                elif 'WqqvsQCD' in histName:
                    xtitle = r'%s $T_{Wqq}^{MD}$'%(obj)

            elif ('mhww' in histName) or ('mww' in histName):
                edges = np.linspace(0,3000,101)
                new_edges = array.array('d',np.linspace(0,3000,31))
                new_edges_np = np.linspace(0,3000,31)
                if 'mhww' in histName:
                    xtitle = r'$m_X$ [GeV]'
                else:
                    xtitle = r'$m_Y$ [GeV]'

            # get test histogram for zeros
            testHist = fTest.Get(histName)
            if Rebin: testHist = rebin(testHist,new_edges)
            testHist = hist2array(testHist)

            data = [np.zeros_like(testHist),'black']
            tt = [np.zeros_like(testHist),'lightcoral']
            st = [np.zeros_like(testHist),'sandybrown']
            vj = [np.zeros_like(testHist),'lightgreen']
            dj = [np.zeros_like(testHist),'cornflowerblue']
            qcd = [np.zeros_like(testHist),'lightyellow']
            dib = [np.zeros_like(testHist),'mediumpurple']
            ttv = [np.zeros_like(testHist),'mediumorchid']
            tth = [np.zeros_like(testHist),'violet']
            vh = [np.zeros_like(testHist),'hotpink']
            hincl = [np.zeros_like(testHist),'pink']

            xy_4000_2500 = [np.zeros_like(testHist),'black']
            xy_3000_1400 = [np.zeros_like(testHist),'red']
            xy_2000_1000 = [np.zeros_like(testHist),'blue']
            xy_1000_500 = [np.zeros_like(testHist),'green']

            # diboson
            for proc in ['WW','ZZ','WZ']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'	Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                dib[0] += a
                f.Close()

            for proc in ['ttbar-allhad','ttbar-semilep','ttbar-alllep']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'	Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                tt[0] += a
                f.Close()

            for proc in ['ST-t','STbar-t','ST-slep','ST-tW','STbar-tW']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                st[0] += a
                f.Close()

            for proc in ['WJetsHT200','WJetsHT400','WJetsHT600','WJetsHT800','WJetsLepHT200','WJetsLepHT400','WJetsLepHT600','WJetsLepHT800','WJetsLepHT1200','WJetsLepHT2500','ZJetsHT400','ZJetsHT600','ZJetsHT800']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'	Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                vj[0] += a
                f.Close()

            for proc in ['DYJetsM4HT100','DYJetsM4HT200','DYJetsM4HT400','DYJetsM4HT600','DYJetsM50HT100','DYJetsM50HT200','DYJetsM50HT400','DYJetsM50HT600','DYJetsM50HT800','DYJetsM50HT1200','DYJetsM50HT2500']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                dj[0] += a
                f.Close()

            for proc in ['XHY-4000-2500']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'	Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                xy_4000_2500[0] += a
                f.Close()

            for proc in ['XHY-3000-1400']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                xy_3000_1400[0] += a
                f.Close()

            for proc in ['XHY-2000-1000']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                xy_2000_1000[0] += a
                f.Close()

            for proc in ['XHY-1000-500']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                xy_1000_500[0] += a
                f.Close()

            for proc in ['QCDHT200','QCDHT300','QCDHT500','QCDHT700','QCDHT1000','QCDHT1500','QCDHT2000']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'	Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                qcd[0] += a
                f.Close()

            for proc in ['ttWJets','ttWJetsLep','ttZJets']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                ttv[0] += a
                f.Close()

            for proc in ['ttHtoBB','ttHnoBB']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                tth[0] += a
                f.Close()

            for proc in ['WpHbb','WmHbb','ZHbbll','ggZHbbll','WpHtoWW','WmHtoWW','ZHtoWW']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                vh[0] += a
                f.Close()

            for proc in ['ggHbb','VBFHbb','ggHWW']:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'        Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                hincl[0] += a
                f.Close()

            datasets = {
                '16': [
                    'SingleMuonDataF','SingleMuonDataG','SingleMuonDataH',
                    'SingleElectronDataF','SingleElectronDataG','SingleElectronDataH',
                    'SinglePhotonDataF','SinglePhotonDataG','SinglePhotonDataH'
                ],
                '16APV': [
                    'SingleMuonDataB1','SingleMuonDataB2','SingleMuonDataC','SingleMuonDataD','SingleMuonDataE','SingleMuonDataF',
                    'SingleElectronDataB','SingleElectronDataC','SingleElectronDataD','SingleElectronDataE','SingleElectronDataF',
                    'SinglePhotonDataB','SinglePhotonDataC','SinglePhotonDataD','SinglePhotonDataE','SinglePhotonDataF'
                ],
                '17': [
                    'SingleMuonDataB','SingleMuonDataC','SingleMuonDataD','SingleMuonDataE','SingleMuonDataF',
                    'SingleElectronDataB','SingleElectronDataC','SingleElectronDataD','SingleElectronDataE','SingleElectronDataF1',#'SingleElectronDataF2',
                    'SinglePhotonDataB','SinglePhotonDataC','SinglePhotonDataD','SinglePhotonDataE','SinglePhotonDataF'
                ],
                '18': [
                    'SingleMuonDataA','SingleMuonDataB','SingleMuonDataC','SingleMuonDataD',
                    'EGammaDataA','EGammaDataB','EGammaDataC','EGammaDataD'
                ]
            }

            for proc in datasets[year]:
                if not fileExists(proc,year):
                    print(f'\t\t{proc} {year} does not exist')
                    continue

                print(f'	Adding histogram for {proc}_{year}')
                f = ROOT.TFile.Open(fname.format(redir=redir,proc=proc,year=year),'READ')
                h = f.Get(histName)
                if Rebin: h = rebin(h,new_edges)
                a = hist2array(h)
                data[0] += a
                f.Close()

            total_withQCD = np.zeros_like(testHist)
            total_noQCD = np.zeros_like(testHist)
            for bkg in [tt,st,vj,dj,dib,ttv,tth,vh,hincl]:
                total_noQCD += bkg[0]
            for bkg in [tt,st,vj,dj,dib,qcd,ttv,tth,vh,hincl]:
                total_withQCD += bkg[0]


            lumis = {
                '16APV': r'$20 fb^{-1}$, 2016APV (13 TeV)',
                '16': r'$17 fb^{-1}$, 2016 (13 TeV)',
                '17': r'$41 fb^{-1}$, 2017 (13 TeV)',
                '18': r'$60 fb^{-1}$, 2018 (13 TeV)'
            }


            bkgHists = OrderedDict(
                [
                    (r'ttV',ttv),
                    (r'ttH',tth),
                    (r'VH',vh),
                    (r'H$\to$ bb/WW',hincl),
                    (r'Diboson',dib),
                    (r'$t\bar{t}$',tt),
                    (r'Single Top',st),
                    (r'QCD',qcd),
                    (r'V+Jets',vj),
                    (r'DY+Jets',dj)
                ]
            )
            bkgHistsNoQCD = OrderedDict(
                [
                    (r'$t\bar{t}$',tt),
                    (r'Single Top',st),
                    (r'V+Jets',vj),
                    (r'DY+Jets',dj),
                    (r'Diboson',dib)
                ]
            )
            sigHists = OrderedDict(
                [
                    (r'$X_{4000}, Y_{2500}$',xy_4000_2500),
                    (r'$X_{3000}, Y_{1400}$',xy_3000_1400),
                    (r'$X_{2000}, Y_{1000}$',xy_2000_1000),
                    (r'$X_{1000}, Y_{500}$',xy_1000_500)
                ]
            )

            # First do with linear scale and no ratio plot
            plot_stack(
                outname=f'plots/distributions/png/{histName}_{year}{"_rebin" if Rebin else ""}.png',
                data=data[0],
                bkgs=bkgHists,
                sigs=sigHists,
                totalBkg=total_withQCD,
                edges=new_edges_np if Rebin else edges,
                title=histName,
                xtitle=xtitle,
                plot_ratio=True,
                lumiText=lumis[year],
                extraText='Work in progress',
                logyFlag=False
            )

            '''
            # Now do it without QCD MC
            plot_stack(
                outname=f'plots/distributions/png/{histName}_{year}{"_rebin" if Rebin else ""}_noQCD.png',
                data=data[0],
                bkgs=bkgHistsNoQCD,
                sigs=sigHists,
                totalBkg=total_noQCD,
                edges=new_edges_np if Rebin else edges,
                title=histName,
                xtitle=xtitle,
                plot_ratio=False,
                lumiText=lumis[year],
                extraText='',
                logyFlag=False
            )


            
            # now do it with log scale and ratio plot
            bkgHists = OrderedDict(
                [
                    (r'Hbb',hbb),
                    (r'HWW',hww),
                    (r'Diboson',dib),
                    (r'ST',st),
                    (r'Z+Jets',zj),
                    (r'W+Jets',wj),
                    (r'$t\bar{t}$',tt),
                    ('QCD',qcd)
                ]
            )
            plot_stack(
                outname=f'plots/{histName}_{year}{"_rebin" if Rebin else ""}_logy.png',
                data=data[0],
                bkgs=bkgHists,
                sigs=sigHists,
                totalBkg=total_withQCD,
                edges=new_edges_np if Rebin else edges,
                title=histName,
                xtitle=xtitle,
                plot_ratio=True,
                lumiText=lumis[year],
                extraText='Work in progress',
                logyFlag=True
            )
            '''
