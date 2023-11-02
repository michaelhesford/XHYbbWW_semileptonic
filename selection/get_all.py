'''
Gathers all rootfiles from the EOS and dumps it into the selection/ directory 
'''
from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
from TIMBER.Analyzer import HistGroup
import ROOT

def combine_ttbar():
    variations = ['','_JES_up','_JES_down','_JER_up','_JER_down','_JMS_up','_JMS_down','_JMR_up','_JMR_down']
    for var in variations:
       for year in ['16','16APV','17','18']:
           semilep = 'selection/raw_selection/XHYbbWWselection_ttbar-semilep_{}{}.root'.format(year,var)
           allhad = 'selection/raw_selection/XHYbbWWselection_ttbar-allhad_{}{}.root'.format(year,var) if not year == '17' else ''
           ExecuteCmd('hadd -f selection/raw_selection/XHYbbWWselection_ttbar_{}{}.root {} {}'.format(year,var,semilep,allhad))

def adjust_negative_bins(hist):
    #currently have some negative bins in some VV MC histos - temporary fix is just to zero them out, or (for better systematics) set them to a slightly positive value
    print('Adjusting negative bins for histogram {}'.format(hist.GetName()))
    for x in range(0,hist.GetNbinsX()):
        for y in range(0,hist.GetNbinsY()):
            if hist.GetBinContent(x,y) < 0.0:
                hist.SetBinContent(x,y,0.002)

def organize_bkg_hists(backgrounds):
    #makes fully hadded histograms for each systematic+year combo as well as a fully nominal histogram, stored in a new rootfile for each setname (format: XHYbbWWselection_<SETNAME>.root)
    #EX: hist 'MXvMY_SR_pass__pileup_16_up' is the sum of the pileup up histo from 2016 + nominal histogram from all other years
    for setname in backgrounds:
        outFile = ROOT.TFile('selection/XHYbbWWselection_{}.root'.format(setname),'UPDATE')
        for region in ['pass_SR','fail_SR','pass_ttCR','fail_ttCR']:
            systs = {
                '16APV' : ['L1PreFiringWeight_up','L1PreFiringWeight_down','Pdfweight_up','Pdfweight_down','Pileup_up','Pileup_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up','JMR_down'],
                '16' : ['L1PreFiringWeight_up','L1PreFiringWeight_down','Pdfweight_up','Pdfweight_down','Pileup_up','Pileup_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up','JMR_down'],
                '17' : ['L1PreFiringWeight_up','L1PreFiringWeight_down','Pdfweight_up','Pdfweight_down','Pileup_up','Pileup_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up','JMR_down'],
                '18' : ['Pileup_up','Pileup_down','Pdfweight_up','Pdfweight_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up','JMR_down']
            }
            '''
            if 'ttbar' in setname:
                systs['16APV'].append('TptReweight_up')
                systs['16APV'].append('TptReweight_down')
                systs['16'].append('TptReweight_up')
                systs['16'].append('TptReweight_down')
                systs['17'].append('TptReweight_up')
                systs['17'].append('TptReweight_down')
                systs['18'].append('TptReweight_up')
                systs['18'].append('TptReweight_down')
            '''
            nom_histgroup = HistGroup('MXvMY_{}__nominal'.format(region))
            for year in systs.keys():
                nom_years = systs.keys()
                nom_years.remove(year) #will retrieve nominal histogram for all years except target year
                '''
                if setname == 'ttbar-allhad':
                    #currently missing ttbar-allhad_17
                    if year == '17':
                        continue
                    else:
                        nom_years.remove('17')
                '''           
                #first get the nominal histogram and store separately
                inFile = ROOT.TFile.Open('selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(setname,year),'READ')
                print('Retreiving nominal histogram for {} region {} and year 20{}'.format(setname,region,year))
                inHist = inFile.Get('MXvMY_{}__nominal'.format(region))
                inHist.SetDirectory(0)
                nom_histgroup.Add(inHist.GetName()+'_'+year,inHist)

                #now make histogram for each systematic and the given year
                for sys in systs[year]:
                    print('Making {} histogram for {} variation {} and year 20{}'.format(region,setname,sys,year))
                    sys_histgroup = HistGroup('MXvMY_{}__{}_{}_{}'.format(region,sys.split('_')[0],year,sys.split('_')[1])) #Ex: MXvMY_fail_SR__Pileup_17_up
                    #first get the histogram with the target systematic variation
                    if not sys.startswith('J'): 
                        sys_path = 'selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(setname,year)
                        sys_hist = 'MXvMY_{}__{}'.format(region,sys)
                    else: #one of the JME variations, have to retrieve differently
                        sys_path = 'selection/raw_selection/XHYbbWWselection_{}_{}_{}.root'.format(setname,year,sys)
                        sys_hist = 'MXvMY_{}__nominal'.format(region)
                    sysFile = ROOT.TFile.Open(sys_path,'READ')
                    sysHist = sysFile.Get(sys_hist)
                    sysHist.SetDirectory(0)
                    sys_histgroup.Add(sysHist.GetName()+'_'+year,sysHist)
                    #next add to the histgroup the corresponding nominal histograms from all other years
                    for nom_year in nom_years:
                        nomFile = ROOT.TFile.Open('selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(setname,nom_year),'READ')
                        nomHist = nomFile.Get('MXvMY_{}__nominal'.format(region))
                        nomHist.SetDirectory(0)
                        sys_histgroup.Add(nomHist.GetName()+'_'+nom_year,nomHist)
                    #print(sys_histgroup.keys())
                    merged_hist = sys_histgroup.Merge()
                    adjust_negative_bins(merged_hist)
                    outFile.cd()
                    merged_hist.Write() #write each systematic histogram to the out file 
            if 'ttbar' in setname:
                #Apply the top pt reweight as a single correction for all the years
                tpt_systs = ['TptReweight_up','TptReweight_down']
                for tpt in tpt_systs:
                    tpt_histgroup = HistGroup('MXvMY_{}__{}'.format(region,tpt)) #Ex: MXvMY_fail_SR__TptReweight_up                
                    for year in systs.keys():
                        tpt_path = 'selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(setname,year)
                        tpt_hist = 'MXvMY_{}__{}'.format(region,tpt)
                        tptFile = ROOT.TFile.Open(tpt_path,'READ')
                        tptHist = tptFile.Get(tpt_hist)
                        tptHist.SetDirectory(0)
                        tpt_histgroup.Add(tptHist.GetName()+'_'+year,tptHist)
                    tpt_hist_merged = tpt_histgroup.Merge()
                    adjust_negative_bins(tpt_hist_merged)
                    outFile.cd()
                    tpt_hist_merged.Write() #write the tpt histogram to the out file
            merged_nom_hist = nom_histgroup.Merge()
            adjust_negative_bins(merged_nom_hist)
            outFile.cd()
            merged_nom_hist.Write() #write the nominal histogram to the out file

def MakeRun2():
    cmd = 'hadd -f selection/XHYbbWWselection_DataRun2.root selection/raw_selection/XHYbbWWselection_Data*'
    ExecuteCmd(cmd)

if __name__ == '__main__':
    
    redirector = 'root://cmseos.fnal.gov/'
    #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/selection/'
    #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/selection_HF/'
    #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/selection_dPhi/'
    eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/selection_combined/'
 
    cmd = 'eos {} ls {}XHYbbWWselection*'.format(redirector,eos_path)
    
    rawFiles = subprocess.check_output(cmd, shell=True)
    files = rawFiles.split('\n')
    for fName in files:
        ExecuteCmd('xrdcp {}{}{} selection/raw_selection/'.format(redirector, eos_path, fName))

    combine_ttbar()
        
    MakeRun2()

    #organize_bkg_hists(['XHY-4000-2500','XHY-3500,1800','XHY-3000-1400','XHY-2000-1000','XHY-1000-500','ttbar'])
    organize_bkg_hists(['XHY-3000-1400','ttbar'])
