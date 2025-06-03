'''
Gather output files from selection script and reorganize them
Thank you for looking at this Amitav :)
'''
from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
from TIMBER.Analyzer import HistGroup
import ROOT
from collections import OrderedDict
from make_QCD import Rebin2D

def combine_ttbar(redirector,eos_path):
    path = redirector+eos_path
    #This function combines semileptonic/hadronic ttbar - so for each year and for each of the JEC's I'm hadding the corresponding semileptonic/hadronic files
    variations = ['','_JES_up','_JES_down','_JER_up','_JER_down','_JMS_up','_JMS_down','_JMR_up','_JMR_down','_UE_up','_UE_down']
    for var in variations:
       for year in ['16','16APV','17','18']:
           
           semilep = '{}XHYbbWWselection_ttbar-semilep_{}{}.root'.format(path,year,var)
           allhad = '{}XHYbbWWselection_ttbar-allhad_{}{}.root'.format(path,year,var)
           alllep = '{}XHYbbWWselection_ttbar-alllep_{}{}.root'.format(path,year,var)
           
           try:
               ExecuteCmd(f'eosrm {eos_path}XHYbbWWselection_ttbar_{year}{var}.root')
           except:
               print(f"file {path}XHYbbWWselection_ttbar_{year}{var}.root doesn't exist: can't rm")
               
           ExecuteCmd('hadd -f {}XHYbbWWselection_ttbar_{}{}.root {} {} {}'.format(path,year,var,semilep,allhad,alllep))
           
           for background in ['WJets','ZJets','DYJets','ST']:
               try:
                   ExecuteCmd(f'eosrm {eos_path}XHYbbWWselection_{background}_{year}{var}.root')
               except:
                   print(f"file {path}XHYbbWWselection_{background}_{year}{var}.root doesn't exist: can't rm")

               cmd = 'eos {} ls {}XHYbbWWselection_{}*{}{}.root'.format(redirector,eos_path,background,year,var)
               rawBkg = subprocess.check_output(cmd, shell=True).decode()
               bkgStr = ''
               for f in rawBkg.split('\n'):
                   if f == '' or 'ST-slep' in f:
                       continue
                   bkgStr += redirector+eos_path+f+' '
               ExecuteCmd('hadd -f {0}XHYbbWWselection_{1}_{2}{3}.root {4}'.format(path,background,year,var,bkgStr))
           ExecuteCmd(f'hadd -f {path}XHYbbWWselection_VJets_{year}{var}.root {path}XHYbbWWselection_WJets_{year}{var}.root {path}XHYbbWWselection_ZJets_{year}{var}.root')

def adjust_negative_bins(hist):
    #currently have some negative bins in some VV MC histos - temporary fix is just to set them to a slightly positive value    
    #NOTE: not currently considering the VV MC in the fits, but this function still gets run (in theory its not doing anything)
    print('Adjusting negative bins for histogram {}'.format(hist.GetName()))
    for x in range(0,hist.GetNbinsX()):
        for y in range(0,hist.GetNbinsY()):
            if hist.GetBinContent(x,y) < 0.0:
                hist.SetBinContent(x,y,0.002)

def organize_bkg_hists(backgrounds,binsX,binsY):
    #makes fully hadded histograms for each systematic+year combo as well as a fully nominal histogram, stored in a new rootfile for each setname (format: XHYbbWWselection_<SETNAME>.root)
    #EX: hist 'MXvMY_pass_SR__pileup_16_up' is the sum of the pileup up histo from 2016 + nominal histogram from all other years
    for setname in backgrounds:
        outFile = ROOT.TFile('selection/XHYbbWWselection_{}.root'.format(setname),'RECREATE')
        cutflow = HistGroup('cutflow_'+setname) # histgroup for cutflow tables
        for year in ['16APV','16','17','18']:
            inFile = ROOT.TFile.Open('selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(setname,year),'READ')
            print('Retreiving cutflow for {} and year 20{}'.format(setname,year))
            inHist = inFile.Get('cutflow')
            inHist.SetDirectory(0)
            cutflow.Add('cutflow_'+year,inHist) #need to add the year so the histograms don't overwrite each other
        merged_cutflow = cutflow.Merge() #merged cutflow histogram
        outFile.cd()
        merged_cutflow.Write()
    
        for region in ['pass_SR','fail_SR','pass_ttCR','fail_ttCR']:
            nom_histgroup = HistGroup('MXvMY_{}__nominal'.format(region)) #In this histgroup I will add all the nominal histograms for the different years for a given process (ttbar), then use the Merge() method to combine them into one histogram at the end
            for year in ['16APV','16','17','18']:
                #first get the nominal histogram from the target year and store separately
                inFile = ROOT.TFile.Open('selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(setname,year),'READ')
                print('Retreiving nominal histogram for {} region {} and year 20{}'.format(setname,region,year))
                inHist = inFile.Get('MXvMY_{}__nominal'.format(region)) #grab the nominal histogram from the nominal file
                inHist.SetDirectory(0)
                nom_histgroup.Add(inHist.GetName()+'_'+year,inHist) #need to add the year so the histograms don't overwrite each other

            merged_nom_hist = nom_histgroup.Merge() #after looping through all the years and grabbing all the nominal histograms, can merge the nominal histgroup
            print(merged_nom_hist.Integral())
            nom_hist_out = Rebin2D(merged_nom_hist,binsX,binsY)
            print(nom_hist_out.Integral())
            print(type(nom_hist_out))
            outFile.cd()
            nom_hist_out.Write() #write the nominal histogram to the out file

def IsoHists(setnames):
    for name in setnames:
        print('Merging isolation histograms for {}'.format(name))
        groups = OrderedDict([
            ('Wtag_SR',HistGroup('Wtag_SR')),
            ('Wtag_ttCR',HistGroup('Wtag_ttCR'))
        ]) 

        for year in ['16APV','16','17','18']:
            print('Grabbing {} histograms'.format(year))
            infile = ROOT.TFile.Open('selection/raw_selection/XHYbbWWselection_{}_{}.root'.format(name,year))
            for region in ['SR','ttCR']:
                tagHist = infile.Get('Wtag_{}'.format(region))
                tagHist.SetDirectory(0)

                groups['Wtag_{}'.format(region)].Add(year,tagHist)

            infile.Close()
        
        outfile = ROOT.TFile.Open('selection/XHYbbWWselection_{}.root'.format(name),'UPDATE')
        outfile.cd()
        for g in groups:
            h_merged = groups[g].Merge()
            h_merged.Write()
        outfile.Close()

def MakeRun2(binsX,binsY):
    #add all the data histograms together
    
    cmd = 'hadd -f selection/XHYbbWWselection_DataRun2_originalBins.root selection/raw_selection/XHYbbWWselection_*Data*'
    ExecuteCmd(cmd)
    data_in = ROOT.TFile.Open('selection/XHYbbWWselection_DataRun2_originalBins.root','READ')
    data_out = ROOT.TFile.Open('selection/XHYbbWWselection_DataRun2.root','RECREATE')
    for hist_name in data_in.GetListOfKeys(): #Rebin the histograms
        inHist = data_in.Get(hist_name.GetName())
        if 'MXvMY' in hist_name.GetName():
            print('Rebin {}'.format(hist_name.GetName()))
            outHist = Rebin2D(inHist,binsX,binsY)
        else:
            outHist = inHist
        data_out.cd()
        outHist.Write()

if __name__ == '__main__':
    redirector = 'root://cmseos.fnal.gov/'
    eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/selection_NEW/'
    cmd = 'eos {} ls {}XHYbbWWselection*'.format(redirector,eos_path)
    
    combine_ttbar(redirector,eos_path)
     
    rawFiles = subprocess.check_output(cmd, shell=True).decode()
    files = rawFiles.split('\n')
    signals = ['XHY-4000-2500','XHY-3500-1800','XHY-3000-1400','XHY-2000-1000','XHY-1600-700','XHY-1000-500']
    for fName in files:
        if fName == '':
            continue
        setname = fName.split('_')[1]
        if 'XHY' in setname and setname not in signals:
            continue
        #Throwing all the output histograms from the selection script into a directory called "raw_selection"
        ExecuteCmd('xrdcp {}{}{} selection/raw_selection/'.format(redirector, eos_path, fName)) 
    
    binsX = [600,1000,1200,1400,2000,4500] #Binning for 2DAlphabet fit, must rebin here so that I can combine with the SR QCD histograms to generate toy data
    binsY = [100,400,600,750,1000,2500]
    MakeRun2(binsX,binsY)
    
    binsX = [600,1000,1200,1400,2000,4500] #Binning for 2DAlphabet fit, must rebin here so that I can combine with the SR QCD histograms to generate toy data
    binsY = [100,400,600,750,1000,2500]
    # 'XHY-4000-2500','XHY-3500-1800','XHY-3000-1400','XHY-2000-1000','XHY-1000-500'
    organize_bkg_hists(['ttbar','DYJets','WJets','ZJets','ST','VJets'],binsX,binsY) 
    IsoHists(setnames = ['ttbar','DYJets','WJets','ZJets','ST','VJets'])
