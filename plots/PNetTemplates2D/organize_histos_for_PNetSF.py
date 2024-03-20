'''
For calculating particleNet scale factors in ttbar
Organizes 2D template histograms for use with JMAR SF tool
'''
from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
from TIMBER.Analyzer import HistGroup
import ROOT

def combine_ttbar(tagger):
    #This function combines semileptonic/hadronic ttbar - so for each year and for each of the JEC's I'm hadding the corresponding semileptonic/hadronic files
    variations = ['','_JES_up','_JES_down','_JER_up','_JER_down','_JMS_up','_JMS_down','_JMR_up','_JMR_down']
    for var in variations:
       for year in ['16','16APV','17','18']:
           semilep = 'plots/PNetTemplates2D/PNetTemplates2D_{}_ttbar-semilep_{}{}.root'.format(tagger,year,var)
           allhad = 'plots/PNetTemplates2D/PNetTemplates2D_{}_ttbar-allhad_{}{}.root'.format(tagger,year,var) if not year == '17' else ''
           ExecuteCmd('hadd -f plots/PNetTemplates2D/PNetTemplates2D_{}_ttbar_{}{}.root {} {}'.format(tagger,year,var,semilep,allhad))

def organize_ttbar_templates(tagger, wp = '0.98'):
    #tagger = 'Wqq' or 'Hbb'
    systs = { #lists of all the systematics applied for the given years, minus the top pt reweighting (deal with this separately later)
        '16APV' : ['nominal','L1PreFiringWeight_up','L1PreFiringWeight_down','Pdfweight_up','Pdfweight_down','Pileup_up','Pileup_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down',
		   'JMR_up','JMR_down','TptReweight_up','TptReweight_down'],
	'16' : ['nominal','L1PreFiringWeight_up','L1PreFiringWeight_down','Pdfweight_up','Pdfweight_down','Pileup_up','Pileup_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up',
                'JMR_down','TptReweight_up','TptReweight_down'],
	'17' : ['nominal','L1PreFiringWeight_up','L1PreFiringWeight_down','Pdfweight_up','Pdfweight_down','Pileup_up','Pileup_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up',
                'JMR_down', 'TptReweight_up','TptReweight_down'],
	'18' : ['nominal','Pileup_up','Pileup_down','Pdfweight_up','Pdfweight_down','JES_up','JES_down','JMS_up','JMS_down','JER_up','JER_down','JMR_up','JMR_down','TptReweight_up','TptReweight_down']
    }

#'ElectronIDWeight_up','ElectronIDWeight_down','ElectronRecoWeight_up','ElectronRecoWeight_down','MuonIDWeight_up','MuonIDWeight_down','MuonRecoWeight_up','MuonRecoWeight_down'
    
    jets = {'Top':'p3','W':'p2','UM':'p1'} #naming conventions for SF tool

    for year in ['16','16APV','17','18']:
        outFile = ROOT.TFile('/uscms/home/mhesford/nobackup/XHYbbWW/ParticleNetSF/templates2D/particlenetmd_tt1l_{}_{}to1._20{}_200to1000_templates.root'.format(tagger,wp,year),'RECREATE')
        for region in ['pass','fail']:            
            for sys in systs[year]:
                if not sys.startswith('J'): #Not one of the JEC's - retreive this from the nominal file
                    sys_path = 'plots/PNetTemplates2D/PNetTemplates2D_{}_ttbar_{}.root'.format(tagger,year)
                else: #one of the JME variations, have to retrieve differently
                    sys_path = 'plots/PNetTemplates2D/PNetTemplates2D_{}_ttbar_{}_{}.root'.format(tagger,year,sys)
            
                sysFile = ROOT.TFile.Open(sys_path,'READ')
                for jet in jets:
                    if not sys.startswith('J'): 
                        sys_hist = '{}_{}__{}'.format(jet,region,sys)
                    else: #one of the JME variations, have to retrieve differently
                        sys_hist = '{}_{}__nominal'.format(jet,region) #retreive the nominal histogram from the file for the target systematic            

                    sysHist = sysFile.Get(sys_hist)
                    print(type(sysHist))
                    #sysHist.SetDirectory(0)
                    if sys == 'nominal':
                        sys_str = ''
                    else:
                        sys_str = sys.split('_')[0]+sys.split('_')[1].capitalize()
                    sysHist.SetName('tt_{}_{}_{}'.format(jets[jet],sys_str,region))
                    outFile.cd()
                    sysHist.Write()
        outFile.Close()

def MakeRun2(tagger, wp):
    #add all the data histograms together
    for year in ['16APV','16','17','18']:
        outFile = ROOT.TFile.Open('/uscms/home/mhesford/nobackup/XHYbbWW/ParticleNetSF/templates2D/particlenetmd_tt1l_{}_{}to1._20{}_200to1000_templates.root'.format(tagger,wp,year),'UPDATE')
        cmd = 'hadd -f plots/PNetTemplates2D/PNetTemplates2D_{0}_Data_{1}.root plots/PNetTemplates2D/PNetTemplates2D_{0}_Data*{1}.root'.format(tagger,year)
        ExecuteCmd(cmd)
        inFile = ROOT.TFile.Open('plots/PNetTemplates2D/PNetTemplates2D_{0}_Data_{1}.root'.format(tagger,year),'READ')
        for hist in ['pass','fail']:
            inHist = inFile.Get(hist+'__nominal')
            inHist.SetName('data_obs_{}'.format(hist))
            outFile.cd()
            inHist.Write()
        outFile.Close()
    
if __name__ == '__main__':
   
     
    redirector = 'root://cmseos.fnal.gov/'
    eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/plots/PNetTemplates2D/'
    
    cmd = 'eos {} ls {}PNetTemplates2D*'.format(redirector,eos_path)
    
    rawFiles = subprocess.check_output(cmd, shell=True).decode()
    files = rawFiles.split('\n')
    for fName in files:
        #Throwing all the output histograms from the selection script into a directory called "raw_selection"
        ExecuteCmd('xrdcp {}{}{} plots/PNetTemplates2D/'.format(redirector, eos_path, fName)) 
    
    tagger_wps = {'Hbb':'0.98','Wqq':'0.80'}
    for tagger, wp in tagger_wps.items():
        combine_ttbar(tagger)
        organize_ttbar_templates(tagger,wp)
        MakeRun2(tagger,wp)
    
