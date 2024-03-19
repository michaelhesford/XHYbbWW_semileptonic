#Final selection script, outputs templates for background estimation

import ROOT, time
from TIMBER.Analyzer import HistGroup, Correction
from TIMBER.Tools.Common import CompileCpp, ExecuteCmd
from collections import OrderedDict
import TIMBER.Tools.AutoJME as AutoJME
from XHYbbWW_class import XHYbbWW

def Make2DHistos(self,variation,tagger,wp):
    #Tagger = "Hbb" or "Wqq"
    #wp used to define pass/fail regions (eg. 0.98)
    #outFile = ROOT.TFile.Open('PNetTemplates2D_{}_{}_{}_{}_{}{}.root'.format(tagger, self.setname, self.year, ijob, njobs, '_'+variation if variation != 'None' else ''),'RECREATE')
    outFile = ROOT.TFile.Open('PNetTemplates2D_{}_{}_{}{}.root'.format(tagger, self.setname, self.year, '_'+variation if variation != 'None' else ''),'RECREATE')
    outFile.cd()

    #Divide probe jet into categories used to calculate SF's separately
        #Top merged: 0    
        #W merged: 1
        #bq merged: 2
        #Unmerged: 3

    if not self.a.isData:
        #Find true MC flavor for fat jet
        self.a.Define('ProbeJet_jetType', 'JetFlavor_ttbar(ProbeJet_eta, ProbeJet_phi, GenPart_eta, GenPart_phi, GenPart_pdgId, GenPart_genPartIdxMother)')
 
        '''
        #Create separate data columns for each type
        self.a.ObjectFromCollection('TopJet','ProbeJet','ProbeJet_jetType == 0')
        self.a.ObjectFromCollection('WJet','ProbeJet','ProbeJet_jetType == 1')
        self.a.ObjectFromCollection('UMJet','ProbeJet','ProbeJet_jetType == 3')
        '''

    start = self.a.GetActiveNode() 

    nodes = OrderedDict()
    #create pass region
    self.a.SetActiveNode(start)
    PR=self.a.Cut('wp_pass','ProbeJet_particleNetMD_{}vsQCD > {}'.format(tagger,wp))
    nodes['pass'] = PR
       
    #create fail region
    self.a.SetActiveNode(start)
    FR=self.a.Cut('wp_fail','ProbeJet_particleNetMD_{}vsQCD < {}'.format(tagger,wp))        
    nodes['fail'] = FR
    
    binsX = [40,50,250] #mass
    binsY = [8,200,1000] #pt
    
    for region in nodes.keys():
        print('Making template histos for {} region'.format(region))
        self.a.SetActiveNode(nodes[region])
        if not self.a.isData:
            for jet, val in {'Top':0, 'W':1, 'BQ':2, 'UM':3}.items():
                self.a.SetActiveNode(nodes[region])
                self.a.Cut('{} jet type'.format(jet),'ProbeJet_jetType == {}'.format(val))
                print('{} jet'.format(jet))
                plot_vars = ['ProbeJet_msoftdrop_corr','ProbeJet_pt_corr']
                templates = self.a.MakeTemplateHistos(ROOT.TH2F('{}_{}'.format(jet,region), '{}_{}'.format(jet,region), binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),plot_vars)
                templates.Do('Write')
        else:
            plot_vars = ['ProbeJet_msoftdrop_corr','ProbeJet_pt_corr']
            templates = self.a.MakeTemplateHistos(ROOT.TH2F(region, region, binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),plot_vars)
            templates.Do('Write')
    

    outFile.Close()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
                        action='store', required=True,
                        help='Setname to process.')
    parser.add_argument('-y', type=str, dest='year',
                        action='store', required=True,
                        help='Year of set (16, 17, 18).')
    parser.add_argument('-v', type=str, dest='variation',
                        action='store', default='None',
                        help='JES_up, JES_down, JMR_up,...')
    parser.add_argument('-t', type=str, dest='tagger',
                        action='store', required=True,
                        help='Hbb or Wqq')
    parser.add_argument('-j', type=int, dest='ijob',required=False,
                        action='store', default=1, help='current job')
    parser.add_argument('-n', type=int, dest='njobs',required=False,
                        action='store', default=1, help='number of jobs')
    args = parser.parse_args()
    setname=args.setname
    year=args.year
    variation = args.variation
    tagger = args.tagger
    ijob = args.ijob
    njobs = args.njobs
    
    filename='PNetSnapshots/{}_{}_PNetSnapshot.txt'.format(setname,year)

    ana = XHYbbWW(filename,ijob,njobs)
    ana.ApplyTrigs(ApplyLeptonTrigs = False)
    ana.ApplyJMECorrections(variation, jets = ['ProbeJet'])
    ana.ApplyStandardCorrections(snapshot=False)
    ana.ApplyTopPtReweight('ProbeJet','BJet',scale = 1, isJet1AK4 = True)
    #ana.ApplyLeptonCorrections()
    ana.a.MakeWeightCols(extraNominal='' if ana.a.isData else 'genWeight*%s'%ana.GetXsecScale())

    wp = ana.cuts['JET']['particleNetMD_{}vsQCD'.format(tagger)]

    Make2DHistos(ana,variation,tagger,wp)
   
 
