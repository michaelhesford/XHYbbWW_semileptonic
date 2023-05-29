'''
Script w/ functions to produce N-1 plots for both the dijet pair and the highest pT isolated lepton. 
'''

import ROOT
from TIMBER.Analyzer import HistGroup, CutGroup
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from XHYbbWW_class import XHYbbWW

cuts = {'mHiggs':[100,150],
        'mWqq':[60,100],
        'Wqq_particleNetMD_WqqvsQCD':0.8,
        'Higgs_particleNetMD_HbbvsQCD':0.8,
        'Lepton_miniPFRelIso_all':0.2,
        'DeltaPhi_Lepton_Higgs':1.571, #pi/2
        'DeltaPhi_Lepton_Wqq':[1,2.5]
}

def XHYbbWW_studies(self):

    ##### NEW VARIABLES FOR LATER USE #####

    #W_leptonic transverse mass
    self.a.Define('W_massTran','TransverseMass(MET_pt,IsoLepton_pt,MET_phi,IsoLepton_phi)') #Transverse W mass
    self.a.Define('W_massTran_genMET','TransverseMass(MET_fiducialGenPt,IsoLepton_pt,MET_fiducialGenPhi,IsoLepton_phi)') #using generator-level MET variables

    #Lorentz 4-vectors
    self.a.Define('MET_vect','hardware::TLvector(MET_pt,0,MET_phi,0)') #neutrino mass negligable, for now assuming MET_eta = 0 (p_z = 0)
    self.a.Define('Lepton_vect','hardware::TLvector(IsoLepton_pt,IsoLepton_eta,IsoLepton_phi,IsoLepton_mass)')
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt,Wqq_eta,Wqq_phi,Wqq_msoftdrop)')
    self.a.Define('Hbb_vect','hardware::TLvector(Higgs_pt,Higgs_eta,Higgs_phi,Higgs_msoftdrop)')

    #Invariant masses of W/Y/X
    self.a.Define('W_massInv','hardware::InvariantMass({MET_vect,Lepton_vect})') #full invariant mass
    self.a.Define('Y_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect})')
    self.a.Define('X_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect,Hbb_vect})')

    studiesPlots = HistGroup('studiesPlots')

    #######################################

    taggers = ['particleNetMD']
    # N-1 - essentially set up mechanics with JetNminus1Group() and LeptonNminus1Group () to automate the process via TIMBER -> a.Nminus1()
    for t in taggers:
        Higgs_tagger = '{}_HbbvsQCD'.format(t)
        Wqq_tagger = '{}_WqqvsQCD'.format(t)

        nminusGroup = self.Nminus1Group(t) # dictionary of nodes
        nminusNodes = self.a.Nminus1(nminusGroup) # TIMBER can now be fed the dict and automatically do N-1
        for n in nminusNodes.keys():
            if n.startswith('m'):    # mass cut
                bins = [25,50,300]
                if n.startswith('mWqq'):
                    var = 'Wqq_msoftdrop'
                elif n.startswith('mHiggs'):
                    var = 'Higgs_msoftdrop'
            elif n == 'full':   # full cuts, has all N of them

                # this key ALWAYS exists in nminusNodes by default. This is the node which has all of the cuts applied
		'''
			|
		       / \
		      /   \
		     /     \
		 [cut1]  [cut2]
		  /   \      \
			      \
			       \
			     [cut3] 
				 \
				['full'] <- this is the node with every cut made
		'''
                # this will effectively plot the Y mass (W1+W2) with EVERY other cut that we've defined
                #	- mH, mWqq, Higgs_tag, Wqq_tag
                var = ['Y_mass','W_massInv','W_massTran','W_massTran_genMET']
            elif n.startswith('DeltaPhi'):
                bins = [30,-3.2,3.2]
                if n.endswith('Wqq_cut'):
                    var = 'DeltaPhi_IsoLepton_Wqq'
                elif n.endswith('Higgs_cut'):
                    var = 'DeltaPhi_IsoLepton_Higgs'
            elif n.startswith('Isolation'):
                bins = [20,0,1]
                var = 'IsoLepton_miniPFRelIso_all'
            else: # tagger cut
                bins = [50,0,1]
                if n.startswith('Higgs'):
                    var = 'Higgs_{}'.format(Higgs_tagger)
                elif n.startswith('Wqq'):
                    var = 'Wqq_{}'.format(Wqq_tagger)
            if (type(var) == list):
                 for v in var:
                    print('N-1: Plotting {} for node {}'.format(v,n))
                    if v == 'Y_mass':
                        if 'XHY' in self.setname:
                            if int(self.setname.split('_')[0].split('-')[2]) > 1000: #expected y mass > 1000   
                                bins = [60,100,3000]
                            else:
                                bins = [30,100,1600]
                        else:
                            bins = [30,100,1600]
                    else: #looking at W masses
                        if 'XHY' in self.setname: #looking at signal
                            if int(self.setname.split('_')[0].split('-')[2]) > 1000: #expected y mass > 1000   
                                bins = [25,0,1000]
                            else:
                                bins = [20,0,800]
                        else: #looking at background
                            bins = [20,0,800]
                    studiesPlots.Add('{}_allCuts'.format(v),nminusNodes[n].DataFrame.Histo1D(('{}_allCuts'.format(v),'{}_allCuts'.format(v),bins[0],bins[1],bins[2]),v,'weight__nominal'))
            else:
                print('N-1: Plotting {} for node {}'.format(var, n))
                studiesPlots.Add(n+'_nminus1',nminusNodes[n].DataFrame.Histo1D((n+'_nminus1',n+'_nminus1',bins[0],bins[1],bins[2]),var,'weight__nominal'))

	# now we want to plot mX vs mY for QCD, ttbar, and signal
        for t in taggers:
            self.ApplyMassCuts()
            start=self.a.GetActiveNode()

            # We use Wqq tagging scores to divide data into two regions: signal (enriched in signal) and control (enriched in background)
            #    - Signal:    Wqq > 0.8
            #    - Control:    Wqq < 0.8
            # We define a pass/fail criteria for the Hbb score within each region 
            #   - Region 1 (fail):      Hbb < 0.94
            #   - Region 2 (pass):      Hbb > 0.94

            SR=self.ApplyWTag('SR',t)
            SR_FP=self.ApplyHiggsTag('SR',t)
            
            self.a.SetActiveNode(start)
            CR=self.ApplyWTag('CR',t)
            CR_FP=self.ApplyHiggsTag('CR',t)
  
            nodes={}
            nodes.update(SR_FP)
            nodes.update(CR_FP)

            bins = [80,0,4500]           

            for node in nodes.keys():
                self.a.SetActiveNode(nodes[node])
                print('MX vs MY: Plotting for {}'.format(node))
                studiesPlots.Add('MXvsMY_{}'.format(node), self.a.DataFrame.Histo2D(('MXvsMY_{}'.format(node), 'X vs Y Invariant Mass - {} {}'.format(node.split('_')[1],node.split('_')[0]), bins[0], bins[1], bins[2], bins[0], bins[1], bins[2]), 'X_mass', 'Y_mass', 'weight__nominal'))           
        outFile = ROOT.TFile.Open('{}_{}_{}_studies.root'.format(self.setname,self.year,self.ijob),'RECREATE')
        outFile.cd()
        studiesPlots.Do('Write') 
        self.a.PrintNodeTree('NodeTree.pdf',verbose=True)
        outFile.Close()    

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
        action='store',help='name of data set to run on')
    parser.add_argument('-y', type=str, dest='year',
        action='store', help='year',required=False)
    parser.add_argument('-j', type=int, dest='ijob',required=False,
        action='store', help='current job')
    parser.add_argument('-n', type=int, dest='njobs',required=False,
        action='store', help='number of jobs')
    args = parser.parse_args()

    setname=args.setname
    year=args.year
    ijob=args.ijob
    njobs=args.njobs

    ana = XHYbbWW(setname,year,ijob,njobs)

    ana.BasicKinematics()
    ana.ApplyStandardCorrections()
    ana.Dijets()
    ana.SignalLepton()

    XHYbbWW_studies(ana)   
