import ROOT
from TIMBER.Analyzer import Correction, CutGroup, VarGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME
from collections import OrderedDict

# Helper file for dealing with .txt files containing NanoAOD file locs
def SplitUp(filename,npieces,nFiles=False):
    '''Take in a txt file name where the contents are root
    file names separated by new lines. Split up the files
    into N lists where N is `npieces` in the case that `nFiles == False`.
    In the case that `nFiles == True`, `npieces` is treated as the
    number of files to have per list.
    '''
    files = open(filename,'r').readlines()
    nfiles = len(files)

    if npieces > nfiles:
        npieces = nfiles
    
    if not nFiles: files_per_piece = float(nfiles)/float(npieces)
    else: files_per_piece = npieces
    
    out = []
    iend = 0
    for ipiece in range(1,npieces+1):
        piece = []
        for ifile in range(iend,min(nfiles,int(ipiece*files_per_piece))):
            piece.append(files[ifile].strip())

        iend = int(ipiece*files_per_piece)
        out.append(piece)
    return out

class XHYbbWW:
    def __init__(self, inputfile, ijob=1, njobs=1):
        #inputfile format is '(raw_nano or snapshots)/SETNAME_YEAR(maybe _+something else).txt'
        infiles = SplitUp(inputfile,njobs)[ijob-1]
        self.setname = inputfile.split('/')[1].split('_')[0]
        self.year = inputfile.split('/')[1].split('_')[1].split('.')[0]
        self.ijob = ijob
        self.njobs = njobs

        self.a = analyzer(infiles)

        if 'Data' in inputfile:
          # self.a = analyzer(infiles)
            self.a.isData = True
        else:
          # MY=self.setname.split('-')[2]
          # self.a = analyzer(infiles,multiSampleStr=MY)
          # self.a.Cut('CorrectMass', 'GenModel_YMass_{} == 1'.format(MY)) 
            self.a.isData = False
 
	# Cuts from config file
        self.config = OpenJSON('XHYbbWWconfig.json')
        self.cuts = self.config['CUTS']

	# triggers
        self.trigs = self.config['TRIGS']

    CompileCpp('HWWmodules.cc') 

    def AddCutflowColumn(self, var, varName):
        '''
        for future reference:
        https://root-forum.cern.ch/t/rdataframe-define-column-of-same-constant-value/34851
        '''
        print('Adding cutflow information...\n\t{}\t{}'.format(varName, var))
        self.a.Define('{}'.format(varName), str(var))
 
    def getNweighted(self):
        if not self.a.isData:
            return self.a.DataFrame.Sum("genWeight").GetValue()
        else:
            return self.a.DataFrame.Count().GetValue()

    def BasicKinematics(self,LT='tight'): #nothing super selective, but should cut out some obviously uninteresting events
        #adding a loose/tight option to help fill out MX vs MY histograms for VV samples
        assert(LT == 'tight' or LT == 'loose')

        self.NSTART = self.getNweighted()
        self.AddCutflowColumn(self.NSTART,'NSTART')

        flags = [
            'Flag_goodVertices',
            'Flag_globalSuperTightHalo2016Filter',
            'Flag_HBHENoiseFilter',
            'Flag_HBHENoiseIsoFilter',
            'Flag_EcalDeadCellTriggerPrimitiveFilter',
            'Flag_BadPFMuonFilter',
            'Flag_BadPFMuonDzFilter',
            'Flag_eeBadScFilter'
        ]
        if self.year == '17' or self.year == '18':
            flags.append('Flag_ecalBadCalibFilter')
        MET_filters = self.a.GetFlagString(flags)       # string valid (existing in RDataFrame node) flags together w logical and
        self.a.Cut('flags', MET_filters)
        self.NFLAGS = self.getNweighted()
        self.AddCutflowColumn(self.NFLAGS, "NFLAGS")

        #at least one lepton passing pre-selection criteria (see HWWmodules.cc)
        if LT == 'tight':
            self.a.Define('isLeptonPre','isLeptonPreselected(nElectron, Electron_pt, Electron_eta, Electron_miniPFRelIso_all, nMuon, Muon_pt, Muon_eta, Muon_miniPFRelIso_all)')
        else:
            self.a.Define('isLeptonPre','isLooseLeptonPreselected(nElectron, Electron_pt, Electron_eta, Electron_miniPFRelIso_all, nMuon, Muon_pt, Muon_eta, Muon_miniPFRelIso_all)')

        self.a.Cut('isLeptonPreselected','isLeptonPre')
        self.LEPPRE=self.getNweighted()
        self.AddCutflowColumn(self.LEPPRE,"LEPPRE")

        #at least two fat jet w/ mass > 50
        if LT == 'tight':
            self.a.Define('FatJetHM_msoftdrop','FatJet_msoftdrop > 50')
            self.a.Cut('nFatJet_HighMass','FatJetHM_msoftdrop.size() > 1')
        else:
            self.a.Define('FatJetHM_msoftdrop','FatJet_msoftdrop > 10')
            self.a.Cut('nFatJet_HighMass','FatJetHM_msoftdrop.size() > 1')            
        self.JETMASS = self.getNweighted()
        self.AddCutflowColumn(self.JETMASS,'JETMASS')

        #high pt fat jets
        if LT == 'tight':
            self.a.Cut('Jet_pt','FatJet_pt[0] > 300 && FatJet_pt[1] > 200')
        else:
            self.a.Cut('Jet_pt','FatJet_pt[0] > 200 && FatJet_pt[1] > 100')
        self.JETPT=self.getNweighted()
        self.AddCutflowColumn(self.JETPT,"JETPT")

        #jet eta
        self.a.Cut('Jet_eta', 'abs(FatJet_eta[0]) < 2.4 && abs(FatJet_eta[1]) < 2.4') # no super forward jets
	self.JETETA = self.getNweighted()
	self.AddCutflowColumn(self.JETETA, "JETETA")

        #MET cut
        if LT == 'tight':
            self.a.Cut('MET_cut','MET_pt > 25') #This looks like it should kill some background for a minimal loss in signal
        else:
            self.a.Cut('MET_cut','MET_pt > 5')
        self.METPT = self.getNweighted()
        self.AddCutflowColumn(self.METPT, "METPT")
 
        return self.a.GetActiveNode()

    def ApplyStandardCorrections(self, snapshot=False):
        # first apply corrections for snapshot phase
        if snapshot:
            if self.a.isData:
                # NOTE: LumiFilter requires the year as an integer 
                lumiFilter = ModuleWorker('LumiFilter','TIMBER/Framework/include/LumiFilter.h',[int(self.year) if 'APV' not in self.year else 16])    # defaults to perform "eval" method 
                self.a.Cut('lumiFilter',lumiFilter.GetCall(evalArgs={"lumi":"luminosityBlock"}))	       # replace lumi with luminosityBlock
                if self.year == '18':
                    HEM_worker = ModuleWorker('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname if 'Muon' not in self.setname else self.setname[10:]])
                    self.a.Cut('HEM','%s[0] > 0'%(HEM_worker.GetCall(evalArgs={"FatJet_eta":"FatJet_eta","FatJet_phi":"FatJet_phi"})))
            else:
                self.a = ApplyPU(self.a, 'XHYbbWWpileup.root', '20{}'.format(self.year), ULflag=True, histname='{}_{}'.format(self.setname,self.year))
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    #self.a.AddCorrection(Correction("Prefire","TIMBER/Framework/include/Prefire_weight.h",[self.year],corrtype='weight'))
                    L1PreFiringWeight = Correction("L1PreFiringWeight","TIMBER/Framework/TopPhi_modules/BranchCorrection.cc",constructor=[],mainFunc='evalWeight',corrtype='weight',columnList=['L1PreFiringWeight_Nom','L1PreFiringWeight_Up','L1PreFiringWeight_Dn'])
                    self.a.AddCorrection(L1PreFiringWeight, evalArgs={'val':'L1PreFiringWeight_Nom','valUp':'L1PreFiringWeight_Up','valDown':'L1PreFiringWeight_Dn'})
                elif self.year == '18':
                    self.a.AddCorrection(Correction('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname],corrtype='corr'))

            self.a = AutoJME.AutoJME(self.a, 'FatJet', '20{}'.format(self.year), self.setname if 'Muon' not in self.setname else self.setname[10:])
            #self.a.MakeWeightCols(extraNominal='' if self.a.isData else 'genWeight*%s'%self.GetXsecScale())	
            #self.a.MakeWeightCols(extraNominal='genWeight' if not self.a.isData else '') #CHANGE THIS

        # now for selection
        else:
            if not self.a.isData:
                self.a.AddCorrection(Correction('Pileup',corrtype='weight'))
                self.a.AddCorrection(Correction('Pdfweight',corrtype='uncert'))
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    #self.a.AddCorrection(Correction('Prefire',corrtype='weight'))
		    # Instead, instantiate ModuleWorker to handle the C++ code via clang. This uses the branches already existing in NanoAODv9
                    self.a.AddCorrection(Correction('L1PreFiringWeight',corrtype='weight'))
                elif self.year == '18':
                    self.a.AddCorrection(Correction('HEM_drop',corrtype='corr'))
                #if 'ttbar' in self.setname:
                    #self.a.AddCorrection(Correction('TptReweight',corrtype='weight'))
        return self.a.GetActiveNode()

    def ApplyJMECorrections(self, variation):
	# for trigger effs
	#self.a.Define('Trijet_vect_trig','hardware::TLvector(Trijet_pt, Trijet_eta, Trijet_phi, Trijet_msoftdrop)')
	#self.a.Define('mhww_trig','hardware::InvariantMass(Trijet_vect_trig)')
	
        # JME variations - we only do this for MC
	if not self.a.isData:
	    # since H, W close enough in mass, we can treat them the same. 
	    # Higgs, W will have same pt and mass calibrations
	    # WARNING --------------------------------------------------------------------------------------------------------------
	    # IS THIS ACTUALLY TRUE?
            # IS IT, AMITAV? IS IT???
	    # ----------------------------------------------------------------------------------------------------------------------
	    pt_calibs, mass_calibs = JMEvariationStr('Higgs',variation)
	    self.a.Define('FatJet_pt_corr','hardware::MultiHadamardProduct(FatJet_pt,{})'.format(pt_calibs))
	    self.a.Define('FatJet_msoftdrop_corr','hardware::MultiHadamardProduct(FatJet_msoftdrop,{})'.format(mass_calibs))
	else:
	    self.a.Define('FatJet_pt_corr','hardware::MultiHadamardProduct(FatJet_pt,{FatJet_JES_nom})')
	    self.a.Define('FatJet_msoftdrop_corr','hardware::MultiHadamardProduct(FatJet_msoftdrop,{FatJet_JES_nom})')
	return self.a.GetActiveNode()

    # for trigger efficiencies
    def ApplyTrigs(self):
        #For now, just apply triggers straight up
        #Eventually, once trigger efficiencies in data measured, can just apply them as a weight correction to MC
        all_trigs = self.trigs[int(self.year) if 'APV' not in self.year else 16]
        trigs = []
        for group in ['HADRONIC','ELECTRON','MUON']:
            trigs.append(trig.encode('ascii','ignore') for trig in all_trigs[group]) #something about unicode	
        self.a.Cut('trigger',self.a.GetTriggerString(trigs))
	return self.a.GetActiveNode()

    def Snapshot(self, node=None):
        #for now, snapshot without applying any corrections (deal w/ this stuff later)
        startNode = self.a.GetActiveNode()
        if node == None:
            node = self.a.GetActiveNode()
              
        columns = [
        'nJet','Jet_btagDeepFlavB','nFatJet','FatJet_eta','FatJet_msoftdrop','FatJet_pt','FatJet_phi',
        'FatJet_JES_nom','FatJet_deepTagMD_HbbvsQCD', 'FatJet_deepTagMD_ZHbbvsQCD',
        'FatJet_deepTagMD_WvsQCD', 'FatJet_deepTag_TvsQCD', 'FatJet_particleNet_HbbvsQCD',
        'FatJet_particleNet_TvsQCD', 'FatJet_particleNet_QCD',
        'FatJet_particleNet_WvsQCD','FatJet_particleNetMD*', 'FatJet_rawFactor', 'FatJet_tau*',
        'FatJet_jetId','nMuon','Muon_*','nElectron','Electron_*','HLT_Mu50','HLT_TkMu*','HLT_IsoMu*','HLT_OldMu100',
        'MET_*','HLT_*MET*','HLT_Photon*','HLT_Ele*','HLT_PFHT*', 'HLT_PFJet*', 'HLT_AK8*','genWeight',
        'event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','LEPPRE','JETMASS',
        'JETPT','JETETA','METPT'
        ]         
       
	#'FatJet_JES_nom', 
        # append to columns list if not Data
        if not self.a.isData:
            columns.extend(['GenPart_*', 'nGenPart', 'genWeight']) 
            columns.extend(['FatJet_JES_up','FatJet_JES_down',
                            'FatJet_JER_nom','FatJet_JER_up','FatJet_JER_down',
                            'FatJet_JMS_nom','FatJet_JMS_up','FatJet_JMS_down',
                            'FatJet_JMR_nom','FatJet_JMR_up','FatJet_JMR_down'])
            columns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__nom','Pdfweight__up','Pdfweight__down'])
            if self.year == '16' or self.year == '17' or 'APV' in self.year:
		#columns.extend(['Prefire__nom','Prefire__up','Prefire__down'])
                columns.extend(['L1PreFiringWeight_Nom', 'L1PreFiringWeight_Up', 'L1PreFiringWeight_Dn'])	# these are the default columns in NanoAODv9
                columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])    # these are the weight columns created by the BranchCorrection module
            elif self.year == '18':
                columns.append('HEM_drop__nom')
            
        # get ready to send out snapshot
        self.a.SetActiveNode(node)
        self.a.Snapshot(columns, 'XHYbbWWsnapshot_{}_{}_{}_{}.root'.format(self.setname,self.year,self.ijob,self.njobs),'Events', openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

    def Dijets(self):

        #Pick highest pt back to back fat jets
        self.a.Define('DijetIdxs','PickDijets(FatJet_eta,FatJet_phi,FatJet_msoftdrop)')
        self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')

        self.nDijets = self.getNweighted()
        self.AddCutflowColumn(self.nDijets, "nDijets")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

        #Define mass-decorrelated Higgs/Wqq taggers
        #Mass-decorrelated W tagger discriminant is defined by inclusion of X->cc
        #See slide 16: https://indico.cern.ch/event/809820/contributions/3632617/attachments/1970786/3278138/MassDecorrelation_ML4Jets_H_Qu.pdf
        #Right now I'm testing some older MC samples for WZ/ZZ-->llqq which do not have particleNet tagger info - run instead using deepTagMD, will treat exactly the same as particleNetMD for these initial tests
        if '2L2Q' in self.setname:
            #Falsely renaming deepTagMD so I don't have to change any of the other code 
            self.a.Define('Dijet_particleNetMD_HbbvsQCD','Dijet_deepTagMD_HbbvsQCD > 0')
            self.a.Define('Dijet_particleNetMD_WqqvsQCD','Dijet_deepTagMD_WvsQCD > 0')
        else:
            self.a.Define('Dijet_particleNetMD_HbbvsQCD','Dijet_particleNetMD_Xbb/(Dijet_particleNetMD_Xbb+Dijet_particleNetMD_QCD)')
            self.a.Define('Dijet_particleNetMD_WqqvsQCD','(Dijet_particleNetMD_Xqq+Dijet_particleNetMD_Xcc)/(Dijet_particleNetMD_Xqq+Dijet_particleNetMD_Xcc+Dijet_particleNetMD_QCD)')
   
        #mark one jet as Higgs and other as Wqq
        self.a.Define('WqqIdx','Dijet_particleNetMD_WqqvsQCD[0] > Dijet_particleNetMD_WqqvsQCD[1] ? 0 : 1')
        self.a.Define('HiggsIdx','WqqIdx == 0 ? 1 : 0')

        self.a.ObjectFromCollection('Higgs','Dijet','HiggsIdx')
        self.a.ObjectFromCollection('Wqq','Dijet','WqqIdx')
 
        return self.a.GetActiveNode()

    def KinematicLepton(self):

        self.a.Define('kinEleIdx','kinElectron(Electron_pt,Electron_eta,Electron_phi,Higgs_phi,Wqq_phi)')
        self.a.Define('kinMuIdx','kinMuon(Muon_pt,Muon_eta,Muon_phi,Higgs_phi,Wqq_phi)')  
        self.a.Cut('kinLepton_cut','kinEleIdx != -1 || kinMuIdx != -1') #at least one good lepton
        self.a.Define('LeptonType','LeptonIdx(kinEleIdx,kinMuIdx,Electron_pt,Muon_pt)') #picks higher pt signal lepton - output = 0 (lepton is electron) or 1 (lepton is muon)

        #self.a.Cut('lepIsolation','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < {} : Muon_miniPFRelIso_all[kinMuIdx] < {}'.format(self.cuts['ELECTRON']['RelIso'],self.cuts['MUON']['RelIso'])) #Right now, save this for signal/control region definitions
        #self.a.Cut('leptonQuality','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]') #"naively" using the standard WP80 boolean for now; eventually must examine more closely efficiency of lepton cuts and find optimal cut value(s) for Electron_mvaFall17V2noIso discriminant

        self.nkinLep = self.getNweighted()
        self.AddCutflowColumn(self.nkinLep,'nkinLep')

        #For ease, merge some lepton columns that will be useful later (for lepton-type specific variables, use LeptonType to determine if electron or muon)
        self.a.Define('Lepton_pt','LeptonType == 1 ? Muon_pt[kinMuIdx] : Electron_pt[kinEleIdx]')
        self.a.Define('Lepton_eta','LeptonType == 1 ? Muon_eta[kinMuIdx] : Electron_eta[kinEleIdx]')
        self.a.Define('Lepton_phi','LeptonType == 1 ? Muon_phi[kinMuIdx] : Electron_phi[kinEleIdx]')
        self.a.Define('Lepton_mass','LeptonType == 1 ? Muon_mass[kinMuIdx] : Electron_mass[kinEleIdx]')

        return self.a.GetActiveNode()

    def JetSelection(self,mW,mH,Wtag,Htag): #ex: mW = [60,100], Wtag = 0.8

        #mass cuts
        self.a.Cut('WqqMass','Wqq_msoftdrop > {} && Wqq_msoftdrop < {}'.format(mW[0],mW[1]))
        self.a.Cut('HiggsMass','Higgs_msoftdrop > {} && Higgs_msoftdrop < {}'.format(mH[0],mH[1]))

        self.JETMASS = self.getNweighted()
        self.AddCutflowColumn(self.JETMASS,'JETMASS')

        #tagging cuts
        self.a.Cut('WqqTag','Wqq_particleNetMD_WqqvsQCD > {}'.format(Wtag))
        self.a.Cut('HiggsTag','Higgs_particleNetMD_HbbvsQCD > {}'.format(Htag))

        self.JETTAG = self.getNweighted()
        self.AddCutflowColumn(self.JETTAG,'JETTAG')

        return self.a.GetActiveNode()

    def GetXsecScale(self):
        lumi = self.config['lumi{}'.format(self.year if 'APV' not in self.year else 16)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        return lumi*xsec/self.a.genEventSumw

    def ApplyMassCuts(self):
        # perform Higgs mass window cut, save cutflow info
	self.a.Cut('mH_{}_cut'.format('window'),'Higgs_msoftdrop > {0} && Higgs_msoftdrop < {1}'.format(*self.cuts['JET']['mh']))
	self.nHiggs = self.getNweighted()
	self.AddCutflowColumn(self.nHiggs, 'nHiggsMassCut')
	# Lead W mass window cut, cutflow
	self.a.Cut('mW_{}_cut'.format('window'),'Wqq_msoftdrop > {0} && Wqq_msoftdrop < {1}'.format(*self.cuts['JET']['mw']))
	self.nWqq = self.getNweighted()
	self.AddCutflowColumn(self.nWqq, 'nWqqMassCut')
	
        return self.a.GetActiveNode()

    def ApplySRorCR(self, SRorCR, tagger): #This is getting frequently updated as I play w/ different definitions
	'''
	SRorCR [str] = 'SR' or 'CR', used to generate cutflow information
	tagger [str] = name of tagger to be used 
	Selection criteria:
		SR: Wqq > 0.8
                    Electron_mvaFall17V2noIso_WP80 || Muon_MediumId  
                    Lepton_miniPFRelIso_all < 0.2
		CR: 0.4 < Wqq < 0.8
                    Electron_mvaFall17V2noIso_WP80 || Muon_MediumId #want good quality leptons even in conrol region
                    0.2 < Lepton_miniPFRelIso_all < 0.5
                ttCR: Wqq < 0.8
                    Electron_mvaFall17V2noIso_WP80 || Muon_MediumId #keep normal lepton cuts
                    Lepton_miniPFRelIso_all < 0.2
                    b-tagged AK4 jet exists
                Note: We want the control region to be void of signal/high in background, but we don't want it to be so different from the signal region. 
                      It should be similar enough to the signal region, but statistically disjunct.
	'''
	assert(SRorCR=='SR' or SRorCR=='CR' or SRorCR=='ttCR')
	# tagger values for each 
	if SRorCR == 'SR':
	    # Signal region
	    self.a.Cut('Wqq_{}_cut_{}'.format(tagger,SRorCR), 'Wqq_{}_WqqvsQCD > {}'.format(tagger,self.cuts['JET']['particleNetMD_WqqvsQCD'][1]))
            self.nWtag_SR = self.getNweighted()
            self.AddCutflowColumn(self.nWtag_SR, 'nWtag_SR')

            self.a.Cut('lepIsolation_{}'.format(SRorCR),'LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < {} : Muon_miniPFRelIso_all[kinMuIdx] < {}'.format(self.cuts['ELECTRON']['RelIso'][0],self.cuts['MUON']['RelIso'][0]))
            self.nlepIso_SR = self.getNweighted()
            self.AddCutflowColumn(self.nlepIso_SR, 'nlepIso_SR')

            self.a.Cut('leptonQuality_{}'.format(SRorCR),'LeptonType == 0? {}[kinEleIdx] : {}[kinMuIdx]'.format(self.cuts['ELECTRON']['mva80'],self.cuts['MUON']['mediumId']))
            self.nlepQ_SR = self.getNweighted()
            self.AddCutflowColumn(self.nlepQ_SR, 'nlepQ_SR')

            #This seems to be quite bad for signal
            #self.a.Define('Jet_btag','Jet_btagDeepFlavB > {}'.format(self.cuts['JET']['btag'])) #we want to veto events with a b-tagged AK4 jet, which we will use to define the ttbar control region
            #self.a.Cut('nbtaggedJet','Jet_btag.size() == 0')
           
	elif SRorCR == 'CR':
	    # Control Region
	    self.a.Cut('Wqq_{}_cut_{}'.format(tagger,SRorCR), 'Wqq_{0}_WqqvsQCD > {1} && Wqq_{0}_WqqvsQCD < {2}'.format(tagger,self.cuts['JET']['particleNetMD_WqqvsQCD'][0],self.cuts['JET']['particleNetMD_WqqvsQCD'][1]))
            self.nWtag_CR = self.getNweighted()
            self.AddCutflowColumn(self.nWtag_CR, 'nWtag_CR')

            self.a.Cut('lepIsolation_{}'.format(SRorCR),'LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] > {} && Electron_miniPFRelIso_all[kinEleIdx] < {} : Muon_miniPFRelIso_all[kinMuIdx] > {} && Muon_miniPFRelIso_all[kinMuIdx] < {}'.format(self.cuts['ELECTRON']['RelIso'][0], self.cuts['ELECTRON']['RelIso'][1], self.cuts['MUON']['RelIso'][0], self.cuts['MUON']['RelIso'][1]))
            self.nlepIso_CR = self.getNweighted()
            self.AddCutflowColumn(self.nlepIso_CR, 'nlepIso_CR')

            self.a.Cut('leptonQuality_{}'.format(SRorCR),'LeptonType == 0? {}[kinEleIdx] : {}[kinMuIdx]'.format(self.cuts['ELECTRON']['mva80'],self.cuts['MUON']['mediumId']))
            self.nlepQ_CR = self.getNweighted()
            self.AddCutflowColumn(self.nlepQ_CR, 'nlepQ_CR')       
 
        else:
            #ttbar control region
            self.a.Cut('Wqq_{}_cut_ttCR'.format(tagger), 'Wqq_{}_WqqvsQCD < {}'.format(tagger,self.cuts['JET']['particleNetMD_WqqvsQCD'][1]))
            self.nWtag_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nWtag_ttCR, 'nWtag_ttCR')

            self.a.Cut('lepIsolation_ttCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < {} : Muon_miniPFRelIso_all[kinMuIdx] < {}'.format(self.cuts['ELECTRON']['RelIso'][0],self.cuts['MUON']['RelIso'][0]))
            self.nlepIso_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nlepIso_ttCR, 'nlepIso_ttCR')

            self.a.Cut('leptonQuality_ttCR'.format(SRorCR),'LeptonType == 0? {}[kinEleIdx] : {}[kinMuIdx]'.format(self.cuts['ELECTRON']['mva80'],self.cuts['MUON']['mediumId']))
            self.nlepQ_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nlepQ_ttCR, 'nlepQ_ttCR')

            self.a.Define('Jet_btag','Jet_btagDeepFlavB > {}'.format(self.cuts['JET']['btag']))
            self.a.Cut('nbtaggedJet','Jet_btag.size() > 0') #must define signal region in advance
            self.nJetB_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nJetB_ttCR, 'nJetB_ttCR')

        return self.a.GetActiveNode()

    def ApplyPassFail(self, SRorCR, tagger): #also getting changed pretty frequently for now
	'''
	Fail:	Hbb < 0.94
	Pass: 	Hbb > 0.94
	'''
	assert(SRorCR=='SR' or SRorCR=='CR' or SRorCR=='ttCR')
	checkpoint = self.a.GetActiveNode()
	FP = OrderedDict()
	# Higgs fail + cutflow info
	FP['fail_'+SRorCR] = self.a.Cut('HbbTag_fail','Higgs_{0}_HbbvsQCD < {1}'.format(tagger, self.cuts['JET']['particleNetMD_HbbvsQCD']))
	if SRorCR=='SR':
	    self.nHF_SR = self.getNweighted()
	    self.AddCutflowColumn(self.nHF_SR, 'higgsF_SR')
	else:
	    self.nHF_CR = self.getNweighted()
	    self.AddCutflowColumn(self.nHF_CR, 'higgsF_CR')

	# Higgs Pass + cutflow
	self.a.SetActiveNode(checkpoint)
	FP['pass_'+SRorCR] = self.a.Cut('HbbTag_pass','Higgs_{0}_HbbvsQCD > {1}'.format(tagger, self.cuts['JET']['particleNetMD_HbbvsQCD']))
	if SRorCR == 'SR':
	    self.nHP_SR = self.getNweighted()
	    self.AddCutflowColumn(self.nHP_SR, 'higgsP_SR')
	else:
	    self.nHP_CR = self.getNweighted()
	    self.AddCutflowColumn(self.nHP_CR, 'higgsP_CR')

	# reset state, return dict
	self.a.SetActiveNode(checkpoint)
	return FP

    def get_standard_int(self,region,binsX,binsY,node):
        #The idea behind this function is to perform all the "standard cuts" for a given region (pass/fail for SR, CR, or ttCR) and produce a histogram
        #Then, we count the total # events in the histogram so that when we apply comparatively looser cuts to the MX vs MY phase space we can scale by (# events passing loose selection)/(# events passing standard selection)
        self.a.SetActiveNode(node)
        #first perform all universal cuts
        self.a.Define('isLeptonPre_{}'.format(region),'isLeptonPreselected(nElectron, Electron_pt, Electron_eta, Electron_miniPFRelIso_all, nMuon, Muon_pt, Muon_eta, Muon_miniPFRelIso_all)')
        self.a.Cut('isLeptonPreselected_{}'.format(region),'isLeptonPre_{}'.format(region))
        self.a.Define('FatJetHM_msoftdrop_{}'.format(region),'FatJet_msoftdrop > 50')
        self.a.Cut('nFatJet_HighMass_{}'.format(region),'FatJetHM_msoftdrop_{}.size() > 1'.format(region))
        self.a.Cut('Jet_pt_{}'.format(region),'FatJet_pt[0] > 300 && FatJet_pt[1] > 200')
        self.a.Cut('MET_cut_{}'.format(region),'MET_pt > 25')
        self.a.Cut('Lepton_pt_{}'.format(region),'Lepton_pt > 25')
        self.a.Cut('mH_{}'.format(region),'Higgs_msoftdrop > 100 && Higgs_msoftdrop < 145')
        self.a.Cut('mW_{}'.format(region),'Wqq_msoftdrop > 65 && Wqq_msoftdrop < 100')

        if region == 'pass_SR':
            self.a.Cut('Wtag_passSR','Wqq_particleNetMD_WqqvsQCD > 0.8')
            self.a.Cut('lepIso_passSR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_passSR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_passSR','Higgs_particleNetMD_HbbvsQCD > 0.94')
        elif region == 'pass_CR':
            self.a.Cut('Wtag_passCR','Wqq_particleNetMD_WqqvsQCD > 0.4 && Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_passCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] > 0.2 && Electron_miniPFRelIso_all[kinEleIdx] < 0.5 : Muon_miniPFRelIso_all[kinMuIdx] > 0.2 && Muon_miniPFRelIso_all[kinMuIdx] < 0.5')
            self.a.Cut('lepQ_passCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_passCR','Higgs_particleNetMD_HbbvsQCD > 0.94')
        elif region == 'fail_SR':
            self.a.Cut('Wtag_failSR','Wqq_particleNetMD_WqqvsQCD > 0.8')
            self.a.Cut('lepIso_failSR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_failSR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_failSR','Higgs_particleNetMD_HbbvsQCD < 0.94')
        elif region == 'fail_CR':
            self.a.Cut('Wtag_failCR','Wqq_particleNetMD_WqqvsQCD > 0.4 && Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_failCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] > 0.2 && Electron_miniPFRelIso_all[kinEleIdx] < 0.5 : Muon_miniPFRelIso_all[kinMuIdx] > 0.2 && Muon_miniPFRelIso_all[kinMuIdx] < 0.5')
            self.a.Cut('lepQ_failCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_failCR','Higgs_particleNetMD_HbbvsQCD < 0.94')
        elif region == 'pass_ttCR':
            self.a.Cut('Wtag_passttCR','Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_passttCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_passttCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Define('Jet_btag_passttCR','Jet_btagDeepFlavB > 0.8')
            self.a.Cut('nbtaggedJet_passttCR','Jet_btag_passttCR.size() > 0')
            self.a.Cut('Htag_passttCR','Higgs_particleNetMD_HbbvsQCD > 0.94')
        else: #region == 'fail_ttCR'
            self.a.Cut('Wtag_failttCR','Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_failttCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_failttCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Define('Jet_btag_failttCR','Jet_btagDeepFlavB > 0.8')
            self.a.Cut('nbtaggedJet_failttCR','Jet_btag_failttCR.size() > 0')
            self.a.Cut('Htag_failttCR','Higgs_particleNetMD_HbbvsQCD < 0.94')

        syst_hists = self.a.MakeTemplateHistos(ROOT.TH2F('MXvMY_{}'.format(region),'temp',binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),['mhwlv','mwlv'])
        stdInts = {}
        for name in syst_hists.keys():
            hist=syst_hists[name]
            hist.SetDirectory(0)
            stdInts[name] = hist.Integral()
        
        self.a.SetActiveNode(node)
        return stdInts

# for use in selection - essentially just creates combinations of all the JME variations
def JMEvariationStr(p, variation):
    base_calibs = ['FatJet_JES_nom','FatJet_JER_nom','FatJet_JMS_nom','FatJet_JMR_nom']
    variationType = variation.split('_')[0]
    pt_calib_vect = '{'
    mass_calib_vect = '{'
    for c in base_calibs:
	if 'JM' in c and p != 'Top':	# WARNING - might need to change this if we treat W, H differently for mass and pt calibrations  
	    mass_calib_vect += '%s,'%('FatJet_'+variation if variationType in c else c)
	elif 'JE' in c:
	    pt_calib_vect += '%s,'%('FatJet_'+variation if variationType in c else c)
	    mass_calib_vect += '%s,'%('FatJet_'+variation if variationType in c else c)
    pt_calib_vect = pt_calib_vect[:-1]+'}'
    mass_calib_vect = mass_calib_vect[:-1]+'}'
    return pt_calib_vect, mass_calib_vect
