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
    def __init__(self, inputfile, ijob, njobs):
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
        self.trigs = {
            16:['HLT_PFHT800','HLT_PFHT900'],
            17:["HLT_PFHT1050","HLT_AK8PFJet500","HLT_AK8PFHT750_TrimMass50","HLT_AK8PFHT800_TrimMass50","HLT_AK8PFJet400_TrimMass30"],
            18:['HLT_AK8PFJet400_TrimMass30','HLT_AK8PFHT850_TrimMass50','HLT_PFHT1050']
        }

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

    def BasicKinematics(self): #nothing super selective, but should cut out some obviously uninteresting events

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
        self.a.Define('isLeptonPre','isLeptonPreselected(nElectron, Electron_pt, Electron_eta, Electron_miniPFRelIso_all, nMuon, Muon_pt, Muon_eta, Muon_miniPFRelIso_all)')
        self.a.Cut('isLeptonPreselected','isLeptonPre')
        self.LEPPRE=self.getNweighted()
        self.AddCutflowColumn(self.LEPPRE,"LEPPRE")

        #at least two fat jet w/ mass > 50
        self.a.Define('FatJetHM_msoftdrop','FatJet_msoftdrop > 50')
        self.a.Cut('nFatJet_HighMass','FatJetHM_msoftdrop.size() > 1')
        self.JETMASS = self.getNweighted()
        self.AddCutflowColumn(self.JETMASS,'JETMASS')

        #high pt fat jets
        self.a.Cut('Jet_pt','FatJet_pt[0] > 300 && FatJet_pt[1] > 200')
        self.JETPT=self.getNweighted()
        self.AddCutflowColumn(self.JETPT,"JETPT")

        #jet eta
        self.a.Cut('eta_cut', 'abs(FatJet_eta[0]) < 2.4 && abs(FatJet_eta[1]) < 2.4') # no super forward jets
	self.JETETA = self.getNweighted()
	self.AddCutflowColumn(self.JETETA, "JETETA")
 
        return self.a.GetActiveNode()

    def GetXsecScale(self):
        lumi = self.config['lumi{}'.format(self.year if 'APV' not in self.year else 16)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        return lumi*xsec/self.a.genEventSumw

 # corrections - used in both snapshots and selection
    def ApplyStandardCorrections(self, post_snapshot=False):
	# first apply corrections for right after snapshot phase
	if post_snapshot:
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

	    self.a = AutoJME.AutoJME(self.a,'FatJet', '20{}'.format(self.year), self.setname if 'Muon' not in self.setname else self.setname[10:])
            #self.a.MakeWeightCols(extraNominal='genWeight' if not self.a.isData else '') #this will have been done separately before the snapshot
            #self.a.MakeWeightCols(extraNominal='' if self.a.isData else 'genWeight*%s'%self.GetXsecScale())

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
                    #self.a.AddCorrection(Correction('TptReweight',corrtype='weight')
	return self.a.GetActiveNode()

    def Snapshot(self, node=None):
        #for now, snapshot without applying any corrections (deal w/ this stuff later)
        startNode = self.a.GetActiveNode()
        if node == None:
            node = self.a.GetActiveNode()
              
        columns = [
        'nJet','nFatJet','FatJet_eta','FatJet_msoftdrop','FatJet_pt','FatJet_phi',
        'FatJet_deepTagMD_HbbvsQCD', 'FatJet_deepTagMD_ZHbbvsQCD',
        'FatJet_deepTagMD_WvsQCD', 'FatJet_deepTag_TvsQCD', 'FatJet_particleNet_HbbvsQCD',
        'FatJet_particleNet_TvsQCD', 'FatJet_particleNet_QCD',
        'FatJet_particleNet_WvsQCD','FatJet_particleNetMD*', 'FatJet_rawFactor', 'FatJet_tau*',
        'FatJet_jetId','nMuon','Muon_*','nElectron','Electron_*','HLT_Mu50','HLT_TkMu*','HLT_IsoMu*','HLT_OldMu100',
        'MET_*','HLT_*MET*','HLT_Photon*','HLT_Ele*','HLT_PFHT*', 'HLT_PFJet*', 'HLT_AK8*','weight__nominal','genWeight',
        'event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','LEPPRE','NJET','JETMASS',
        'JETPT','JETETA','Pileup*','fixedGridRhoFastjetAll','L1PreFiringWeight*'
        ]         
       
	#'FatJet_JES_nom', 
        # append to columns list if not Data
        if not self.a.isData:
            columns.extend(['GenPart_*', 'nGenPart', 'genWeight'])
            
            columns.extend(['FatJet_JES_up','FatJet_JES_down',
                            'FatJet_JER_nom','FatJet_JER_up','FatJet_JER_down',
                            'FatJet_JMS_nom','FatJet_JMS_up','FatJet_JMS_down',
                            'FatJet_JMR_nom','FatJet_JMR_up','FatJet_JMR_down'])
            columns.extend(['Pileup__nom','Pileup__up','Pileup__down'])
                            #'Pdfweight__nom','Pdfweight__up','Pdfweight__down'])
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

        self.NDIJETS = self.getNweighted()
        self.AddCutflowColumn(self.NDIJETS, "NDIJETS")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

        #Define mass-decorrelated Higgs/Wqq taggers
        self.a.Define('Dijet_particleNetMD_HbbvsQCD','Dijet_particleNetMD_Xbb/(Dijet_particleNetMD_Xbb+Dijet_particleNetMD_QCD)')
        self.a.Define('Dijet_particleNetMD_WqqvsQCD','Dijet_particleNetMD_Xqq/(Dijet_particleNetMD_Xqq+Dijet_particleNetMD_QCD)')

        #mark one jet as Higgs and other as Wqq
        self.a.Define('HiggsIdx','Dijet_particleNetMD_HbbvsQCD[0] > Dijet_particleNetMD_HbbvsQCD[1] ? 0 : 1')
        self.a.Define('WqqIdx','HiggsIdx == 0 ? 1 : 0')

        self.a.ObjectFromCollection('Higgs','Dijet','HiggsIdx')
        self.a.ObjectFromCollection('Wqq','Dijet','WqqIdx')

        return self.a.GetActiveNode()

    def SignalLepton(self):

        self.a.Define('kinEleIdx','kinElectron(Electron_pt,Electron_eta,Electron_phi,Higgs_phi,Wqq_phi)')
        self.a.Define('kinMuIdx','kinMuon(Muon_pt,Muon_eta,Muon_phi,Higgs_phi,Wqq_phi)')  
        self.a.Cut('kinLepton_cut','kinEleIdx != -1 || kinMuIdx != -1') #at least one good lepton
        self.a.Define('LeptonType','LeptonIdx(kinEleIdx,kinMuIdx,Electron_pt,Muon_pt)') #picks higher pt signal lepton - output = 0 (lepton is electron) or 1 (lepton is muon)

        self.a.Cut('lepIsolation','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < {} : Muon_miniPFRelIso_all[kinMuIdx] < {}'.format(self.cuts['ELECTRON']['RelIso'],self.cuts['MUON']['RelIso']))
        self.a.Cut('leptonQuality','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]') #"naively" using the standard WP80 boolean for now; eventually must examine more closely efficiency of lepton cuts and find optimal cut value(s) for Electron_mvaFall17V2noIso discriminant

        self.SIGLEP = self.getNweighted()
        self.AddCutflowColumn(self.SIGLEP,'SIGLEP')

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
	self.nLeadW = self.getNweighted()
	self.AddCutflowColumn(self.nLeadW, 'nWqqMassCut')
	
        return self.a.GetActiveNode()

    def ApplySRorCR(self, SRorCR, tagger): #This is getting frequently updated as I play w/ different definitions
	'''
	SRorCR [str] = 'SR' or 'CR', used to generate cutflow information
	tagger [str] = name of tagger to be used 
	Jet tagging criteria:
		SR: Wqq > 0.8
                    Electron_mvaFall17V2noIso_WP80 || Muon_MediumId  
                    Lepton_miniPFRelIso_all < 0.2
		CR: Wqq < 0.8
                    !Electron_mvaFall17V2noIso_WP80 && !Muon_MediumId
                    Lepton_miniPFRelIso_all > 0.2
	'''
	assert(SRorCR=='SR' or SRorCR=='CR')
	# tagger values for each 
	if SRorCR == 'SR':
	    # Signal region
	    self.a.Cut('Wqq_{}_cut_{}'.format(tagger,SRorCR), 'Wqq_{0}_WqqvsQCD > {1}'.format(tagger,self.cuts['JET']['particleNetMD_WqqvsQCD']))
            self.a.Cut('lepIsolation_{}'.format(SRorCR),'LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < {} : Muon_miniPFRelIso_all[kinMuIdx] < {}'.format(self.cuts['ELECTRON']['RelIso'],self.cuts['MUON']['RelIso']))
            self.a.Cut('leptonQuality_{}'.format(SRorCR),'LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')

	else:
	    # Control Region
	    self.a.Cut('Wqq_{}_cut_{}'.format(tagger,SRorCR), 'Wqq_{0}_WqqvsQCD < {1}'.format(tagger,self.cuts['JET']['particleNetMD_WqqvsQCD']))
            self.a.Cut('lepIsolation_{}'.format(SRorCR),'LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] > {} : Muon_miniPFRelIso_all[kinMuIdx] > {}'.format(self.cuts['ELECTRON']['RelIso'],self.cuts['MUON']['RelIso']))
            self.a.Cut('leptonQuality_{}'.format(SRorCR),'LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] == false : Muon_mediumId[kinMuIdx] == false')
	# save cutflow info
	self.nWTag = self.getNweighted()
	self.AddCutflowColumn(self.nWTag, 'nWTag_{}'.format(SRorCR))
	return self.a.GetActiveNode()

    def ApplyPassFail(self, SRorCR, tagger): #also getting changed pretty frequently for now
	'''
	Fail:	Hbb < 0.94
	Pass: 	Hbb > 0.94
	'''
	assert(SRorCR=='SR' or SRorCR=='CR')
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
