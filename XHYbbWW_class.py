import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME

AutoJME.AK8collection = 'Trijet'

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
        # format: (raw_nano/data_type/setname_year.txt
        infiles = SplitUp(inputfile,njobs)[ijob-1]
        self.setname=inputfile.split('/')[-1].split('_')[0]
        self.year=inputfile.split('/')[-1].split('_')[-1].split('.')[0]
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
 
        self.ijob = ijob
        self.njobs = njobs
     
	# Cuts from config file
        self.config = OpenJSON('XHYbbWW_config.json')
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

    def Preselection(self):
        self.NPROC = self.getNweighted()
        self.AddCutflowColumn(self.NPROC, "NPROC")

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
        MET_filters = self.a.GetFlagString(flags)	# string valid (existing in RDataFrame node) flags together w logical and
        self.a.Cut('flags', MET_filters)
        self.NFLAGS = self.getNweighted()
        self.AddCutflowColumn(self.NFLAGS, "NFLAGS")

        #Lepton number/pt/eta cuts

        self.a.Cut('nLepton','nElectron > 1 && nMuon > 1')
        self.NLEPTON=self.getNweighted()
        self.AddCutflowColumn(self.NLEPTON,"NLEPTON")

        self.a.Cut('Lepton_pt','Electron_pt[0] > 10 || Muon_pt[0] > 10')
        self.LEPPT=self.getNweighted()
        self.AddCutflowColumn(self.LEPPT,"LEPPT")

        self.a.Cut('Lepton_eta', 'abs(Electron_eta[0]) < 2.4 || abs(Muon_eta[0]) < 2.5')
        self.LEPETA=self.getNweighted()
        self.AddCutflowColumn(self.LEPETA,"LEPETA")

        #FatJet number/pt/eta/mass cuts

        self.a.Cut('nFatJet','nFatJet > 1')
        self.NJETS=self.getNweighted()
        self.AddCutflowColumn(self.NJETS,"NJETS")
 
        self.a.Cut('Jet_pt','FatJet_pt[0] > 400 && FatJet_pt[1] > 200')
        self.JETPT=self.getNweighted()
        self.AddCutflowColumn(self.JETPT,"JETPT")
        self.a.Cut('Jet_eta', 'abs(FatJet_eta[0]) < 2.4 && abs(FatJet_eta[1]) < 2.4')
        self.JETETA=self.getNweighted()
        self.AddCutflowColumn(self.JETETA,"JETETA")

        self.a.Cut('Jet_mass', 'FatJet_msoftdrop[0] > 50 && FatJet_msoftdrop[1] > 40')
        self.JETMASS=self.getNweighted()
        self.AddCutflowColumn(self.JETMASS,"JETMASS")

        return self.a.GetActiveNode()

    def Snapshot(self, node=None):
        startNode = self.a.GetActiveNode()
        if node == None:
            node = self.a.GetActiveNode()
        
        columns = [
        'FatJet_J*',	# this will collect all the JME variations created during snapshotting and used in selection 
        'FatJet_lsf3','FatJet_eta','FatJet_msoftdrop','FatJet_pt','FatJet_phi',
        'FatJet_deepTagMD_HbbvsQCD', 'FatJet_deepTagMD_ZHbbvsQCD',
        'FatJet_deepTagMD_WvsQCD', 'FatJet_deepTag_TvsQCD', 'FatJet_particleNet_HbbvsQCD',
        'FatJet_particleNet_TvsQCD', 'FatJet_particleNetMD.*', 'FatJet_rawFactor', 'FatJet_tau*',
        'FatJet_jetId', 'nFatJet', 'FatJet_JES_nom','FatJet_particleNetMD_Xqq',
        'FatJet_particleNetMD_Xcc', 'FatJet_particleNet_QCD','FatJet_particleNet_WvsQCD',
        'Electron_eta','Electron_pt','Electron_mass','Electron_phi','Electron_miniPFRelIso_all'
        'Electron_dxy','Electron_dz','Electron_sip3d','Electron_hoe','Electron_eInvMinusPInv','Electron_deltaEtaSC','Electron_sieie',
        'Muon_eta','Muon_pt','Muon_mass','Muon_phi','Muon_miniPFRelIso_all','Muon_dxy','Muon_dz, Muon_sip3d',
        'Muon_mediumPromptId','Muon_segmentComp','Muon_isGlobal','HLT_PFHT.*', 'HLT_PFJet.*', 'HLT_AK8.*', 'HLT_Mu50',
        'event', 'eventWeight', 'luminosityBlock', 'run','NPROC', 'NJETS', 'NPT', 'NETA', 'NMSD']
        
        # append to columns list if not Data
        if not self.a.isData:
            columns.extend(['GenPart_.*', 'nGenPart', 'genWeight', 'GenModel*'])
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
        self.a.Snapshot(columns, 'HWWsnapshot_{}_{}_{}of{}.root'.format(self.setname,self.year,self.ijob,self.njobs),'Events', openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

    def SignalLepton(self):

        self.a.Define('ElectronIdx','LeptonPre(Electron_pt,Electron_eta,Electron_miniPFRelIso_all)')
        self.a.Define('MuonIdx','LeptonPre(Muon_pt,Muon_eta,Muon_miniPFRelIso_all)')
        self.a.Define('LeptonIdxs','LeptonIdx(ElectronIdx,MuonIdx,Electron_pt,Muon_pt)')
        self.a.Cut('Lepton_cut','LeptonIdxs[0] != -1 || LeptonIdxs[1] != -1') 

        self.SIGLEP = self.getNweighted()
        self.AddCutflowColumn(self.SIGLEP,'SIGLEP')

        return self.a.GetActiveNode()

    def Dijets(self):

        self.a.Define('DijetIdxs','PickDijets(FatJet_eta,FatJet_phi,FatJet_msoftdrop,LeptonIdxs,Electron_phi,Muon_phi)')
        self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')
       
        self.NDIJETS = self.getNweighted()
        self.AddCutflowColumn(self.NDIJETS, "NDIJETS")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

        return self.a.GetActiveNode()

    def JetSelection(self,mW,mH,Wtag,Htag):
    # eg: mH = [100,150], mW = [60,100], Htag = 0.8, Wtag = 0.8 
        self.a.Define('Dijet_particleNetMD_HbbvsQCD','Dijet_particleNetMD_Xbb/(Dijet_particleNetMD_Xbb+Dijet_particleNetMD_QCD)') #mass-decorrelated H tagger
        self.a.Define('Dijet_particleNetMD_WqqvsQCD','Dijet_particleNetMD_Xqq/(Dijet_particleNetMD_Xqq+Dijet_particleNetMD_QCD)') #do same for W

        self.a.Cut('Mass_W','Dijet_msoftdrop[0] > {} && Dijet_msoftdrop[0] < {}'.format(mW[0],mW[1]))
        self.a.Cut('Mass_H','Dijet_msoftdrop[1] > {} && Dijet_msoftdrop[1] < {}'.format(mH[0],mH[1]))

        self.a.Cut('Wtag','Dijet_particleNetMD_WqqvsQCD[0] > {}'.format(Wtag))
        self.a.Cut('Htag','Dijet_particleNetMD_HbbvsQCD[1] > {}'.format(Htag))

        self.JETSELECTION = self.getNweighted()
        self.AddCutflowColumn(self.JETSELECTION,'JETSELECTION')

        return self.a.GetActiveNode()
