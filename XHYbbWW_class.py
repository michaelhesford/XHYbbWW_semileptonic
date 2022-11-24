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
        if inputfile.endswith('.txt'):
            infiles = SplitUp(inputfile, njobs)[ijob-1]
        else:
            infiles = inputfile

        if inputfile.endswith('.txt'):
	    # we're looking at signal
            if 'XYH' in inputfile:
		# format: (raw_nano/XHY-<XMASS>-<YMASS>_year.txt
                prefix = inputfile.split('/')[-1].split('-')   # [XHY, <XMASS>, <YMASS>_year.txt]
                self.year = prefix[2].split('_')[1].split('.')[0]
                self.setname = ('MX_' + prefix[1] + '_MY_' + prefix[2].split('_')[0])  # MX_XMASS_MY_YMASS
	        # create an analyzer module with the proper multiSampleStr argument
                self.a = analyzer(infiles,multiSampleStr=self.setname[3])
		# ensure we're working with the proper YMass
                self.a.Cut('CorrectMass', 'GenModel_YMass_{} == 1'.format(self.setname[3]))
            else:
		# Looking at data, format: (raw_nano/setname_era.txt)
                self.setname = inputfile.split('/')[-1].split('_')[0]
                self.year= inputfile.split('/')[-1].split('_')[1].split('.')[0]
                self.a = analyzer(infiles)
        else:	# not likely to encounter this for this analysis
            self.setname = inputfile.split('/')[-1].split('_')[1]
            self.a = analyzer(infiles)
        
#       self.year = str(year)
        self.ijob = ijob
        self.njobs = njobs
        
	# get config from JSON
        self.config = OpenJSON('XHYbbWWconfig.json')
        self.cuts = self.config['CUTS']

	# triggers for various years
        self.trigs = {
            #'16':['HLT_PFHT800','HLT_PFHT900'],
            #'17':['HLT_PFHT1050','HLT_AK8PFJet500'],
            #'18':['HLT_AK8PFJet400_TrimMass30','HLT_AK8PFHT850_TrimMass50','HLT_PFHT1050']
            16:['HLT_PFHT800','HLT_PFHT900'],
            17:["HLT_PFHT1050","HLT_AK8PFJet500","HLT_AK8PFHT750_TrimMass50","HLT_AK8PFHT800_TrimMass50","HLT_AK8PFJet400_TrimMass30"],
            18:['HLT_AK8PFJet400_TrimMass30','HLT_AK8PFHT850_TrimMass50','HLT_PFHT1050']
        }

        # check if data or sim
        if 'Data' in inputfile:
            self.a.isData = True
        else:
            self.a.isData = False
    
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

    # CompileCPP('~/public/XHYbbWW_analysis/CMSSW_11_1_4/XHYbbWW_semileptonic') # compile C++ code

    def KinematicSnap(self):
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


        #Lepton number/pt/eta cuts

        self.a.Cut('nLepton','nElectron > 0 || nMuon > 0')
        self.NLEPTON=self.getNweighted()
        self.AddCutflowColumn(self.NLEPTON,"NLEPTON")

        self.a.Cut('Lepton_pt','Electron_pt[0] > 10 || Muon_pt[0] > 10')
        self.LEPPT=self.getNweighted()
        self.AddCutflowColumn(self.LEPPT,'LEPPT')

        self.a.Cut('Lepton_eta', 'abs(Electron_eta[0]) < 2.4 || abs(Muon_eta[0]) < 2.5')
        self.LEPETA=self.getNweighted()
        self.AddCutflowColumn(self.LEPETA,'LEPETA')
'''
        self.a.Cut('Lepton_mass','Electron_mass[0] > ## || Muon_mass[0] > ##') # No need to do mass cuts on leptons, right?
        self.LEPMASS=self.getNweighted()
        self.AddCutflowColumn(self.LEPMASS,'LEPMASS')
'''


'''
    def LeptonISL(self,param = 'miniPFRelIso_all'):
        Cuts=[0.3,0.1] #[miniPFRelIso_all,FatJet_lsf3] note: find a value that makes sense for FatJet_lsf3
        if param == 'miniPFRelIso_all':
            self.a.Cut('LepISL','Electron_{0}[0] < {1} || Muon_{0}[0] < {1}'.format(param,Cuts[0]))
        else: # using FatJet_lsf3
            self.a.Cut('LepISL','FatJet_lsf3[0] < {0}'.format(cuts[1]))
        self.LEPISL=self.getNweighted()
        self.AddCutflowColumn(self.LEPISL,'LEPISL')
# I doubt it make sense to impose an isolation requirement on the highest pT lepton. Consolidate everything into one C++ function to pick good leptons
'''

    def Dijets(self):
        self.a.Define('DijetIdxs','PickDijets(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_msoftdrop)')
        self.a.Cut('DijetsExist','DijetIdxs[0] > -1 && DijetIdxs[1] > -1')
        self.NDIJETS = self.getNweighted()
        self.AddCutflowColumn(self.NDIJETS, "NDIJETS")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

    def GoodLeptons(self):
        self.a.Define('ElectronIdx', 'PickIsolatedLeptons(FatJet_phi,Electron_phi,Electron_pt,Electron_eta,Electron_dxy,Electron_dz,Electron_miniPFRelIso_all)')
        self.a.Define('MuonIdx', 'PickIsolatedLeptons(FatJet_phi,Muon_phi,Muon_pt,Muon_eta,Muon_dxy,Muon_dz,Muon_miniPFRelIso_all)')

        self.a.Cut('LeptonExists','ElectronIdx[0] != -1 || MuonIdx[0] != -1')
        self.GOODLEP = self.getNweighted()
        self.AddCutflowColumn(self.GOODLEP, "GOODLEP")

        self.a.Define('Lepton_type','MuonIdx[0] == -1 ? 0 : ElectronIdx[0] == -1 ? 1 : Electron_pt[ElectronIdx[0]] > Muon_pt[MuonIdx[0]] ? 0 : 1)') # From now on, 0 = Electron and 1 = Muon
        self.a.SubCollection('LeptonISL','Lepton_type == 0? Electron : Muon','Lepton_type == 0? ElectronIdx : MuonIdx',useTake=True) #No clue if this works 
   
