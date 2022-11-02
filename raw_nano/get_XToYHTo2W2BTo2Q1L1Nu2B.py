rom TIMBER.Tools.Common import ExecuteCmd
import sys

# dictionary containing the year as keys, with each value being another dictionary of {sampleName: DAS location}
das = {
    '16': {
	'XHY': '/NMSSM_XToYHTo2W2BTo2Q1L1Nu2B_MX-{}_MY-{}_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'ttbar-allhad': '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'ttbar-semilep': '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'QCDHT700': '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'QCDHT1000': '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'QCDHT1500': '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'QCDHT2000': '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
        'ZJetsHT200': '/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'ZJetsHT400': '/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'ZJetsHT600': '/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'ZJetsHT800': '/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'WJetsHT200': '/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'WJetsHT400': '/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'WJetsHT600': '/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'WJetsHT800': '/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2/NANOAODSIM',
        'DataF':  '/JetHT/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataG':  '/JetHT/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataH':  '/JetHT/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataF': '/SingleMuon/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataG': '/SingleMuon/Run2016G-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataH': '/SingleMuon/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD'
    },
    '16APV': {
	'XHY': '/NMSSM_XToYHTo2W2BTo2Q1L1Nu2B_MX-{}_MY-{}_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
        'ttbar-allhad': '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
        'ttbar-semilep': '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
        'QCDHT700': '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
        'QCDHT1000': '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',        'QCDHT1500': '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',        'QCDHT2000': '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
        'ZJetsHT200': '/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'ZJetsHT400': '/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'ZJetsHT600': '/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'ZJetsHT800': '/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'WJetsHT200': '/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'WJetsHT400': '/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'WJetsHT600': '/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'WJetsHT800': '/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2/NANOAODSIM',
        'DataB':  '/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'DataC':  '/JetHT/Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'DataD':  '/JetHT/Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'DataE':  '/JetHT/Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'DataF':  '/JetHT/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SingleMuonDataB': '/SingleMuon/Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SingleMuonDataC': '/SingleMuon/Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SingleMuonDataD': '/SingleMuon/Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SingleMuonDataE': '/SingleMuon/Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SingleMuonDataF': '/SingleMuon/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD'
    },
    '17': {
	'XHY': '/NMSSM_XToYHTo2W2BTo2Q1L1Nu2B_MX-1000_MY-100_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'ttbar-allhad': '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'ttbar-semilep': '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'QCDHT700': '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'QCDHT1000': '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'QCDHT1500': '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'QCDHT2000': '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM',
        'ZJetsHT200': '/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'ZJetsHT400': '/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'ZJetsHT600': '/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'ZJetsHT800': '/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'WJetsHT200': '/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'WJetsHT400': '/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'WJetsHT600': '/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'WJetsHT800': '/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM',
        'DataB': '/JetHT/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataC': '/JetHT/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataD': '/JetHT/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataE': '/JetHT/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataF': '/JetHT/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataB': '/SingleMuon/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataC': '/SingleMuon/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataD': '/SingleMuon/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataE': '/SingleMuon/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataF': '/SingleMuon/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD'
    },
    '18': {
	'XHY': '/NMSSM_XToYHTo2W2BTo2Q1L1Nu2B_MX-1000_MY-125_TuneCP5_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
        'ttbar-allhad': '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
        'ttbar-semilep': '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
        'QCDHT700': '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'QCDHT1000': '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
        'QCDHT1500': '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
        'QCDHT2000': '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
        'ZJetsHT200': '/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'ZJetsHT400': '/ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'ZJetsHT600': '/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'ZJetsHT800': '/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'WJetsHT200': '/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'WJetsHT400': '/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'WJetsHT600': '/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'WJetsHT800': '/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM',
        'DataA': '/JetHT/Run2018A-UL2018_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'DataB': '/JetHT/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataC': '/JetHT/Run2018C-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'DataD': '/JetHT/Run2018D-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SingleMuonDataA': '/SingleMuon/Run2018A-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD',
        'SingleMuonDataB': '/SingleMuon/Run2018B-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD',
        'SingleMuonDataC': '/SingleMuon/Run2018C-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD',
        'SingleMuonDataD': '/SingleMuon/Run2018D-UL2018_MiniAODv2_NanoAODv9_GT36-v1/NANOAOD'
    }
}

def GetFiles(das_name, setname, year):
    '''
    das_name (str) = name of the dataset on DAS for query
    setname  (str) = name of the dataset for output file
    year     (str) = year of dataset (16, 16PAV, 17, 18)

    The CMS Data Aggregation System (DAS - https://cmsweb.cern.ch/das/) is a service used to
    aggregate all of the data and Monte Carlo files centrally. We will perform DAS queries to 
    get the remote location of all MC signal/background files, as well as data from the JetHT
    and SingleMuon datasets. The JetHT dataset contains events collected at the CMS detector 
    using hadronic jet-based triggers and it is in this data that you'd perform the hadronic
    search. The SingleMuon dataset contains data collected by muon triggers. This dataset is 
    orthogonal to (statistically independent from) the JetHT dataset and can therefore be used 
    to measure our trigger efficiencies across the different years without worrying about 
    biasing due to the use of jet-based triggers. 
    '''
    # first, execute the command-line DAS query and redirect (>) the results into a temp file which we'll read
    ExecuteCmd('dasgoclient -query "file dataset={}" > {}_{}_temp.txt'.format(das_name, setname, year))
    # read the temp file
    f = open('{}_{}_temp.txt'.format(setname, year), 'r')
    # open a new file to write to, with the remote redirector in front of each line
    fout = open('{}_{}.txt'.format(setname, year), 'w')
    for l in f.readlines():
	# add redirector so that you (and TIMBER analyzer) can open the files remotely 
	fout.write('root://cmsxrootd.fnal.gov/'+l)
    f.close()
    fout.close()
    # clean up
    ExecuteCmd('rm {}_{}_temp.txt'.format(setname, year))

if __name__=='__main__':
    # We'll have the program save a LaTeX formatted table of the DAS locations for all the datasets, in case you need that later
    latex_lines = {k:[] for k in das.keys()}
    # loop over our DAS dataset dictionary
    for year in das.keys():
	for setname in das[year].keys():
	    if 'XHY' in setname:    # we're looking at signal
		for mX in [240, 280, 300, 320, 360, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2500, 2600, 2800, 3000, 3500, 4000]:
		    for mY in [60, 70, 80, 90, 100, 125,150, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2500, 2600, 2800]:
			if mY < mX+125:		# Y mass must be less than X mass plus Higgs mass
			    # the XHY key in the das sub-dictionary has {} format placeholders for use right now
			    das_name = das[year][setname].format(mX, mY)
			    # let's modify the setname so it comes out looking nicer in the txt file: XHY-XMASS-YMASS
			    setname_mod = '{}-{}-{}'.format(setname, mX, mY)
			    # now let's get the files and save them all to their respective .txt file
			    GetFiles(das_name, setname_mod, year)
			    # append to our latex output
			    latex_lines[year].append('| %s | %s |'%(setname_mod.replace('-', ' ')+' GeV',das_name))
	    # otherwise, we're looking at data or SM MC bkgs
	    else:
		GetFiles(das[year][setname], setname, year)
		latex_lines[year].append('| %s | %s |'%(setname,das[year][setname]))
    # print latex lines to stdout
    for y in sorted(latex_lines.keys()):
	print('\n20%s'%y)
	print('| Setname | DAS location |')
	print('|---------|--------------|')
	for l in latex_lines[y]: print(l)


