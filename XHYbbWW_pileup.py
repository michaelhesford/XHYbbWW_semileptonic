import ROOT, glob
from TIMBER.Analyzer import TIMBERPATH, analyzer, Correction
from TIMBER.Tools.AutoPU import MakePU

'''
to be run *before* creating snapshots
'''

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
                        action='store', required=True,
                        help='Setname to process.')
    parser.add_argument('-y', type=str, dest='era',
                        action='store', required=True,
                        help='Year of set (16, 17, 18).')
    parser.add_argument('-j', type=int, dest='ijob',
                        action='store', required=False,
                        default=1, help='Current job to process')
    parser.add_argument('-n', type=int, dest='njobs',
                        action='store', required=False,
                        default=1, help='Total number of jobs')
    args = parser.parse_args()
    setname=args.setname
    era=args.era
    ijob=args.ijob
    njobs=args.njobs

    # setname: DataX, XMASS-YMASS
    fullname = '%s_%s'%(setname,era)

    filename = 'raw_nano/{}_{}.txt'.format(setname,era)

    bkgs = [
        'ttWJets',
        'ttWJetsLep',
        'ttZJets',
        'ttHtoBB',
        'ttHnoBB',
        'WpHbb',
        'WmHbb',
        'ZHbbll',
        'VBFHbb',
        'ggHbb',
        'ggZHbbll',
        'ggHWW',
        'WpHtoWW',
        'WmHtoWW',
        'ZHtoWW'
    ]

    rootFiles=open(filename,'r').readlines()
    #inFile = rootFiles[ijob-1][:-1] if not njobs == 1 else filename #requires we work with all files separately or all together
    inFile = filename
    print(inFile)
    print(rootFiles)
    '''
    if setname in bkgs:
        for i in range(len(rootFiles)):
            rootFiles[i] = rootFiles[i].replace('cmsxrootd.fnal.gov','cms-xrd-global.cern.ch').strip()
        print('USING GLOBAL REDIRECTOR: ',rootFiles)
        a = analyzer(rootFiles)
    else:
        a = analyzer(inFile)
    ''' 
    a = analyzer(inFile)

    # get pointer to histogram
    #hptr = MakePU(a, '20%sUL'%args.era, fullname+'.root')
    hptr = MakePU(a, era, ULflag=True, filename=fullname+'.root')	# update to latest TIMBER version
    hout = hptr.Clone()
    
    out_name = 'XHYbbWW_pileup_{}_{}_{}.root'.format(fullname,ijob,njobs) if njobs > 1 else 'XHYbbWW_pileup_{}.root'.format(fullname)
    out = ROOT.TFile.Open(out_name, 'RECREATE')
    out.WriteTObject(hout, fullname)
    out.Close()
