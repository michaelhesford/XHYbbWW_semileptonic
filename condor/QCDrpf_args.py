import glob

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--recycle', dest='recycle',
                    action='store_true', default=False,
                    help='Recycle existing files and just plot.')
args = parser.parse_args()

out = open('condor/QCDrpf_args.txt','w')
for f in glob.glob('snapshots/*.txt'):
    # filename format is "setname_year_snapshot.txt"
    filename = f.split('/')[-1].split('.')[0]
    setname = filename.split('_')[0]
    year = filename.split('_')[1]
    #if setname.startswith('Data') or 'XHY' in setname: 
    #     continue
    toSkip = [
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
        'ZHtoWW',
        'WW',
        'WZ',
        'ZZ'
        'ZJetsHT200',
        'ST-shad'
    ] 
    if setname.startswith('Data') or setname in toSkip: continue
    if 'XHY' in setname and setname != 'XHY-3500-1800': continue
    if 'Data' not in setname:
        out.write('-s {} -y {} -v None\n'.format(setname, year))
        systs = ['JES','JER','JMS','JMR','UE']
        for syst in systs:
            for v in ['up','down']:
                out.write('-s {} -y {} -v {}_{}\n'.format(setname, year, syst, v))
        
    else:
        out.write('-s {} -y {} -v None\n'.format(setname, year))

out.close()
