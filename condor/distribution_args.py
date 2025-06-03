import glob,os

out = open('condor/distribution_args.txt','w')

for f in glob.glob('snapshots/*snapshot.txt'):
    if os.path.getsize(f) == 0 or f.split('/')[-1].startswith('Data'):
        print ('File %s is empty or JetHT data... Skipping.'%(f))
        continue
    signals = ['XHY-4000-2500','XHY-3500-1800','XHY-3000-1400','XHY-2000-1000','XHY-1600-700','XHY-1000-500']
    # filename will look like: snapshots/setname_year_snapshot.txt
    filename = f.split('/')[-1].split('.')[0]
    nfiles = len(open(f,'r').readlines())
    setname = filename.split('_')[0]
    if 'XHY' in setname and setname not in signals:
        continue
    year = filename.split('_')[1]
    out.write('-s %s -y %s -j 1 -n 1\n'%(setname,year))

out.close()
