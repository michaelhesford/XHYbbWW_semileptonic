import glob,os

out = open('condor/distribution_args.txt','w')

for t in ['Signal','Background']:
    if t == 'Signal':
        for f in glob.glob('raw_nano/Signal/*'):
            if os.path.getsize(f) == 0:
                print ('File %s is empty... Skipping.'%(f))
                continue

            # filename will look like: raw_nano/Signal/setname_year.txt
            filename = f.split('/')[-1].split('.')[0]
            nfiles = len(open(f,'r').readlines())
            setname = filename.split('_')[0]
            year = filename.split('_')[1]

            njobs = int(nfiles/2)
            if njobs == 0:	# this occurs when nfiles = 1
	        njobs += 1
            for i in range(1,njobs+1):
  	        out.write('-s %s -y %s -j %s -n %s \n'%(setname,year,i,njobs))

    elif t == 'Background': # Maybe i can figure out a way to do this in one loop...
        for d in glob.glob('raw_nano/Background/*'):
            for f in glob.glob(d+'/*'):
                if os.path.getsize(f) == 0:
                    print ('File %s is empty... Skipping.'%(f))
                    continue

                # filename will look like: raw_nano/Background/setname_year.txt
                filename = f.split('/')[-1].split('.')[0]
                nfiles = len(open(f,'r').readlines())
                setname = filename.split('_')[0]
                year = filename.split('_')[1]

                njobs = int(nfiles/2)
                if njobs == 0:      # this occurs when nfiles = 1
                    njobs += 1
                for i in range(1,njobs+1):
                    out.write('-s %s -y %s -j %s -n %s \n'%(setname,year,i,njobs))

out.close()
