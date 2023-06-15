import glob,os

out = open('condor/MXvsMYstudies_args.txt','w')
for f in glob.glob('snapshots/*.txt'):
    # filename format is "setname_year_snapshot.txt"
    filename = f.split('/')[-1].split('.')[0]
    nfiles = len(open(f,'r').readlines())
    setname = filename.split('_')[0]
    year = filename.split('_')[1]

    # now write to file
    if 'XHY' in filename:
	#for signal, there's only one data file per .txt file 
	njobs = 1
    else:
        njobs = int(nfiles/4)
    if njobs == 0: 	# this occurs when nfiles = 1
	njobs += 1
    for i in range(1, njobs+1):
	out.write('-s {} -y {} -j {} -n {}\n'.format(setname, year, i, njobs))

