import glob,os

out = open('broken_template_jobs.txt','w')
cluster = '60749928'
for f in glob.glob('logs/output_{}*.stderr'.format(cluster)):
    lines = open(f,'r').readlines()
    for l in lines:
        if 'Traceback' in l:
            job = f.split('_')[-1].split('.')[0]
            out.write('%s\n'%job)
out.close() 
           
