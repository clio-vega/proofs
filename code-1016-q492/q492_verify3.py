import sys; sys.path.insert(0,'/home/clio/projects/scratch')
exec(open('/home/clio/projects/scratch/q492_verify.py').read().split("print('=== CLEAN RUN ===')")[0])
print('=== CLEAN RUN N=3 ===')
run(3,tag='clean')
