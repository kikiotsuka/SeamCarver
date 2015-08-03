import os, pathlib

p = pathlib.Path()

os.chdir('./images')

carvelist = list(p.rglob('*.jpg')) + \
            list(p.rglob('*.jpeg')) + \
            list(p.rglob('*.png')) + \
            list(p.rglob('*.gif'))
carvelist = [str(f) for f in carvelist]

os.chdir('../')

fout = open('toseamcarve.txt', 'w')
for f in carvelist:
    fout.write(str(f) + '\n')
fout.close()
