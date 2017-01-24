import numpy as np

infile = '/Users/alexialewis/research/galbase/code/adam/MyAllSkyTable_akleroy.csv'
newcommandfile = '/Users/alexialewis/research/galbase/ngc2976_wget_commands.txt'


testfiles = ['AIS_74_sg02-fd', 'AIS_74_sg04-fd', 'AIS_74_sg05-fd', 'AIS_74_sg08-fd', 'AIS_74_sg92-fd', 'AIS_74_sg93-fd', 'GI1_071001_M81-fd', 'GI2_024001_NGC3077_Stream-fd', 'GI2_024002_NGC2976_stream-fd', 'GI3_061016_KK77-fd', 'NGA_NGC2976-fd']


f = open(infile, 'r')
lines = f.readlines()
f.close()

goodlines = lines[1:]
allfiles = [goodlines[i].split(' ')[-1].strip('"').strip('\n').split('/')[-1].strip('"').split('.')[0].rstrip('-flags').rstrip('-int') for i in range(len(goodlines))]

blah = np.where(np.in1d(allfiles, testfiles))
names = np.asarray(allfiles)[blah[0]]


from collections import defaultdict
D = defaultdict(list)
for i,item in enumerate(names):
    D[item].append(i)
D = {k:v for k,v in D.items() if len(v)>1}

inds = []
for i in range(len(names)):
    inds.append(D[names[i]][0])

unique_inds = np.unique(inds)
goodinds = blah[0][unique_inds]

filestarts = np.asarray(goodlines)[goodinds]
filestarts = [f.replace('-int', '-flags')for f in filestarts]

newlines1 = [f.replace('-flags', '-cnt')for f in filestarts]
newlines2 = [f.replace('-flags', '-exp')for f in filestarts]
newlines3 = [f.replace('-flags', '-rrhr')for f in filestarts]
newlines4 = [f.replace('-flags', '-wt')for f in filestarts]

g = open(newcommandfile, 'w')
for nl in [newlines1, newlines2, newlines3, newlines4]:
    g.writelines(nl)

g.close()
