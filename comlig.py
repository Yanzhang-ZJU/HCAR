"""
MIT License

Copyright (c) 2023 Kun, Xi @ Zhejiang University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import mdtraj as md
import numpy as np
import gromacs.run
import gromacs.tools
import gromacs

def _format_83(f):
	"""Format a single float into a string of width 8, with ideally 3 decimal
	places of precision. If the number is a little too large, we can
	gracefully degrade the precision by lopping off some of the decimal
	places. If it's much too large, we throw a ValueError"""
	if -999.999 < f < 9999.999:
		return '%8.3f' % f
	if -9999999 < f < 99999999:
		return ('%8.3f' % f)[:8]
	raise ValueError('coordinate "%s" could not be represented '
		'in a width-8 field' % f)

xtcName = 'traj.xtc'
topName = 'temp.gro'
comIndName = 'index.ndx'

cInd = np.loadtxt(comIndName, dtype=np.int32)

ori = md.load(xtcName, top=topName)
align = ori.superpose(ori, 0)

fInd = open(comIndName, 'r')
fIndr = fInd.readlines()

fw = open('comlig.pdb', 'w')

posIndex = 1
for i in range(len(align)): 
    xc = 0.0
    yc = 0.0
    zc = 0.0
    
    line = fIndr[0].split()
    for j in range(len(line)):
        xc += float( align.xyz[i, int(line[j])-1, :][0] )/len(line)*10
        yc += float( align.xyz[i, int(line[j])-1, :][1] )/len(line)*10
        zc += float( align.xyz[i, int(line[j])-1, :][2] )/len(line)*10
    line = "ATOM  %5d  H  HEX%5d    %s%s%s  1.00  0.00           C" % (posIndex, posIndex, _format_83( xc ),
							_format_83( yc ), _format_83( zc ) )
    posIndex += 1
    print(line, file=fw)
print('TER', file=fw)
    
    

fw.close()
