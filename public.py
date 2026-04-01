import sys
import os
import numpy as np
from astropy.table import Table, setdiff, Column, vstack

###############################################################################
workdir='/data/munka/SBNAF/irdatabase'

os.chdir(workdir)

infile='hso_out.csv'
outfile='public_hso.csv'
reffile='public_ref.dat'

inp=Table.read(infile,data_start=2,format='ascii.csv')
ref=Table.read(reffile,data_start=1,format='ascii.csv')

tt=Table()
i=0
for reference in ref['ref']:
    t=inp[np.where(inp['documents_references'] == reference)]
    print(reference,len(t))
    if i == 0:
        tt=t
    else:
        tt=(vstack([tt, t]))
    i=i+1
tt.write(outfile, format='ascii.csv',overwrite=True)
