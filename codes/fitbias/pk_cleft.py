import numpy as np
import sys, os

path_to_cleft = '/home/chirag/Research/codes/CLEFT_GSM/ps_py3_package/'
sys.path.append(path_to_cleft)

import cleftpool as cpool

pkfile = '~/Research/Projects/qso_kappacmb/data/pklin_ds14.dat'
pkfile = '/home/chirag/Research/Projects/qso_kappacmb/data/pklin_ds14.dat'
pklin = np.loadtxt(pkfile).T

print(pklin.shape)

#if os.path.isfile('../../data-generated/q_kernel.dat'):
#    qfile = '../../data-generated/q_kernel.dat'
#else: qfile = None
#if os.path.isfile('../../data-generated/r_kernel.dat'):
#    rfile = '../../data-generated/r_kernel.dat'
#else: rfile = None

qfile, rfile = None, None

#cl = cpool.CLEFT(k=pklin[0], p=pklin[1], npool=32, qfile=qfile, rfile=rfile)
cl = cpool.CLEFT(pfile = pkfile, npool=8, qfile=qfile, rfile=rfile)

if qfile is None: cpool.save_qkernel(cl, '../../data-generated/q_kernel.dat')
if rfile is None: cpool.save_rkernel(cl, '../../data-generated/r_kernel.dat')
cpool.save_qfunc(cl, '../../data-generated/q_func.dat')

for zz in np.arange(0, 2.5, 0.1):
    print('For z = %0.2f'%zz)
    pk= cpool.make_table(cl, nk=100, npool=32, M=0.295, z=zz)
    cpool.save_pk(pk, '../../data-generated/pkcleft_zz%03d.dat'%(zz*100))
