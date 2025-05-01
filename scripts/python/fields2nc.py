#!/usr/bin/python3

import netCDF4 as nc
import numpy as np
import struct
import sys

# etype = ">" # big-endian
etype = "<" # little-endian

dtype = "d" # floating-point number
sizeofdata = 8 # in bytes

# dtype = "f" # floating-point number
# sizeofdata = 4 # in bytes

# dtype = 'B' # unsigned character, for gate files
# sizeofdata = 1 # in bytes

nx = 0 # number of points in Ox; if 0, then search tlab.ini
ny = 0 # number of points in Oy; if 0, then search tlab.ini
nz = 0 # number of points in Oz; if 0, then search tlab.ini

# do not edit below this line

# getting grid size from tlab.ini, if necessary
if ( nx == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("imax="):
            nx = int(line.split("=",1)[1])
            break

if ( ny == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("jmax="):
            ny = int(line.split("=",1)[1])
            break

if ( nz == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("kmax="):
            nz = int(line.split("=",1)[1])
            break

print("Grid size is {}x{}x{}.".format(nx,ny,nz))

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [3d,xy,xz,yz] list-of-files.")
    quit()

datatype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

if   ( datatype == '3d' ):          # full 3d field
    sizeofmask = 0
elif ( datatype == 'xy' ):          # plane of constant z
    sizeofmask = len(str(nz))
elif ( datatype == 'xz' ):          # plane of constant y
    sizeofmask = len(str(ny))
elif ( datatype == 'yz' ):          # plane of constant x
    sizeofmask = len(str(nx))
else:
    print("Usage: python $0 [3d,xy,xz,yz] list-of-files.")
    quit()

# read grid information
fin = open('grid', 'rb')
#
fin.seek(56,0)
raw = fin.read(nx*8)
x = struct.unpack(etype+'{}d'.format(nx), raw)
#
fin.seek(8,1)
raw = fin.read(ny*8)
y = struct.unpack(etype+'{}d'.format(ny), raw)
#
fin.seek(8,1)
raw = fin.read(nz*8)
z = struct.unpack(etype+'{}d'.format(nz), raw)
#
fin.close()

# crop grid data in case of planes
if   ( datatype == 'xy' ):
    nz = 1
    z = np.array([0.0], dtype=np.float32)
elif ( datatype == 'xz' ):
    ny = 1
    y = np.array([0.0], dtype=np.float32)
elif ( datatype == 'yz' ):
    nx = 1
    x = np.array([0.0], dtype=np.float32)

# process files
for file in setoffiles:
    # reading data
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    fin.seek(-nx*ny*nz*sizeofdata, 2)
    raw = fin.read()
    a = np.array(struct.unpack((etype+'{}'+dtype).format(int(nx*ny*nz)), raw))
    a = a.reshape((nz,ny,nx))
    fin.close()

    # creating netcdf
    file_dst = nc.Dataset(file+'.nc', 'w')

    # create dimensions for destiny nc-file
    file_dst.createDimension('x',len(x))
    file_dst.createDimension('y',len(y))
    file_dst.createDimension('z',len(z))

    # create and write independent variables in destiny nc-file using single precision
    x_dst = file_dst.createVariable('x', 'f4', ('x',))
    y_dst = file_dst.createVariable('y', 'f4', ('y',))
    z_dst = file_dst.createVariable('z', 'f4', ('z',))
    x_dst[:] = x[:]
    y_dst[:] = y[:]
    z_dst[:] = z[:]

    var_dst = file_dst.createVariable(file, 'f4', ('z','y','x',))
    var_dst[:] = a[:]

    file_dst.close()
