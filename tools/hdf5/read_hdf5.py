import numpy as np
import h5py




dir="/data/s5/yhwu/data/TNG/snap/snap_099.0.hdf5"
f=h5py.File(dir,'r')
print f.keys()
print f['Header'].attrs.items()
header=dict(f['Header'].attrs.items())
print f['Parameters']
print f['PartType1'].keys()


print f['PartType1']['Coordinates']
print f['PartType1']['Velocities']
print f['PartType1']['ParticleIDs']
