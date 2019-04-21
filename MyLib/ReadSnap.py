'''
Standalone execution: ./snapread.py SNAPSHOT
'''

import sys
import time
import numpy as np
import pdb
#pdb.set_trace()

# For reading and storing the header.
class header:
    def __init__(self, snapshot):
        #initialize_block(snapshot)
        read_dummy(snapshot,1)
        self.n_part = np.fromfile(snapshot, 'int32', 6)
        self.mass = np.fromfile(snapshot, 'float64', 6)
        self.time = np.fromfile(snapshot, 'float64', 1)
        self.redshift = np.fromfile(snapshot, 'float64', 1)
        self.flag_sfr = np.fromfile(snapshot, 'int32', 1)
        self.flag_feedback = np.fromfile(snapshot, 'int32', 1)
        self.n_part_total = np.fromfile(snapshot, 'int32', 6)
        self.flag_cooling = np.fromfile(snapshot, 'int32', 1)
        self.num_files = np.fromfile(snapshot, 'int32', 1)
        self.box_size = np.fromfile(snapshot, 'float64', 1)
        self.omega0 = np.fromfile(snapshot, 'float64', 1)
        self.omega_lambda = np.fromfile(snapshot, 'float64', 1)
        self.hubble_param = np.fromfile(snapshot, 'float64', 1)
        self.flag_stellarage =np.fromfile(snapshot,'int32',1)
        self.flag_metals=np.fromfile(snapshot,'int32',1)
        self.n_part_total_HW=np.fromfile(snapshot,'int32',6)
        if self.n_part_total[1]>0:
            self.n_part_TOTAL=int(bin(self.n_part_total_HW[0]).replace('0b','')+bin(self.n_part_total[1]).replace('0b',''),2)
        else:
            self.n_part_TOTAL=int(bin(self.n_part_total[2]).replace('0b','')+(bin(self.n_part_total[1] & 0b11111111111111111111111111111111).replace('0b','')).zfill(32),2)
        self.fill = np.fromfile(snapshot, 'int8', 64)
        read_dummy(snapshot, 1)


# These declared variables are for clarity purposes.
class particle:
    pos = None
    vel = None
    mass = 0.0
    U = None
    rho = None
    smoothing = None
    ID = None
def snapshot_list(basename="snapshot_000",files="10"):
	snapshot_list=[]
	if files==1:
		snapshot_list.append(basename)
		return snapshot_list
	else:
		for i in range(files):
			snapshot_list.append(basename+'.'+str(i))
	return snapshot_list
	

# For reading the leading and trailing ints that are present in each block.
def read_dummy(snapshot, n_dummies):
    for i in np.arange(n_dummies):
        dummy = np.fromfile(snapshot, 'int32', 1)


def initialize_block(snapshot):
    read_dummy(snapshot, 1)
    block_ID = np.fromfile(snapshot, 'int8', 4)
    read_dummy(snapshot, 3)
    return ''.join([chr(i) for i in block_ID])

    
def read_data(snapshot):
    h=header(snapshot)
    p_list = []
    n_part = sum(h.n_part)
    range_ = np.arange(n_part)
    
    # Positions
    #initialize_block(snapshot)
    read_dummy(snapshot,1)
    for i in range_:
        p = particle()
        p.pos = np.fromfile(snapshot, 'float32', 3)
        p_list.append(p)
    read_dummy(snapshot, 2)
    return p_list
    # Velocities
    #initialize_block(snapshot)
    for i in range_:
        p_list[i].vel = np.fromfile(snapshot, 'float32', 3)
    read_dummy(snapshot, 2)

    # IDs
    #initialize_block(snapshot)
    for i in range_:
        p_list[i].ID = np.fromfile(snapshot, 'int32', 1)[0]
    read_dummy(snapshot, 2)

    # Variable masses, which are read in case the mass of the
    # particle of type 'i' is declared as 0, in the header.
    cur = 0
    read_something = 0
    for i in np.arange(6):
        if(h.mass[i] != 0):
            cur += h.n_part[i]
            continue
        else:
            if(read_something == 0):
                read_something = 1
                #initialize_block(snapshot)
            for j in np.arange(h.n_part[i]):
                p_list[cur].mass = np.fromfile(snapshot, 'float32', 1)[0]
                cur += 1

    # The variable masses block might not exist.
    if(read_something):
        read_dummy(snapshot, 1)
    
    # Blocks related to the internal energies, densities and smoothing
    # lengths of the gas particles, in case there is any.
    if(h.n_part[0] > 0):
        range_ = np.arange(h.n_part[0])

        # First the energies
        #initialize_block(snapshot)
        for i in range_:
            p_list[i].U = np.fromfile(snapshot, 'float32', 1)[0]
        read_dummy(snapshot, 2)

        # Then the densities
        #initialize_block(snapshot)
        for i in range_:
            p_list[i].rho = np.fromfile(snapshot, 'float32', 1)[0]
        read_dummy(snapshot, 2)

        # And the smoothing lengths
        #initialize_block(snapshot)
        for i in range_:
            p_list[i].smoothing = np.fromfile(snapshot, 'float32', 1)[0]
        read_dummy(snapshot, 2)
    
    # There are some optional blocks that are yet to be implemented. So,
    # in case any of them is present, send a warning.
    chunk = snapshot.read()
    if chunk:
        print "There still were things to be read..."
    return p_list

#An example of simple application of the script.

def main():

    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    filename="/mnt/d1/P512_B300_Z100_w7_LCDM/snapshot_000.0"
    snapshot = open(filename, 'r')
    #h = header(snapshot)
    #print (h.n_part)
    #print (h.mass)
    #print (h.time)
    #print (h.redshift)
    #print (h.flag_sfr)
    #print (h.flag_feedback)
    #print (h.n_part_total)
    #print (h.n_part_TOTAL)
    #print (h.flag_cooling)
    #print (h.num_files)
    #print (h.box_size)
    #print (h.omega0)
    #print (h.omega_lambda)
    #print (h.hubble_param)
    #print (h.n_part_total_HW)
    #print (h.fill)
    p_list = read_data(snapshot)
    snapshot.close()
    # From here on, the data for all the particles is accessible
    # As an example, showing the positions of all the particles
    i=p_list[0]
    print (i.vel[0], i.vel[1], i.vel[2])
    print (i.pos[0], i.pos[1], i.pos[2])
    print (i.ID)
    print (i.rho)
    print (i.smoothing)

    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
if __name__ == '__main__':
    main()

