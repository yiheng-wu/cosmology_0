#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <inttypes.h>

int main()
{
	char *dir="/data/s5/yhwu/data/hdf5/fof_subhalo_tab_000.0.hdf5";
	hid_t hdf=H5Fopen(dir);



	return 0;
}
