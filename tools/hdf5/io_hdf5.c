#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>
#include "unistd.h"
#include "/data/s5/yhwu/soft/rockstar/universal_constants.h"
#include "/data/s5/yhwu/soft/rockstar/check_syscalls.c"
#include "/data/s5/yhwu/soft/rockstar/config_vars.h"
#include "/data/s5/yhwu/soft/rockstar/config.h"
#include "/data/s5/yhwu/soft/rockstar/particle.h"

hid_t check_H5Fopen(char *filename) {
  hid_t HDF_FileID = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (HDF_FileID < 0) {
    fprintf(stderr, "[Error] Failed to open HDF5 file %s!\n", filename);
    exit(1);
  }
  return HDF_FileID;
}

hid_t check_H5Gopen(hid_t HDF_FileID, char *gid, char *filename) {
  hid_t HDF_GroupID = H5Gopen(HDF_FileID, gid,H5P_DEFAULT);
  if (HDF_GroupID < 0) {
    fprintf(stderr, "[Error] Failed to open group %s in HDF5 file %s!\n", gid, filename);
    exit(1);
  }
  return HDF_GroupID;
}

hid_t check_H5Dopen(hid_t HDF_GroupID, char *dataid, char *gid, char *filename){
  hid_t HDF_DatasetID = H5Dopen(HDF_GroupID, dataid,H5P_DEFAULT);
  if (HDF_DatasetID < 0) {
    fprintf(stderr, "[Error] Failed to open dataset %s/%s in HDF5 file %s!\n", gid, dataid, filename);
    exit(1);
  }
  return HDF_DatasetID;
}

hid_t check_H5Dget_space(hid_t HDF_DatasetID) {
  hid_t HDF_DataspaceID = H5Dget_space(HDF_DatasetID);
  if (HDF_DataspaceID < 0) {
    fprintf(stderr, "[Error] Failed to get HDF5 dataspace!\n");
    exit(1);
  }
  return HDF_DataspaceID;
}

void check_H5Dread(hid_t HDF_DatasetID, hid_t type, void *buffer, char *dataid, char *gid, char *filename) {
  if (H5Dread(HDF_DatasetID, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0){
    fprintf(stderr, "[Error] failed to read dataspace %s/%s in HDF5 file %s\n", gid, dataid, filename);
    exit(1);
  }
}

hid_t check_H5Aopen_name(hid_t HDF_GroupID, char *dataid, char *gid, char *filename){
  hid_t HDF_AttrID = H5Aopen_name(HDF_GroupID, dataid);
  if (HDF_AttrID < 0) {
    fprintf(stderr, "[Error] Failed to open attribute %s/%s in HDF5 file %s!\n", gid, dataid, filename);
    exit(1);
  }
  return HDF_AttrID;
}

hid_t check_H5Aget_space(hid_t HDF_AttrID) {
  hid_t HDF_DataspaceID = H5Aget_space(HDF_AttrID);
  if (HDF_AttrID < 0) {
    fprintf(stderr, "[Error] Failed to get HDF5 dataspace!\n");
    exit(1);
  }
  return HDF_DataspaceID;
}

void check_H5Aread(hid_t HDF_AttrID, hid_t type, void *buffer, char *dataid, char *gid, char *filename) {
  if (H5Aread(HDF_AttrID, type, buffer) < 0) {
    fprintf(stderr, "[Error] failed to read attribute %s/%s in HDF5 file %s\n", gid, dataid, filename);
    exit(1);
  }
}

void check_H5Sselect_all(hid_t HDF_DataspaceID) {
  if (H5Sselect_all(HDF_DataspaceID)<0 || H5Sselect_valid(HDF_DataspaceID)<=0) {
    fprintf(stderr, "[Error] Failed to select all elements in HDF5 dataspace!\n");
    exit(1);
  }
}

int64_t check_H5Sget_simple_extent_ndims(hid_t HDF_DataspaceID) {
  int64_t ndims = H5Sget_simple_extent_ndims(HDF_DataspaceID);
  if (ndims < 0) {
    fprintf(stderr, "[Error] Failed to get number of dimensions for HDF5 dataspace!\n");
    exit(1);
  }
  return ndims;
}

void check_H5Sget_simple_extent_dims(hid_t HDF_DataspaceID, hsize_t *dimsize) {
  if (H5Sget_simple_extent_dims(HDF_DataspaceID, dimsize, NULL)<0) {
    fprintf(stderr, "[Error] Failed to get dimensions for HDF5 dataspace!\n");
    exit(1);
  }
}

float tng_readheader_float(hid_t groupid, char *filename, char *objname)
{
	char *group="Header";
	hid_t attrid=check_H5Aopen_name(groupid,objname,group,filename);
	hid_t dataspaceid=check_H5Aget_space(attrid);

	check_H5Sselect_all(dataspaceid);

	float data=0.0;
	check_H5Aread(attrid,H5T_NATIVE_FLOAT,&data,objname,group,filename);

	H5Sclose(dataspaceid);
	H5Aclose(attrid);
}
void tng_readheader_array(hid_t groupid, char *filename,char *objname,hid_t type, void *data)
{
	char *group="Header";
	hid_t attrid=check_H5Aopen_name(groupid,objname,group,filename);
	hid_t dataspaceid=check_H5Aget_space(attrid);
	check_H5Sselect_all(dataspaceid);

	int ndims=check_H5Sget_simple_extent_ndims(dataspaceid);
	hsize_t dimsize=0;
	check_H5Sget_simple_extent_dims(dataspaceid,&dimsize);

	check_H5Aread(attrid,type,data,objname,group,filename);

	H5Sclose(dataspaceid);
	H5Aclose(attrid);
}

int main()
{
	
	char *dm="PartType1";
	char *loc="PartType1/Coordinates";
	char *vel="PartType1/Velocities";
	char *id="PartType1/ParticleIDs";
	char *dir="/data/s5/yhwu/data/TNG/snap/snap_099.0.hdf5";
	hid_t fid=check_H5Fopen(dir);
	hid_t dmid=check_H5Gopen(fid,dm,dir);

	hid_t headid=check_H5Gopen(fid,"Header",dir);
	float boxsize=tng_readheader_float(headid,dir,"MassTable");

	int64_t NumPart[6];
	tng_readheader_array(headid,dir,"NumPart_Total",H5T_NATIVE_UINT32,NumPart);

	printf("BoxSize=%d\n",NumPart[0]);
	printf("BoxSize=%d\n",NumPart[1]);
	printf("BoxSize=%d\n",NumPart[2]);
	printf("BoxSize=%d\n",NumPart[3]);


	return 0;


























	herr_t status = H5Lget_info(fid,dm, NULL, H5P_DEFAULT);
	if (status == 0)
	{
		//Position
		struct location{ double pos[3];}; 
		struct location *ploc;

		hid_t locid=check_H5Dopen(fid,loc,dm,dir);
		hid_t velid=check_H5Dopen(fid,vel,dm,dir);

		hid_t locspace=check_H5Dget_space(locid);
		hid_t velspace=check_H5Dget_space(velid);
		int ndims=check_H5Sget_simple_extent_ndims(locspace);
		check_H5Sselect_all(locspace);
		hsize_t *num_p;
		status=H5Sget_simple_extent_dims(locspace,num_p,NULL);
		ploc=(struct location *)malloc((*num_p)*sizeof(struct location));
		check_H5Dread(locid, H5T_NATIVE_DOUBLE,ploc,loc,dm,dir);
		H5Dclose(locid);
		printf("%d\n",*num_p);

		printf("%f\n",ploc[*num_p-1].pos[1]);
/*		FILE *fp=fopen("pos.txt","w");
		int i;
		for (i=0;i<*num_p;i+)
		{
			fprintf(fp,"%f %f %f\n",ploc[i].pos[0],ploc[i].pos[1],ploc[i].pos[2]);
		}
*/		
		free(ploc);
		//Velocity
		struct velocity{ double vel[3];}; 
		struct velocity *pvel;

/*

		pvel=(struct velocity *)malloc((*num_p)*sizeof(struct velocity));
		check_H5Dread(locid, H5T_NATIVE_DOUBLE,pvel,vel,dm,dir);
		H5Dclose(velid);

		printf("%f\n",pvel[*num_p].vel[0]);
		int i;
		FILE *fp=fopen("vel.txt","w");
		for (i=0;i<*num_p;i++)
		{
			fprintf(fp,"%f %f %f\n",pvel[i].vel[0],pvel[i].vel[1],pvel[i].vel[2]);
		}
*/

//		hid_t velid=check_H5Dopen(fid,vel,dm,dir);
//		hid_t idid=check_H5Dopen(fid,id,dm,dir);

//		hid_t velspace=check_H5Dget_space(velid);
//		hid_t idspace=check_H5Dget_space(idid);




	//	*ploc=(struct location *)check_realloc(*ploc,(10)*sizeof(struct location),"Allocating particle location.");

/*
		float *vel[3];

		*loc=


		int j;
		for(j=0;j<*num_p;j++)
		{
			//p[j].pos[1]=
		}
*/		



//		H5Dclose(velid);
//		H5Dclose(idid);
	}
	else
		printf("The group either does NOT exist or some other error occurred.\n");
	H5Fclose(fid);
	return 0;
}

