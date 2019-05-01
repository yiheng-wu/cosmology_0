#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>

#define BOXSIZE 1200
#define NP 28991029248
#define NGRID 1200
#define timed_output(...) {print_time(); fprintf(stderr, __VA_ARGS__);}
void print_time(void)
{
	int64_t time_now=time(NULL);
	fprintf(stderr,"[%6"PRId64"s]",time_now);
}

int pro_data[NGRID][NGRID];

typedef union FLOAT_CHAR
{
	float f;
	char c[4];
}float_char;
float FSWAP(float x)
{
	float_char f1,f2;
	f1.f=x;
	f2.c[0]=f1.c[3];
	f2.c[1]=f1.c[2];
	f2.c[2]=f1.c[1];
	f2.c[3]=f1.c[0];
	return f2.f;
}

struct header_ljing
{
	int64_t np;
	int64_t ips;
	float ztp;
	float omegat;
	float lamdat;
	float boxsize;
	float xscale;
	float vscale;
}	headerl;

struct particle_data
{
	float Pos[3];
}	*P;


int main(int argc, char **argv)
{
	char *f1="/data/s1/simu/Jing6620/pos6620.5000.01";
	char *f2="/data/s1/simu/Jing6620/pos6620.5000.02";
	char *f3="/data/s1/simu/Jing6620/pos6620.5000.03";
	char *f4="/data/s1/simu/Jing6620/pos6620.5000.04";

	allocate_memory(NP/4);
	initial_profile();
	load_ljing(f1);
	profile();
	load_ljing(f2);
	profile();
	load_ljing(f3);
	profile();
	load_ljing(f4);
	profile();
	slice();	
	
	return 0;
}

int load_ljing(char *fname)
{
	FILE *fd;
	unsigned int dummy;

#define SKIP fread(&dummy, sizeof(dummy),1,fd);

	if(!(fd=fopen(fname,"r")))
	{
		printf("can't open file '%s'\n",fname);
		exit(0);
	}

	SKIP;
	if (dummy==40)
	{
		fread(&headerl, sizeof(headerl),1,fd);
		SKIP;
		SKIP;
	}
	else 
	{
		headerl.boxsize=BOXSIZE;
		headerl.np=NP;
	}

	int64_t np;
	float tem;
	int64_t N_b=0;//读取的字节数
	int64_t N0=2147483639;//fortran单块字节数
	int xyz;

	for(np=0;np<headerl.np/4;np++)
	{
		for(xyz=0;xyz<3;xyz++)
		{
			fread(&tem,sizeof(tem),1,fd);
			P[np].Pos[xyz]=FSWAP(tem);
			N_b+=4;
			if(N_b==N0-3)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				fread(&joint.c[1],1,1,fd);
				fread(&joint.c[2],1,1,fd);
				SKIP;
				SKIP;
				fread(&joint.c[3],1,1,fd);
				P[np].Pos[2]=FSWAP(joint.f);
				xyz++;
				N_b+=4;
			}
			else if(N_b==N0*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				fread(&joint.c[1],1,1,fd);
				SKIP;
				SKIP;
				fread(&joint.c[2],1,1,fd);
				fread(&joint.c[3],1,1,fd);
				P[np].Pos[2]=FSWAP(joint.f);
				xyz++;
				N_b+=4;
			}
			else if(N_b==N0*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				SKIP;
				SKIP;
				fread(&joint.c[1],1,1,fd);
				fread(&joint.c[2],1,1,fd);
				fread(&joint.c[3],1,1,fd);
				P[np].Pos[2]=FSWAP(joint.f);
				xyz++;
				N_b+=4;
			}
			else if(N_b==N0*4)
			{
				SKIP;
				SKIP;
				N_b=0;
			}
		}
	}
	fclose(fd);
}

int allocate_memory(int64_t np)
{
  printf("allocating memory...\n");

  if(!(P = malloc(np * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  //P--;				/* start with offset 1 */

  printf("allocating memory...done\n");
}
int initial_profile(void)
{
	int i,j;
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
		{
			pro_data[i][j]=0;
		}
	}
}
int profile(void)
{
	int64_t i;
	
	for(i=0;i<headerl.np/4;i++)
	{
		if (P[i].Pos[0]>0.5 && P[i].Pos[0]<0.5+0.5/headerl.boxsize)
		{
			pro_data[(int)(P[i].Pos[1]*NGRID)][(int)(P[i].Pos[2]*NGRID)]+=1;
		}
	}
}
int slice(void)
{
	int64_t i,j;
	FILE *new=fopen("/data/s5/yhwu/data/SLICE","w");
	int64_t All=0;
	for(i=0;i<NGRID;i++)
	{
		for(j=0;j<NGRID;j++)
		{
			All+=pro_data[i][j];
			fprintf(new,"%d\n",pro_data[i][j]);
		}
	}
	fclose(new);
	printf("%ld\n",All);
}
