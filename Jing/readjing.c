#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGRID 3000
int pro_data[NGRID][NGRID];

//int大小端转换
unsigned int UISWAP(int x)
{
	return \
	(x&0xff000000)>>24|\
	(x&0x00ff0000)>>8|\
	(x&0x0000ff00)<<8|\
	(x&0x000000ff)<<24;
}
//64位int大小端转换
int64_t LISWAP(int64_t x)
{
	int l,h;
	int64_t t;
	l=x&0xFFFFFFFF;
	h=(x>>32)&0xFFFFFFFF;

	t=UISWAP(l);

	return (t<<32)|UISWAP(h);
}
//float大小端转换
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


struct header_jing
{
	int np;
	int ips;
	float ztp;
	float omegat;
	float lamdat;
	float boxsize;
	float xscale;
	float vscale;
}	header0;

struct particle_data
{
	float Pos[3];
}	*P;


int main(int argc, char **argv)
{
	char *f="/data/s6/data_Jing/6411/simu/pos6411.0600";
	initial_profile();
	load_jing(f);

	profile();
	
	return 0;
}

int load_jing(char *fname)
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
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
	fread(&header0, sizeof(header0),1,fd);
	SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
	//Jing类型需要做大小端转换，LJing不需要转换
	header0.np=UISWAP(header0.np);
	header0.ips=UISWAP(header0.ips);
	header0.ztp=FSWAP(header0.ztp);
	header0.omegat=FSWAP(header0.omegat);
	header0.lamdat=FSWAP(header0.lamdat);
	header0.boxsize=FSWAP(header0.boxsize);
	header0.xscale=FSWAP(header0.xscale);
	header0.vscale=FSWAP(header0.vscale);
	printf("%ld\n",header0.np);
	printf("%ld\n",header0.ips);
	printf("%f\n",header0.ztp);
	printf("%f\n",header0.omegat);
	printf("%f\n",header0.lamdat);
	printf("%f\n",header0.boxsize);
	allocate_memory(header0.np);

	SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
	int np;
	float tem;
	int64_t N_b=0;//读取的字节数
	int64_t N0=2147483639;//fortran单块字节数
	int xyz;

	for(xyz=0;xyz<3;xyz++)
	{
		for(np=0;np<header0.np;np++)
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
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				fread(&joint.c[3],1,1,fd);
				P[np+1].Pos[xyz]=FSWAP(joint.f);
				np++;
				N_b+=4;
			}
			else if(N_b==N0*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				fread(&joint.c[1],1,1,fd);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				fread(&joint.c[2],1,1,fd);
				fread(&joint.c[3],1,1,fd);
				P[np+1].Pos[xyz]=FSWAP(joint.f);
				np++;
				N_b+=4;
			}
			else if(N_b==N0*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				fread(&joint.c[1],1,1,fd);
				fread(&joint.c[2],1,1,fd);
				fread(&joint.c[3],1,1,fd);
				P[np+1].Pos[xyz]=FSWAP(joint.f);
				np++;
				N_b+=4;
			}
			else if(N_b==N0*4)
			{
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
				SKIP;
	fprintf(stderr,"SKIP=%d\n",UISWAP(dummy));
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
	int i,j;
	for(i=0;i<header0.np;i++)
	{
		if (P[i].Pos[0]>0.5 && P[i].Pos[0]<0.5005)
		{
			pro_data[(int)(P[i].Pos[1]*NGRID)][(int)(P[i].Pos[2]*NGRID)]+=1;
		}
	}
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
