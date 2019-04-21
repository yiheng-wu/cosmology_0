#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



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
	//float Vel[3];
	//int ID;
}	*P;


int main(int argc, char **argv)
{
	char path[200],filename[200],wholename[200];
	char pathl[200],filenamel[200],wholenamel[200];


	sprintf(path,"/data/s6/data_Jing/6411/simu/");
	sprintf(filename,"pos6411.0600");

	sprintf(pathl,"/data/s1/simu/Jing6620/");
	sprintf(filenamel,"pos6620.0247.01");

	sprintf(wholename,"%s%s",path,filename);
	sprintf(wholenamel,"%s%s",pathl,filenamel);

//	load_jing(wholename);
	load_ljing(wholenamel);

	profile();

/*
	printf("%lf %lf %lf\n",P[178956969].Pos[0],P[178956969].Pos[1],P[178956969].Pos[2]);
	int i=0;
	for (i=0;i<10;i++)
	{
		printf("%lf %lf %lf\n",P[178956969+i].Pos[0],P[178956969+i].Pos[1],P[178956969+i].Pos[2]);
		printf("%lf %lf %lf\n",P[178956969*2+i].Pos[0],P[178956969*2+i].Pos[1],P[178956969*2+i].Pos[2]);
		printf("%lf %lf %lf\n",P[178956969*3+i].Pos[0],P[178956969*3+i].Pos[1],P[178956969*3+i].Pos[2]);
	}
*/
	
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
	fprintf(stderr,"SKIP=%d\n",dummy);
	fread(&headerl, sizeof(headerl),1,fd);
	SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
	printf("%ld\n",headerl.np);
	printf("%ld\n",headerl.ips);
	printf("%f\n",headerl.ztp);
	printf("%f\n",headerl.omegat);
	printf("%f\n",headerl.lamdat);
	printf("%f\n",headerl.boxsize);
	allocate_memory(headerl.np/4);

	SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
	int64_t np;
	float tem;
	int64_t N_b=0;//读取的字节数
	int64_t N0=2147483639;//fortran单块字节数
	int xyz;

/*
	while(!feof(fd))//判断文件末尾
	{
		fseek(fd,abs(dummy),1);
		SKIP;
		fprintf(stderr,"SKIP=%d\n",dummy);
		if(feof(fd)) fprintf(stderr,"END");
		SKIP;
		fprintf(stderr,"SKIP=%d\n",dummy);
		N_b++;
		if(feof(fd)) fprintf(stderr,"END");
	}
	
		fprintf(stderr,"N=%d\n",N_b);
*/

	for(np=0;np<headerl.np/4;np++)
	{
		for(xyz=0;xyz<3;xyz++)
		{
			fread(&tem,sizeof(tem),1,fd);
			P[np].Pos[xyz]=tem;
			N_b+=4;
			if(N_b==N0-3)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				fread(&joint.c[1],1,1,fd);
				fread(&joint.c[2],1,1,fd);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				fread(&joint.c[3],1,1,fd);
				P[np].Pos[2]=joint.f;
				xyz++;
				N_b+=4;
			}
			else if(N_b==N0*2-2)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				fread(&joint.c[1],1,1,fd);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				fread(&joint.c[2],1,1,fd);
				fread(&joint.c[3],1,1,fd);
				P[np].Pos[2]=joint.f;
				xyz++;
				N_b+=4;
			}
			else if(N_b==N0*3-1)
			{
				float_char joint;
				fread(&joint.c[0],1,1,fd);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				fread(&joint.c[1],1,1,fd);
				fread(&joint.c[2],1,1,fd);
				fread(&joint.c[3],1,1,fd);
				P[np].Pos[2]=joint.f;
				xyz++;
				N_b+=4;
			}
			else if(N_b==N0*4)
			{
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
				SKIP;
	fprintf(stderr,"SKIP=%d\n",dummy);
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
int profile(void)
{
	int pro_data[1200][1200];
	int i,j,k;
	for(i=0;i<1200;i++)
	{
		for(j=0;j<1200;j++)
		{
			pro_data[i][j]=0;
		}
	}
	for(i=0;i<header0.np;i++)
	{
		if (P[i].Pos[0]>0.5 && P[i].Pos[0]<0.501)
		{
			pro_data[(int)(P[i].Pos[1]*1200)][(int)(P[i].Pos[2]*1200)]+=1;
		}
	}
	FILE *new=fopen("SLICE","w");
	int64_t All=0;
	for(i=0;i<1200;i++)
	{
		for(j=0;j<1200;j++)
		{
			All+=pro_data[i][j];
			fprintf(new,"%d\n",pro_data[i][j]);
		}
	}
	fclose(new);
	printf("%ld\n",All);
}
