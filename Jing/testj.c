#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>

#define timed_output(...) {print_time(); fprintf(stderr, __VA_ARGS__);}

int main()
{
	char *f="/data/s1/simu/Jing6611/pos6611.2448.01";
	int start,end;
	FILE *fd;

	fd=fopen(f,"r");
	while(feof(fd)==0)
	{
		fread(&start,sizeof(start),1,fd);
		fseek(fd,abs(start),1);
		fread(&end,sizeof(end),1,fd);
		printf("start=%d,end=%d\n",start,end);
	}
	return 0;
}
