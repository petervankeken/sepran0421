#include <stdio.h>
#include <sys/types.h>
#include <sys/file.h>
#include <math.h>
#include <string.h>

static  char name[80];

void rasout(NTOT,PIX,fname)
int *NTOT;
int PIX[];
char fname[];
{
  int n,isc,i;
  unsigned char fj;
  FILE *fp, *fopen();

  
/* Copy Fortran name to C character array to get rid of the spaces */
  for (i=0;i<80;i++){
      if (fname[i] != 32){ 
              name[i]=fname[i];
      }
      else {
         name[i]=0;
         break;
      }
  }
  fp = fopen(name,"w");
  n  = *NTOT;
  for (i=0;i<n;i++){
      fj = (unsigned char) PIX[i];
      putc(fj,fp);
  }
  fclose(fp);
}

void rasout_(NTOT,PIX,fname)
int *NTOT;
int PIX[];
char fname[];
{
  int n,isc,i;
  unsigned char fj;
  FILE *fp, *fopen();

  
/* Copy Fortran name to C character array to get rid of the spaces */
  for (i=0;i<80;i++){
      if (fname[i] != 32){ 
              name[i]=fname[i];
      }
      else {
         name[i]=0;
         break;
      }
  }
  fp = fopen(name,"w");
  n  = *NTOT;
  for (i=0;i<n;i++){
      fj = (unsigned char) PIX[i];
      putc(fj,fp);
  }
  fclose(fp);
}

RASOUT(NTOT,PIX,fname)
int *NTOT;
int PIX[];
char fname[];
{
  int n,isc,i;
  unsigned char fj;
  FILE *fp, *fopen();

  
/* Copy Fortran name to C character array to get rid of the spaces */
  for (i=0;i<80;i++){
      if (fname[i] != 32){ 
              name[i]=fname[i];
      }
      else {
         name[i]=0;
         break;
      }
  }
  fp = fopen(name,"w");
  n  = *NTOT;
  for (i=0;i<n;i++){
      fj = (unsigned char) PIX[i];
      putc(fj,fp);
  }
  fclose(fp);
}
