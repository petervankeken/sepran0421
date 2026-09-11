#include <sys/param.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

csecond_(cptime)
float *cptime;
{ 
  struct tms buffer;
  int utime,stime;

  times(&buffer);
  utime = buffer.tms_utime;
  stime = buffer.tms_stime;
  *cptime = utime*1.0/HZ;
}
csecond (cptime)
float *cptime;
{ 
  struct tms buffer;
  int utime,stime;

  times(&buffer);
  utime = buffer.tms_utime;
  stime = buffer.tms_stime;
  *cptime = utime*1.0/HZ;
}

/*second_(cptime)
float *cptime;
{ 
  struct tms buffer;
  int utime,stime;

  times(&buffer);
  utime = buffer.tms_utime;
  stime = buffer.tms_stime;
  *cptime = utime*1.0/HZ;
}
second (cptime)
float *cptime;
{ 
  struct tms buffer;
  int utime,stime;

  times(&buffer);
  utime = buffer.tms_utime;
  stime = buffer.tms_stime;
  *cptime = utime*1.0/HZ;
}
*/
