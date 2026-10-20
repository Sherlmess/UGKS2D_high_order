#include <stdint.h>
#include <unistd.h>
#include <sys/resource.h>
#include <stdio.h>
#ifdef __APPLE__
#include <mach/mach.h>
#endif
int64_t mpi2_peak_rss(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1;
#ifdef __APPLE__
    return (int64_t)usage.ru_maxrss;
#else
    return (int64_t)usage.ru_maxrss * 1024;
#endif
}
int64_t mpi2_current_rss(void) {
#ifdef __APPLE__
    struct mach_task_basic_info info;
    mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
    if(task_info(mach_task_self(), MACH_TASK_BASIC_INFO,(task_info_t)&info,&count)!=KERN_SUCCESS)return -1;
    return (int64_t)info.resident_size;
#else
    long pages=0, resident=0;
    FILE *f=fopen("/proc/self/statm","r");
    if(!f)return -1;
    int valid=fscanf(f,"%ld %ld",&pages,&resident);fclose(f);
    if(valid!=2)return -1;
    return (int64_t)resident*sysconf(_SC_PAGESIZE);
#endif
}
/* Generation checksum: unsigned arithmetic has defined wraparound. Zero means I/O failure. */
int64_t s2_file_hash(const char *path) {
    FILE *f=fopen(path,"rb");
    if(!f)return 0;
    uint64_t h=UINT64_C(14695981039346656037);
    unsigned char buffer[16384];size_t n;
    while((n=fread(buffer,1,sizeof buffer,f))>0) {
        for(size_t i=0;i<n;i++) {h^=buffer[i];h*=UINT64_C(1099511628211);}
    }
    int failed=ferror(f);if(fclose(f)!=0)failed=1;
    return failed?0:(int64_t)h;
}
