#include <stdio.h>
#include <stdlib.h>

#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
// unistd.h for sysconf
#include <unistd.h>
#elif defined(__APPLE__) && defined(__MACH__)
// sys/resource.h for getrusage
#include <sys/resource.h>
//#include <sys/time.h>
#endif


long get_memory_usage(){
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    long rss = 0L;  // resident set size ( in memory page )
    FILE* fmem = NULL;
    if ( (fmem = fopen( "/proc/self/statm", "r" )) == NULL )
        return -1;  // file read error
    if ( fscanf( fmem, "%*ld %ld %*ld %*ld %*ld %*ld %*ld", &rss ) != 1 )
    {
        fclose( fmem );
        return -1;  // fscanf read error
    }
    fclose( fmem );
    return rss * sysconf( _SC_PAGESIZE);    // in bytes

#elif defined(__APPLE__) && defined(__MACH__)
    struct rusage ruse;
    if(0 == getrusage(RUSAGE_SELF, &ruse))
        return ruse.ru_maxrss; // in bytes
    else
        return -1;  // getrusage error
#else
    return 0;  //unsupported
#endif
}

void get_memory_usage_str(char * mem_usage_str){
    // output format:
    //       ${original string} xxx.xxx GiB
    long mem_usage_bytes;

    mem_usage_bytes = get_memory_usage();
    if ( mem_usage_bytes > 1024*1024*1024 ){
        sprintf(mem_usage_str, "%s%.3f GiB", mem_usage_str, ((double) mem_usage_bytes) /1024/1024/1024 );
    }
    else if ( mem_usage_bytes > 1024*1024 ){
        sprintf(mem_usage_str, "%s%.3f MiB", mem_usage_str, ((double) mem_usage_bytes) /1024/1024 );
    }
    else {
        sprintf(mem_usage_str, "%s%.3f KiB", mem_usage_str, ((double) mem_usage_bytes) /1024 );
    }
    return;
}
