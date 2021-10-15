#ifndef USEFULIO_H
#define USEFULIO_H

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

long get_memory_usage();
void get_memory_usage_str(char * mem_usage_str);

#endif
