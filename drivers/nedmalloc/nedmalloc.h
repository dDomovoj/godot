#ifdef NEDMALLOC_ENABLED

/* nedalloc, an alternative malloc implementation for multiple threads without
lock contention based on dlmalloc v2.8.3. (C) 2005-2009 Niall Douglas

Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/

#ifndef NEDMALLOC_H
#define NEDMALLOC_H

#include "typedefs.h"
#define MALLOC_ALIGNMENT DEFAULT_ALIGNMENT

#ifdef PSP_ENABLED
#define USE_LOCKS 0
#define HAVE_MMAP 0
#endif

/* See malloc.c.h for what each function does.

REPLACE_SYSTEM_ALLOCATOR on POSIX causes nedalloc's functions to be called
malloc, free etc. instead of nedmalloc, nedfree etc. You may or may not want
this. On Windows it causes nedmalloc to patch all loaded DLLs and binaries
to replace usage of the system allocator.

NO_NED_NAMESPACE prevents the functions from being defined in the nedalloc
namespace when in C++ (uses the global namespace instead).

NEDMALLOCEXTSPEC can be defined to be __declspec(dllexport) or
__attribute__ ((visibility("default"))) or whatever you like. It defaults
to extern unless NEDMALLOC_DLL_EXPORTS is set as it would be when building
nedmalloc.dll.

USE_LOCKS can be 2 if you want to define your own MLOCK_T, INITIAL_LOCK,
ACQUIRE_LOCK, RELEASE_LOCK, TRY_LOCK, IS_LOCKED and NULL_LOCK_INITIALIZER.

NEDMALLOC_DEBUG can be defined to cause DEBUG to be set differently for nedmalloc
than for the rest of the build. Remember to set NDEBUG to disable all assertion
checking too.

USE_MAGIC_HEADERS causes nedalloc to allocate an extra three sizeof(size_t)
to each block. nedpfree() and nedprealloc() can then automagically know when
to free a system allocated block. Enabling this typically adds 20-50% to
application memory usage.

ENABLE_TOLERANT_NEDMALLOC is automatically turned on if REPLACE_SYSTEM_ALLOCATOR
is set or the Windows DLL is being built. This causes nedmalloc to detect when a
system allocator block is passed to it and to handle it appropriately. Note that
without USE_MAGIC_HEADERS there is a very tiny chance that nedmalloc will segfault
on non-Windows builds (it uses Win32 SEH to trap segfaults on Windows and there
is no comparable system on POSIX).

USE_ALLOCATOR can be one of these settings (it defaults to 1):
  0: System allocator (nedmalloc now simply acts as a threadcache).
     WARNING: Intended for DEBUG USE ONLY - not all functions work correctly.
  1: dlmalloc

ENABLE_LARGE_PAGES enables support for requesting memory from the system in large
(typically >=2Mb) pages if the host OS supports this. These occupy just a single
TLB entry and can significantly improve performance in large working set applications.

ENABLE_FAST_HEAP_DETECTION enables special logic to detect blocks allocated
by the system heap. This avoids 1.5%-2% overhead when checking for non-nedmalloc
blocks, but it assumes that the NT and glibc heaps function in a very specific
fashion which may not hold true across OS upgrades.
*/

#include <stddef.h>   /* for size_t */

#ifndef NEDMALLOCEXTSPEC
 #ifdef NEDMALLOC_DLL_EXPORTS
  #ifdef WIN32
   #define NEDMALLOCEXTSPEC extern __declspec(dllexport)
  #elif defined(__GNUC__)
   #define NEDMALLOCEXTSPEC extern __attribute__ ((visibility("default")))
  #endif
  #ifndef ENABLE_TOLERANT_NEDMALLOC
   #define ENABLE_TOLERANT_NEDMALLOC 1
  #endif
 #else
  #define NEDMALLOCEXTSPEC extern
 #endif
#endif

#if __STDC_VERSION__ >= 199901L		/* C99 or better */
 #define RESTRICT restrict
#else
 #if defined(_MSC_VER) && _MSC_VER>=1400
  #define RESTRICT __restrict
 #endif
 #ifdef __GNUC__
  #define RESTRICT __restrict
 #endif
#endif
#ifndef RESTRICT
 #define RESTRICT
#endif

#if defined(_MSC_VER) && _MSC_VER>=1400
 #define NEDMALLOCPTRATTR __declspec(restrict)
 #define NEDMALLOCNOALIASATTR __declspec(noalias)
#endif
#ifdef __GNUC__
 #define NEDMALLOCPTRATTR __attribute__ ((malloc))
#endif
#ifndef NEDMALLOCPTRATTR
 #define NEDMALLOCPTRATTR
#endif
#ifndef NEDMALLOCNOALIASATTR
 #define NEDMALLOCNOALIASATTR
#endif

#ifndef USE_MAGIC_HEADERS
 #define USE_MAGIC_HEADERS 0
#endif

#ifndef USE_ALLOCATOR
 #define USE_ALLOCATOR 1 /* dlmalloc */
#endif

#if !USE_ALLOCATOR && !USE_MAGIC_HEADERS
#error If you are using the system allocator then you MUST use magic headers
#endif

#ifdef REPLACE_SYSTEM_ALLOCATOR
 #if USE_ALLOCATOR==0
  #error Cannot combine using the system allocator with replacing the system allocator
 #endif
 #ifndef ENABLE_TOLERANT_NEDMALLOC
  #define ENABLE_TOLERANT_NEDMALLOC 1
 #endif
 #ifndef WIN32	/* We have a dedicated patcher for Windows */
  #define nedmalloc               malloc
  #define nedcalloc               calloc
  #define nedrealloc              realloc
  #define nedfree                 free
  #define nedmemalign             memalign
  #define nedmallinfo             mallinfo
  #define nedmallopt              mallopt
  #define nedmalloc_trim          malloc_trim
  #define nedmalloc_stats         malloc_stats
  #define nedmalloc_footprint     malloc_footprint
  #define nedindependent_calloc   independent_calloc
  #define nedindependent_comalloc independent_comalloc
  #ifdef _MSC_VER
   #define nedblksize              _msize
  #endif
 #endif
#endif

#if defined(__cplusplus)
extern "C" {
#endif
struct nedmallinfo {
  size_t arena;    /* non-mmapped space allocated from system */
  size_t ordblks;  /* number of free chunks */
  size_t smblks;   /* always 0 */
  size_t hblks;    /* always 0 */
  size_t hblkhd;   /* space in mmapped regions */
  size_t usmblks;  /* maximum total allocated space */
  size_t fsmblks;  /* always 0 */
  size_t uordblks; /* total allocated space */
  size_t fordblks; /* total free space */
  size_t keepcost; /* releasable (via malloc_trim) space */
};
#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
 #if !defined(NO_NED_NAMESPACE)
namespace nedalloc {
 #else
extern "C" {
 #endif
 #define THROWSPEC throw()
#else
 #define THROWSPEC
#endif

/* These are the global functions */

/* Gets the usable size of an allocated block. Note this will always be bigger than what was
asked for due to rounding etc. Optionally returns 1 in isforeign if the block came from the
system allocator - note that there is a small (>0.01%) but real chance of segfault on non-Windows
systems when passing non-nedmalloc blocks if you don't use USE_MAGIC_HEADERS.
*/
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR size_t nedblksize(int *RESTRICT isforeign, void *RESTRICT mem) THROWSPEC;

NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR void nedsetvalue(void *v) THROWSPEC;

NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR NEDMALLOCPTRATTR void * nedmalloc(size_t size) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR NEDMALLOCPTRATTR void * nedcalloc(size_t no, size_t size) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR NEDMALLOCPTRATTR void * nedrealloc(void *mem, size_t size) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR void   nedfree(void *mem) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR NEDMALLOCPTRATTR void * nedmemalign(size_t alignment, size_t bytes) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR struct nedmallinfo nedmallinfo(void) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR int    nedmallopt(int parno, int value) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR void*  nedmalloc_internals(size_t *granularity, size_t *magic) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR int    nedmalloc_trim(size_t pad) THROWSPEC;
NEDMALLOCEXTSPEC void   nedmalloc_stats(void) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR size_t nedmalloc_footprint(void) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR NEDMALLOCPTRATTR void **nedindependent_calloc(size_t elemsno, size_t elemsize, void **chunks) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCNOALIASATTR NEDMALLOCPTRATTR void **nedindependent_comalloc(size_t elems, size_t *sizes, void **chunks) THROWSPEC;

/* Destroys the system memory pool used by the functions above.
Useful for when you have nedmalloc in a DLL you're about to unload.
If you call ANY nedmalloc functions after calling this you will
get a fatal exception!
*/
NEDMALLOCEXTSPEC void neddestroysyspool() THROWSPEC;

/* These are the pool functions */
struct nedpool_t;
typedef struct nedpool_t nedpool;

/* Creates a memory pool for use with the nedp* functions below.
Capacity is how much to allocate immediately (if you know you'll be allocating a lot
of memory very soon) which you can leave at zero. Threads specifies how many threads
will *normally* be accessing the pool concurrently. Setting this to zero means it
extends on demand, but be careful of this as it can rapidly consume system resources
where bursts of concurrent threads use a pool at once.
*/
NEDMALLOCEXTSPEC NEDMALLOCPTRATTR nedpool *nedcreatepool(size_t capacity, int threads) THROWSPEC;

/* Destroys a memory pool previously created by nedcreatepool().
*/
NEDMALLOCEXTSPEC void neddestroypool(nedpool *p) THROWSPEC;

/* Returns a zero terminated snapshot of threadpools existing at the time of call. Call
nedfree() on the returned list when you are done. Returns zero if there is only the
system pool in existence.
*/
NEDMALLOCEXTSPEC nedpool **nedpoollist() THROWSPEC;

/* Sets a value to be associated with a pool. You can retrieve this value by passing
any memory block allocated from that pool.
*/
NEDMALLOCEXTSPEC void nedpsetvalue(nedpool *p, void *v) THROWSPEC;

/* Gets a previously set value using nedpsetvalue() or zero if memory is unknown.
Optionally can also retrieve pool. You can detect an unknown block by the return
being zero and *p being unmodifed.
*/
NEDMALLOCEXTSPEC void *nedgetvalue(nedpool **p, void *mem) THROWSPEC;

/* Trims the thread cache for the calling thread, returning any existing cache
data to the central pool. Remember to ALWAYS call with zero if you used the
system pool. Setting disable to non-zero replicates neddisablethreadcache().
*/
NEDMALLOCEXTSPEC void nedtrimthreadcache(nedpool *p, int disable) THROWSPEC;

/* Disables the thread cache for the calling thread, returning any existing cache
data to the central pool. Remember to ALWAYS call with zero if you used the
system pool.
*/
NEDMALLOCEXTSPEC void neddisablethreadcache(nedpool *p) THROWSPEC;


NEDMALLOCEXTSPEC NEDMALLOCPTRATTR void * nedpmalloc(nedpool *p, size_t size) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCPTRATTR void * nedpcalloc(nedpool *p, size_t no, size_t size) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCPTRATTR void * nedprealloc(nedpool *p, void *mem, size_t size) THROWSPEC;
NEDMALLOCEXTSPEC void   nedpfree(nedpool *p, void *mem) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCPTRATTR void * nedpmemalign(nedpool *p, size_t alignment, size_t bytes) THROWSPEC;
NEDMALLOCEXTSPEC struct nedmallinfo nedpmallinfo(nedpool *p) THROWSPEC;
NEDMALLOCEXTSPEC int    nedpmallopt(nedpool *p, int parno, int value) THROWSPEC;
NEDMALLOCEXTSPEC int    nedpmalloc_trim(nedpool *p, size_t pad) THROWSPEC;
NEDMALLOCEXTSPEC void   nedpmalloc_stats(nedpool *p) THROWSPEC;
NEDMALLOCEXTSPEC size_t nedpmalloc_footprint(nedpool *p) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCPTRATTR void **nedpindependent_calloc(nedpool *p, size_t elemsno, size_t elemsize, void **chunks) THROWSPEC;
NEDMALLOCEXTSPEC NEDMALLOCPTRATTR void **nedpindependent_comalloc(nedpool *p, size_t elems, size_t *sizes, void **chunks) THROWSPEC;

#if defined(__cplusplus)
}
#endif

#endif

#endif
