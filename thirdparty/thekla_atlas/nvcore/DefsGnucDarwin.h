#ifndef NV_CORE_H
#error "Do not include this file directly."
#endif

#include <stdint.h> // uint8_t, int8_t, ... uintptr_t
#include <stddef.h> // operator new, size_t, NULL

// Function linkage
#define DLL_IMPORT
#if __GNUC__ >= 4
#	define DLL_EXPORT __attribute__((visibility("default")))
#	define DLL_EXPORT_CLASS DLL_EXPORT
#else
#	define DLL_EXPORT
#	define DLL_EXPORT_CLASS
#endif

// Function calling modes
#if NV_CPU_X86
#	define NV_CDECL 	__attribute__((cdecl))
#	define NV_STDCALL	__attribute__((stdcall))
#else
#	define NV_CDECL 
#	define NV_STDCALL
#endif

#define NV_FASTCALL		__attribute__((fastcall))
#define NV_FORCEINLINE	__attribute__((always_inline)) inline
#define NV_DEPRECATED   __attribute__((deprecated))
#if NV_OS_IOS
#define NV_THREAD_LOCAL // @@ IC: Looks like iOS does not have support for TLS declarations.
#else
#define NV_THREAD_LOCAL __thread
#endif

#if __GNUC__ > 2
#define NV_PURE     __attribute__((pure))
#define NV_CONST    __attribute__((const))
#else
#define NV_PURE
#define NV_CONST
#endif

#define NV_NOINLINE __attribute__((noinline))

// Define __FUNC__ properly.
#if __STDC_VERSION__ < 199901L
#	if __GNUC__ >= 2
#		define __FUNC__ __PRETTY_FUNCTION__	// __FUNCTION__
#	else
#		define __FUNC__ "<unknown>"
#	endif
#else
#	define __FUNC__ __PRETTY_FUNCTION__
#endif

#define restrict    __restrict__
