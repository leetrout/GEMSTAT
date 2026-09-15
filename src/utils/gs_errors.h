#ifndef GS_ERRORS_H
#define GS_ERRORS_H

/**
 * Inspired by : http://stackoverflow.com/a/3767904
 *
 */

#include <cstdio>
#include <cstdlib>

/**
 * Check a condition that depends on user input or on program state that a
 * release build must still verify.  Unlike assert(), this is never compiled
 * out: with NDEBUG the old definition became a no-op and the program went on
 * to crash (or worse, run) on the invalid state it was meant to reject.
 */
#define ASSERT_MESSAGE(condition, message)\
	do { if( !( condition ) ) { fprintf(stderr, "%s\n", message); fprintf(stderr, "Failed check: %s (%s:%d)\n", #condition, __FILE__, __LINE__); exit(1); } } while(0)

#endif
