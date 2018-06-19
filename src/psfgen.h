
/* psfgen.h
 * Defines set of data structures used in creation of molecule structures
 * Exported here so that new modules can be written to interface with psfgen
*/ 

#ifndef PSFGEN_H
#define PSFGEN_H

#include "topo_defs.h"
#include "topo_mol.h"
#include "stringhash.h"

/* psfgen-specific data */
struct psfgen_data {
  int id, in_use, all_caps;
  topo_defs *defs;
  topo_mol *mol;
  stringhash *aliases;
  FILE* outstream;
};
typedef struct psfgen_data psfgen_data;

// Emit a compile error if int is not 32 bits
(void)sizeof(char[1 - 2*!!(sizeof(int) != 4)]);

// Some utility functions
char* strtoupper(const char *str, int all_caps) {
  char *s, *tmp;
  tmp = strdup(str);
  if ( all_caps ) {
    s=tmp;
    while ( *s ) { *s = toupper(*s); ++s; }
  }
  return tmp;
}

char* splitcolon(char *s) {
  if ( s ) {
    while ( *s && *s != ':' ) { ++s; }
    if ( *s ) *(s++) = 0; else s = 0;
  }
  return s;
}

#endif /* PSFGEN_H */
