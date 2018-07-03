#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "python_psfgen.h"
#include "psfgen.h"
#include "hasharray.h"
#include "pdb_file_extract.h"
#include "psf_file_extract.h"
#include "topo_defs.h"
#include "topo_mol_struct.h"
#include "topo_mol_output.h"
#include "charmm_parse_topo_defs.h"
#include "stringhash.h"
#include "extract_alias.h"

/* Helper functions */
void python_msg(void *v, const char *msg)
{
    v = (FILE*)v; // print to this file
    fprintf(v, msg);
    fprintf(v, "\n");
}

/* Wrapper function for getting char* from a python string */
static char* as_charptr(PyObject *target)
{
#if PY_MAJOR_VERSION >=3
    return PyUnicode_AsUTF8(target);
#else
    return PyString_AsString(target);
#endif
}

/* Wrapper function for turning char* into a python string */
static PyObject* as_pystring(char *target)
{
    PyObject *result;
    if (!target) {
        PyErr_SetString(PyExc_ValueError, "cannot convert null string");
        return NULL;
    }
#if PY_MAJOR_VERSION >= 3
    result = PyUnicode_FromString(target);
#else
    result = PyString_FromString(target);
#endif

    if (!result || PyErr_Occurred()) {
        PyErr_Format(PyExc_ValueError, "cannot convert char* '%s'", target);
        return NULL;
    }
    return result;
}

/* Similar wrapper function for turning int into a python int/long */
static PyObject* as_pyint(int target)
{
    PyObject *result;
#if PY_MAJOR_VERSION >= 3
    result = PyLong_FromLong((long) target);
#else
    result = PyInt_FromLong((long) target);
#endif
    if (!result || PyErr_Occurred()) {
        PyErr_Format(PyExc_ValueError, "cannot convert int %d", target);
        return NULL;
    }
    return result;
}

/* Initialization / destruction functions */
static PyObject* py_init_mol(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "outfd", NULL};
    PyObject *capsule;
    psfgen_data *data;
    int outfd = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|i:__init__", kwnames,
                                     &outfd)) {
        return NULL;
    }

    data = malloc(sizeof(psfgen_data));

    // Initialize topologies
    data->defs = topo_defs_create();

    // Initialize aliases
    data->aliases = stringhash_create();
    data->mol = topo_mol_create(data->defs);

    // Initialize other stuffs
    data->id = 0; // Doesn't matter since data is per class instance
    data->in_use = 0;
    data->all_caps = 1;

    /*
     * Handle output file argument.. Default to stdout
     * This doesn't need to be closed in del_mol because the file
     * descriptor is closed in the Python destructor if it is used
     */
    data->outstream = outfd ? fdopen(outfd, "a") : stdout;
    topo_defs_error_handler(data->defs, data->outstream, python_msg);
    topo_mol_error_handler(data->mol, data->outstream, python_msg);

    // Encapsulate psf_data object and return it
    capsule = PyCapsule_New(data, NULL, NULL);
    if (!capsule || PyErr_Occurred())
        return NULL;
    return capsule;
}

static PyObject* py_del_mol(PyObject *self, PyObject *stateptr)
{
    psfgen_data *data;

    // Unpack molecule capsule
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
       return NULL;

    // Invoke cleanup functions
    topo_mol_destroy(data->mol);
    topo_defs_destroy(data->defs);
    stringhash_destroy(data->aliases);
    free(data);

    Py_INCREF(Py_None);
    return Py_None;
}

/* Aliases and names and stuff */
static PyObject* py_alias(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "type", (char*) "name",
                       (char*) "newname", (char*) "resname", NULL};
    char *type, *name, *newname, *resname = NULL;
    PyObject *stateptr;
    psfgen_data *data;

    // Parse arguments and put them in the right format
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Osss|s:alias", kwnames,
                                     &stateptr, &type, &name, &newname,
                                     &resname)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    name = strtoupper(name, data->all_caps);
    newname = strtoupper(newname, data->all_caps);

    if (!strcasecmp(type, "residue")) {
        fprintf(data->outstream, "Aliasing residue %s to %s\n", name, newname);
        if(extract_alias_residue_define(data->aliases, name, newname)) {
            PyErr_SetString(PyExc_ValueError, "failed on residue alias");
            return NULL;
        }
    } else if (!strcasecmp(type, "atom")) {
        if (!resname) {
            PyErr_SetString(PyExc_ValueError,
                            "resname must be provided when aliasing atoms");
            return NULL;
        }
        resname = strtoupper(name, data->all_caps);
        fprintf(data->outstream, "Aliasing residue %s atom %s to %s\n",
                resname, name, newname);
        if(extract_alias_atom_define(data->aliases, resname, name, newname)) {
            PyErr_SetString(PyExc_ValueError, "failed on atom alias");
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "alias type must be either 'atom' or 'residue'");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_set_allcaps(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "allcaps", NULL};
    PyObject *stateptr;
    psfgen_data *data;
    int allcaps = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|p:set_allcaps", kwnames,
                                     &stateptr, &allcaps)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    data->all_caps = allcaps ? 1 : 0;

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_regenerate(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "task", NULL};
    PyObject *stateptr;
    psfgen_data *data;
    char *task;
    int rc;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:regenerate", kwnames,
                                     &stateptr, &task)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    if (!strcasecmp(task, "angles")) {
        rc = topo_mol_regenerate_angles(data->mol);
    } else if (!strcasecmp(task, "dihedrals")) {
        rc = topo_mol_regenerate_dihedrals(data->mol);
    } else if (!strcasecmp(task, "resids")) {
        rc = topo_mol_regenerate_resids(data->mol);
    } else {
        PyErr_Format(PyExc_ValueError,
                     "regenerate must be [angles,resids,dihedrals], got '%s'",
                     task);
        return NULL;
    }
    if (rc) {
        PyErr_Format(PyExc_ValueError, "%s regeneration failed", task);
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/* IO functions */
static PyObject* py_write_namdbin(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "filename",
                       (char*) "velocity_filename", NULL};
    PyObject *velobj = NULL;
    FILE *velfile = NULL;
    PyObject *stateptr;
    psfgen_data *data;
    char *velfilename;
    char *filename;
    FILE *pfile;
    int rc;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|O:write_namdbin", kwnames,
                                     &stateptr, &filename, &velfilename)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // Handle velocity filename potentially being Py_None
    velfilename = (velobj && velobj != Py_None) ? as_charptr(velobj) : NULL;

    // Open files for writing, with some error checking
    pfile = fopen(filename, "wb");
    if (!pfile) {
        PyErr_Format(PyExc_ValueError, "Cannot open file '%s' for writing",
                     filename);
        return NULL;
    }

    if (velfilename) {
        velfile = fopen(velfilename, "wb");
        if (!velfile) {
            fclose(pfile);
            PyErr_Format(PyExc_ValueError, "Cannot open file '%s' for writing",
                         velfilename);
            return NULL;
        }
    }

    // Write the file, clean up, then errocr check
    rc = topo_mol_write_namdbin(data->mol, pfile, velfile,
                                data->outstream, python_msg);
    fclose(pfile);
    if (velfile)
        fclose(velfile);

    if (rc) {
        PyErr_Format(PyExc_ValueError, "Failed writing namdbin file '%s'",
                     filename);
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_write_pdb(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "filename", NULL};
    PyObject *stateptr;
    psfgen_data *data;
    char *filename;
    FILE *fd;
    int rc;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:write_pdb", kwnames,
                                     &stateptr, &filename)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    fd = fopen(filename, "w");
    if (!fd) {
        PyErr_Format(PyExc_OSError, "cannot open pdb file '%s' for writing",
                     filename);
        return NULL;
    }

    rc = topo_mol_write_pdb(data->mol, fd, data->outstream, python_msg);
    fclose(fd);
    if (rc) {
        PyErr_Format(PyExc_ValueError, "cannot write pdb '%s'", filename);
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_write_psf(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "filename",
                       (char*) "type", NULL};
    char *filename, *type;
    PyObject *stateptr;
    psfgen_data *data;
    int rc, charmmfmt;
    FILE *fd;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oss:write_psf", kwnames,
                                     &stateptr, &filename, &type)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    fd = fopen(filename, "w");
    if (!fd) {
        PyErr_Format(PyExc_OSError, "cannot open psf file '%s' for writing",
                     filename);
        return NULL;
    }
    if (!strcasecmp(type, "charmm")) {
        charmmfmt = 1;
    } else if (!strcasecmp(type, "x-plor")) {
        charmmfmt = 0;
    } else {
        PyErr_Format(PyExc_ValueError, "psf format '%s' not in [charmm,x-plor]",
                     type);
        return NULL;
    }

    rc = topo_mol_write_psf(data->mol, fd, charmmfmt, 0, 0,
                            data->outstream, python_msg);
    fclose(fd);
    if (rc) {
        PyErr_Format(PyExc_ValueError, "cannot write psf '%s'", filename);
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_read_psf(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "filename",
                       (char*) "pdbfile", (char*) "namdbinfile",
                       (char*) "velnamdbinfile", NULL};
    PyObject *pdbobj = NULL, *namdobj = NULL, *velobj = NULL;
    FILE *pdb = NULL, *psf = NULL, *namd = NULL, *vel = NULL;
    char *pdbfile = NULL, *namdfile = NULL, *velfile = NULL;
    PyObject *stateptr;
    psfgen_data *data;
    char *psffile;
    int rc;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|OOO:read_psf", kwnames,
                                     &stateptr, &psffile, &pdbobj,
                                     &namdobj, &velobj)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // Check arguments for pdbfile, namdbinfile, velnamdbinfile aren't Py_None
    pdbfile = (pdbobj && pdbobj != Py_None) ? as_charptr(pdbobj) : NULL;
    namdfile = (namdobj && namdobj != Py_None) ? as_charptr(namdobj) : NULL;
    velfile = (velobj && velobj != Py_None) ? as_charptr(velobj) : NULL;

    // Open files as psf_file_extract takes file pointers
    psf = fopen(psffile, "rb");
    if (pdbfile)
        pdb = fopen(pdbfile, "rb");
    if (namdfile)
        namd = fopen(namdfile, "rb");
    if (velfile)
        vel = fopen(velfile, "rb");

    // Error check all files at once, nicer code but less specific error msg
    if (!psf || (!pdb && pdbfile) || (!namd && namdfile) || (!vel && velfile)) {
        if (psf)
            fclose(psf);
        if (pdb)
            fclose(pdb);
        if (namd)
            fclose(namd);
        if (vel)
            fclose(vel);
        PyErr_Format(PyExc_ValueError, "Cannot open files for read psf '%s'",
                     psffile);
        return NULL;
    }

    // Actually do the parsing, and then clean up
    rc = psf_file_extract(data->mol, psf, pdb, namd, vel,
                          data->outstream, python_msg);
    fclose(psf);
    if (pdb)
        fclose(pdb);
    if (namd)
        fclose(namd);
    if (vel)
        fclose(vel);

    if (rc) {
        PyErr_Format(PyExc_ValueError, "Failed to parse psf file '%s'",
                     psffile);
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_read_coords(PyObject *self, PyObject *args, PyObject *kwargs) 
{
    char *kwnames[] = {(char*) "psfstate", (char*) "filename",
                       (char*) "segid", NULL};
    char *filename, *segid;
    PyObject *stateptr;
    psfgen_data *data;
    FILE *fd;
    int rc;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oss:read_coords", kwnames,
                                     &stateptr, &filename, &segid)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    fd = fopen(filename, "r");
    if (!fd) {
        PyErr_Format(PyExc_FileNotFoundError, "cannot open coordinate file '%s'",
                     filename);
        return NULL;
    }

    rc = pdb_file_extract_coordinates(data->mol, fd, NULL, segid, data->aliases,
                                      data->all_caps, data->outstream,
                                      python_msg);
    fclose(fd);
    if (rc) {
        PyErr_Format(PyExc_ValueError,
                     "cannot read coordinates '%s' into segment '%s'", filename,
                     segid);
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/* Segment functions */
static PyObject* py_add_segment(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "segid", (char*) "pdbfile",
                       (char*) "first", (char*) "last", (char*) "auto_angles",
                       (char*) "auto_dihedrals", (char*) "residues", (char*)
                       "mutate", NULL};
    char *first = NULL, *last = NULL, *filename = NULL;
    PyObject *mutate = NULL, *residues = NULL;
    int autoang = 1, autodih = 1;
    PyObject *stateptr;
    psfgen_data *data;
    char *segname;
    FILE *fd;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|sssppOO:add_segment",
                                     kwnames, &stateptr, &segname, &filename,
                                     &first, &last, &autoang, &autodih,
                                     &residues, &mutate)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) return NULL;

    // Sanity check segment name
    segname = strtoupper(segname, data->all_caps);
    if (strlen(segname) > 7) {
        PyErr_Format(PyExc_ValueError, "segment name '%s' more than 7 characters",
                     segname);
        return NULL;
    }
    if (topo_mol_segment(data->mol, segname))
        return NULL;

    // Set first and last, if present
    if (first) {
        if (topo_mol_segment_first(data->mol, first)) {
            PyErr_Format(PyExc_ValueError, "Cannot set first patch in segment "
                         "'%s' to '%s'", segname, first);
            return NULL;
        }
    }

    if (last) {
        if (topo_mol_segment_last(data->mol, last)) {
            PyErr_Format(PyExc_ValueError, "Cannot set last patch in segment "
                         "'%s' to '%s'", segname, last);
            return NULL;
        }
    }

    // Set auto angles and dihedrals
    if (topo_mol_segment_auto_angles(data->mol, autoang)) {
        PyErr_Format(PyExc_ValueError, "Failed setting angle autogen for "
                     "segment %s", segname);
        return NULL;
    }

    if (topo_mol_segment_auto_dihedrals(data->mol, autodih)) {
        PyErr_Format(PyExc_ValueError, "Failed setting dihedral autogen for "
                     "segment %s", segname);
        return NULL;
    }

    // If pdb file, do that before finishing the segment
    if (filename) {
        int rc;
        fd = fopen(filename, "r");
        if (!fd) {
            PyErr_Format(PyExc_FileNotFoundError,
                         "cannot open coordinate file '%s'", filename);
            return NULL;
        }
        rc = pdb_file_extract_residues(data->mol, fd, data->aliases,
                                       data->all_caps, data->outstream,
                                       python_msg);
        fclose(fd);
        if (rc) {
            PyErr_Format(PyExc_ValueError, "cannot read pdb file '%s'",
                    filename);
            return NULL;
        }
    }

    // Add residues to end, if desired
    if (residues && (residues != Py_None)) {
        char *resid, *resname, *chain;
        PyObject *residue;
        int n;

        if (!PyList_Check(residues)) {
            PyErr_SetString(PyExc_ValueError, "residues must be a list!");
            return NULL;
        }

        for (int i = 0; i < (int)PyList_Size(residues); i++) {
            residue = PyList_GetItem(residues, i);

            if (!PyTuple_Check(residue)) {
                PyErr_SetString(PyExc_ValueError,
                                "residues must be list of tuple");
                return NULL;
            }

            n = (int) PyTuple_Size(residue);
            if ((n != 2) && (n != 3)) {
                PyErr_SetString(PyExc_ValueError, "residues must be a list of "
                                "2 or 3 tuples");
                return NULL;
            }

            // Unpack tuple arguments, with chain being optional
            resid = as_charptr(PyTuple_GetItem(residue, 0));
            resname = as_charptr(PyTuple_GetItem(residue, 1));
            chain = (n == 3) ? as_charptr(PyTuple_GetItem(residue, 2)) : "";

            if (topo_mol_residue(data->mol, resid, resname, chain)) {
                    PyErr_Format(PyExc_ValueError,
                                 "Failed to add residue '%s:%s'",
                                 resname, resid);
                    return NULL;
            }
        }
    }

    // Do mutations, if provided
    if (mutate && (mutate != Py_None)) {
        char *resid, *resname;
        PyObject *residue;

        if (!PyList_Check(mutate)) {
            PyErr_SetString(PyExc_ValueError, "mutate must be a list!");
            return NULL;
        }

        for (int i = 0; i < (int)PyList_Size(mutate); i++) {
            residue = PyList_GetItem(mutate, i);
            if ( (!PyTuple_Check(residue)) || (int)PyTuple_Size(residue) != 2) {
                PyErr_SetString(PyExc_ValueError, "mutate must be a list of "
                                "2 or 3-tuples");
                return NULL;
            }

            // Unpack tuple
            resid = as_charptr(PyTuple_GetItem(residue, 0));
            resname = as_charptr(PyTuple_GetItem(residue, 1));

            if (topo_mol_mutate(data->mol, resid, resname)) {
                PyErr_Format(PyExc_ValueError,
                             "Failed to mutate residue '%s:%s'",
                             resname, resid);
                return NULL;
            }
        }
    }

    // Check result
    if (topo_mol_end(data->mol)) {
        PyErr_Format(PyExc_ValueError, "failed building segment '%s'", segname);
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_query_segment(PyObject *self, PyObject *args,
                                  PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "task", (char*) "segid",
                       (char*) "resid", NULL};
    char *segid = NULL, *resid = NULL;
    PyObject *stateptr, *objid;
    topo_mol_segment_t *seg;
    PyObject *result = NULL;
    psfgen_data *data;
    char *task;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|ss:query_segment", kwnames,
                                     &stateptr, &task, &segid, &resid)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // Ensure task argument is valid
    if (strcasecmp(task, "first") && strcasecmp(task, "last") \
     && strcasecmp(task, "resids") && strcasecmp(task, "residue") \
     && strcasecmp(task, "segids")) {
        PyErr_Format(PyExc_ValueError, "Unknown segment query '%s'", task);
        return NULL;
    }

    // Extract segid argument for non-segid requests
    if ((!segid) && strcasecmp(task, "segids")) {
            PyErr_Format(PyExc_ValueError, "segid argument must be passed for "
                         "segment task '%s'", task);
            return NULL;

    } else if (strcasecmp(task, "segids")) {
        // Identify the segment
        int segidx = (data->mol ? hasharray_index(data->mol->segment_hash, segid)
                                : HASHARRAY_FAIL);
        if (segidx == HASHARRAY_FAIL) {
            PyErr_Format(PyExc_ValueError,
                         "segid '%s' doesn't exist", segid);
            return NULL;
        }
        seg = data->mol->segment_array[segidx];
    }

    // Now build the result, depending on what the query is
    if (!strcasecmp(task, "segids")) {
        result = PyList_New(0);
        if (data->mol) {
            for (int i = 0; i < hasharray_count(data->mol->segment_hash); i++) {
                if (hasharray_index(data->mol->segment_hash,
                                    data->mol->segment_array[i]->segid)
                  != HASHARRAY_FAIL) {
                    objid = as_pystring(data->mol->segment_array[i]->segid);

                    if (PyList_Append(result, objid))
                        return NULL;
               }
           }
        }

    } else if (!strcasecmp(task, "first")) {
        result = as_pystring(seg->pfirst);

        // If patch is "none", return the Python None object
        if (!strcasecmp(seg->pfirst, "none")) {
            Py_INCREF(Py_None);
            result = Py_None;
        }

    } else if (!strcasecmp(task, "last")) {
        result = as_pystring(seg->plast);

        if (!strcasecmp(seg->plast, "none")) {
            Py_INCREF(Py_None);
            result = Py_None;
        }

    } else if (!strcasecmp(task, "resids")) {
        result = PyList_New(0);
        for (int i = 0; i < hasharray_count(seg->residue_hash); i++) {
            if (hasharray_index(seg->residue_hash, seg->residue_array[i].resid)
             != HASHARRAY_FAIL) {
                objid = as_pystring(seg->residue_array[i].resid);
                if (PyList_Append(result, objid))
                    return NULL;
            }
        }

    } else if (!strcasecmp(task, "residue")) {
        int residx;
        if (!resid) {
                PyErr_Format(PyExc_ValueError, "resid argument must be passed "
                             "for segment task '%s'", task);
                return NULL;
        }
        residx = hasharray_index(seg->residue_hash, resid);

        if (residx == HASHARRAY_FAIL) {
            PyErr_Format(PyExc_ValueError, "invalid resid '%s' for segment "
                         "'%s'", resid, segid);
            return NULL;
        }
        result = as_pystring(seg->residue_array[residx].name);
    }
    return result;
}

/* Module topology functions */
static PyObject* py_query_system(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "task", NULL};
    PyObject *stateptr, *result;
    psfgen_data *data = NULL;
    topo_defs *defs;
    char *task;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:query_system", kwnames,
                                     &stateptr, &task)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    result = PyList_New(0);
    if (!result)
        return NULL;

    defs = data->mol->defs;

    if (!strcasecmp(task, "topologies")) {
        for (int i = 0; i < hasharray_count(defs->topo_hash); i++) {
            if (PyList_Append(result, as_pystring(defs->topo_array[i].filename))) {
                PyErr_SetString(PyExc_ValueError, "cannot enumerate topologies");
                return NULL;
            }
        }

    } else if (!strcasecmp(task, "patches") || !strcasecmp(task, "residues")) {
        for (int i = 0; i < hasharray_count(defs->residue_hash); i++) {
            if ((!strcasecmp(task, "patches") && defs->residue_array[i].patch)
            ||  (!strcasecmp(task, "residues") && !defs->residue_array[i].patch)) {
                if (PyList_Append(result,
                                  as_pystring(defs->residue_array[i].name))) {
                    PyErr_SetString(PyExc_ValueError, "cannot enumerate residues");
                    return NULL;
                }
            }
        }

    } else {
        PyErr_Format(PyExc_ValueError, "Task '%s' invalid system query", task);
        return NULL;
    }
    return result;
}

static PyObject* py_parse_topology(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "filename", NULL};
    char *filename;
    PyObject *stateptr;
    psfgen_data *data = NULL;
    FILE *fd = NULL;
    int rc;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:parse_topology", kwnames,
                                     &stateptr, &filename)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    fd = fopen(filename, "r");
    if (!fd) {
        PyErr_Format(PyExc_FileNotFoundError, "cannot open topology file '%s'",
                     filename);
        return NULL;
    }
    rc = charmm_parse_topo_defs(data->defs, fd, data->all_caps,
                                data->outstream, python_msg);
    fclose(fd);
    if (rc) {
        PyErr_Format(PyExc_ValueError, "error parsing topology file '%s'",
                    filename);
        return NULL;
    }

    topo_defs_add_topofile(data->defs, filename);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_get_patches(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "listall", NULL};
    topo_mol_patchres_t *patchres;
    PyObject *stateptr, *result;
    topo_mol_patch_t *patch;
    psfgen_data *data;
    int listall = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|p:get_patches", kwnames,
                                     &stateptr, &listall)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    result = PyList_New(0);
    if (!result)
        return NULL;

    for (patch = data->mol->patches; patch; patch = patch->next) {

        PyObject *patchtuple = NULL;
        if (patch->deflt && !listall)
            continue;

        // Find resids this patch is applied to
        for (patchres = patch->patchresids; patchres; patchres = patchres->next) {
            // If invalid patch for this residue, keep going
            if (!topo_mol_validate_patchres(data->mol, patch->pname,
                                            patchres->segid, patchres->resid)) {
                break;
            }

            // Generate a tuple of (patchname, segid, resid)
            patchtuple = PyTuple_Pack(3, as_pystring(patch->pname),
                                      as_pystring(patchres->segid),
                                      as_pystring(patchres->resid));
            if (!patchtuple || PyList_Append(result, patchtuple)) {
                PyErr_Format(PyExc_ValueError, "Patch %s %s:%s failed",
                             patch->pname, patchres->segid, patchres->resid);
                return NULL;
            }
        }
    }
    return result;
}

static PyObject* py_patch(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "patchname",
                       (char*) "targets", NULL};
    PyObject *stateptr, *targlist, *target;
    topo_mol_ident_t *targets;
    psfgen_data *data;
    int rc, ntargets;
    char *patchname;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OsO:patch", kwnames,
                                     &stateptr, &patchname, &targlist)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // If targets are lists, just turn them into tuples here
    if (PyList_Check(targlist))
        targlist = PyList_AsTuple(targlist);

    if (!PyTuple_Check(targlist)) {
        PyErr_SetString(PyExc_ValueError,
                        "patch targets must be a list or tuple!");
        return NULL;
    }
    ntargets = (int)PyTuple_Size(targlist);

    // Construct targets
    // It's a list/tuple of list/tuples
    targets = malloc(ntargets*sizeof(topo_mol_ident_t));

    for (int i = 0; i < ntargets; ++i) {
        target = PyTuple_GetItem(targlist, i);

        if (PyList_Check(target)) {
            target = PyList_AsTuple(target);
        }
        if (!PyTuple_Check(target) || (int)PyTuple_Size(target) != 2) {
            PyErr_SetString(PyExc_ValueError,
                            "patch target must be a list or tuple (segid, resid)");
            return NULL;
        }

        targets[i].segid = as_charptr(PyTuple_GetItem(target, 0));
        targets[i].resid = as_charptr(PyTuple_GetItem(target, 1));
        targets[i].aname = NULL;

        if (PyErr_Occurred()) {
            free(targets);
            return NULL;
        }
    }

    // Actually do the work
    rc = topo_mol_patch(data->mol, targets, ntargets, patchname, 0, 0, 0, 0);
    free(targets);
    if (rc) {
        PyErr_Format(PyExc_ValueError, "Cannot apply patch %s", patchname);
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_query_atoms(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "segid", (char*) "resid",
                       (char*) "task", NULL};
    PyObject *stateptr, *result, *atomresult;
    char *segid, *resid, *task;
    topo_mol_segment_t *seg;
    topo_mol_atom_t *atoms;
    int segidx, residx;
    psfgen_data* data;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Osss:query_atoms", kwnames,
                                     &stateptr, &segid, &resid, &task)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // Get segment index from segid
    segidx = hasharray_index(data->mol->segment_hash, segid);
    if (segidx == HASHARRAY_FAIL) {
        PyErr_Format(PyExc_ValueError, "Segment '%s' does not exist", segid);
        return NULL;
    }
    seg = data->mol->segment_array[segidx];

    // Get residue index from segid and resid
    residx = hasharray_index(seg->residue_hash, resid);
    if (residx == HASHARRAY_FAIL) {
        PyErr_Format(PyExc_ValueError, "No resid '%s' in segment '%s'",
                     resid, segid);
        return NULL;
    }

    // Loop through atoms in residue and add names to list
    result = PyList_New(0);
    atoms = seg->residue_array[residx].atoms;
    while (atoms) {
        if (!strcasecmp(task, "name")) {
            atomresult = as_pystring(atoms->name);

        } else if (!strcmp(task, "coordinates")) {
            atomresult = PyTuple_Pack(3, PyFloat_FromDouble(atoms->x),
                                      PyFloat_FromDouble(atoms->y),
                                      PyFloat_FromDouble(atoms->z));

        } else if (!strcmp(task, "velocities")) {
            atomresult = PyTuple_Pack(3, PyFloat_FromDouble(atoms->vx),
                                      PyFloat_FromDouble(atoms->vy),
                                      PyFloat_FromDouble(atoms->vz));

        } else if (!strcmp(task, "mass")) {
            atomresult = PyFloat_FromDouble(atoms->mass);

        } else if (!strcmp(task, "charge")) {
            atomresult = PyFloat_FromDouble(atoms->charge);

        } else if (!strcmp(task, "atomid")) {
            atomresult = as_pyint(atoms->atomid);

        } else {
            PyErr_Format(PyExc_ValueError, "invalid atom task '%s'", task);
            return NULL;
        }

        if (PyList_Append(result, atomresult)) {
            PyErr_SetString(PyExc_ValueError, "cannot gather atoms");
            return NULL;
        }
        atoms = atoms->next;
    }
    return result;
}

static PyObject* py_delete_atoms(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "segid", (char*) "resid",
                       (char*) "aname", NULL};
    PyObject *resobj = NULL, *atomobj = NULL;
    char *resid = NULL, *aname = NULL;
    topo_mol_ident_t target;
    PyObject *stateptr;
    psfgen_data *data;
    char *segid;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|OO:delete_atoms", kwnames,
                                     &stateptr, &segid, &resobj, &atomobj)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // Handle None for resid and atom name fields
    if (resobj && resobj != Py_None)
        resid = as_charptr(resobj);

    if (atomobj && atomobj != Py_None)
        aname = as_charptr(atomobj);

    // Build target object
    target.segid = segid;
    target.resid = resid;
    target.aname = aname;

    if (topo_mol_delete_atom(data->mol, &target)) {
        PyErr_SetString(PyExc_ValueError, "failed to delete atoms");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_set_coord(PyObject *self, PyObject *args, PyObject *kwargs)
{
    char *kwnames[] = {(char*) "psfstate", (char*) "segid", (char*) "resid",
                       (char*) "aname", (char*) "position", NULL};
    PyObject *position, *stateptr;
    char *segid, *aname, *resid;
    topo_mol_ident_t target;
    psfgen_data *data;
    double x, y, z;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OsssO:set_coords", kwnames,
                                     &stateptr, &segid, &resid, &aname,
                                     &position)) {
        return NULL;
    }

    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
        return NULL;

    // Unpack position tuple to x, y, z
    if (!PyTuple_Check(position) || PyTuple_Size(position) != 3) {
        PyErr_SetString(PyExc_ValueError, "position must be a 3-tuple");
        return NULL;
    }
    x = PyFloat_AsDouble(PyTuple_GetItem(position, 0));
    y = PyFloat_AsDouble(PyTuple_GetItem(position, 1));
    z = PyFloat_AsDouble(PyTuple_GetItem(position, 2));
    if (PyErr_Occurred()) {
        return NULL;
    }

    // Build target object
    target.segid = segid;
    target.resid = resid;
    target.aname = aname;

    // Do the work
    if(topo_mol_set_xyz(data->mol, &target, x, y, z)) {
        PyErr_SetString(PyExc_ValueError, "failed to set coordinates");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_guess_coords(PyObject *self, PyObject *stateptr)
{
    psfgen_data* data;

    // Unpack molecule capsule
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred())
       return NULL;

    if (topo_mol_guess_xyz(data->mol)) {
        PyErr_SetString(PyExc_ValueError, "failed to guess coordinates");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/* Method definitions */
static PyMethodDef methods[] = {
    {(char *) "set_coord", (PyCFunction)py_set_coord, METH_VARARGS | METH_KEYWORDS},
    {(char *) "init_mol", (PyCFunction)py_init_mol, METH_VARARGS | METH_KEYWORDS},
    {(char *) "del_mol", (PyCFunction)py_del_mol, METH_O},
    {(char *) "alias_residue", (PyCFunction)py_alias, METH_VARARGS | METH_KEYWORDS},
    {(char *) "query_system", (PyCFunction)py_query_system, METH_VARARGS | METH_KEYWORDS},
    {(char *) "query_segment", (PyCFunction)py_query_segment, METH_VARARGS | METH_KEYWORDS},
    {(char *) "parse_topology", (PyCFunction)py_parse_topology, METH_VARARGS | METH_KEYWORDS},
    {(char *) "read_coords", (PyCFunction)py_read_coords, METH_VARARGS | METH_KEYWORDS},
    {(char *) "patch", (PyCFunction)py_patch, METH_VARARGS | METH_KEYWORDS},
    {(char *) "regenerate", (PyCFunction)py_regenerate, METH_VARARGS | METH_KEYWORDS},
    {(char *) "write_psf", (PyCFunction)py_write_psf, METH_VARARGS | METH_KEYWORDS},
    {(char *) "write_pdb", (PyCFunction)py_write_pdb, METH_VARARGS | METH_KEYWORDS},
    {(char *) "add_segment", (PyCFunction)py_add_segment, METH_VARARGS | METH_KEYWORDS},
    {(char *) "guess_coords", (PyCFunction)py_guess_coords, METH_O},
    {(char *) "get_patches", (PyCFunction)py_get_patches, METH_VARARGS | METH_KEYWORDS},
    {(char *) "query_atoms", (PyCFunction)py_query_atoms, METH_VARARGS | METH_KEYWORDS},
    {(char *) "delete_atoms", (PyCFunction)py_delete_atoms, METH_VARARGS | METH_KEYWORDS},
    {(char *) "set_allcaps", (PyCFunction)py_set_allcaps, METH_VARARGS | METH_KEYWORDS},
    {NULL, NULL, 0, NULL}
};

/* Module initialization functions
 * Python 2 and 3 are totally separate here for clarity
 */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef psfgendef = {
    PyModuleDef_HEAD_INIT,
    "_psfgen",
    NULL,
    -1,
    methods,
};

PyMODINIT_FUNC PyInit__psfgen(void)
{
    PyObject *m = PyModule_Create(&psfgendef);
    return m;
}
#else
PyMODINIT_FUNC init_psfgen(void)
{
    PyObject *m = Py_InitModule("_psfgen", methods);
}
#endif

