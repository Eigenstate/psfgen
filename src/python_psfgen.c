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
void python_msg(void *v, const char *msg) {
    v = (FILE*)v; // print to this file
    fprintf(v, msg);
    fprintf(v, "\n");
}

/* Initialization / destruction functions */
static PyObject* py_init_mol(PyObject *self) {
    psfgen_data * data = malloc(sizeof(psfgen_data));

    // Initialize topologies
    data->defs = topo_defs_create();
    //topo_defs_error_handler(data->defs, inter??);

    // Initialize aliases
    data->aliases = stringhash_create();
    data->mol = topo_mol_create(data->defs);

    // Initialize other stuffs
    data->id = 0; // Doesn't matter since data is per class instance
    data->in_use = 0;
    data->all_caps = 1;
    data->outstream = stdout; // TODO: does this even work
    topo_mol_error_handler(data->mol, data->outstream, python_msg);

    // Encapsulate psf_data object and return it
    PyObject * capsule = PyCapsule_New(data, NULL, NULL);
    if (!capsule || PyErr_Occurred()) {
        return NULL;
    }
    return capsule;
}

static PyObject* py_del_mol(PyObject *self, PyObject *stateptr) {
    psfgen_data* data;

    // Unpack molecule capsule
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) {
       return NULL;
    }

    // Do the work
    topo_mol_destroy(data->mol);
    topo_defs_destroy(data->defs);
    stringhash_destroy(data->aliases);
    free(data);

    Py_INCREF(Py_None); // Don't mess up reference counts
    return Py_None;
}

/* Aliases and names and stuff */
static PyObject* py_alias(PyObject *self, PyObject *args, PyObject *kwargs) {
    char *type, *name, *newname;
    char *resname = NULL;
    PyObject *stateptr;
    psfgen_data *data;
    int rc;

    // Parse arguments and put them in the right format
    char *kwnames[] = {(char*) "psfstate", (char*) "type", (char*) "name",
                       (char*) "newname", (char*) "resname", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Osss|s:alias", kwnames,
                                     &stateptr, &type, &name, &newname,
                                     &resname)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

    name = strtoupper(name, data->all_caps);
    newname = strtoupper(newname, data->all_caps);

    if (!strcasecmp(type, "residue")) {
        fprintf(data->outstream, "Aliasing residue %s to %s\n", name, newname);
        rc = extract_alias_residue_define(data->aliases, name, newname);
        if (rc) {
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
        rc = extract_alias_atom_define(data->aliases, resname, name, newname);
        if (rc) {
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

static PyObject* py_regenerate(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *stateptr;
    char *task;
    psfgen_data *data;
    int rc;
    char *kwnames[] = {(char*) "psfstate", (char*) "task", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:regenerate", kwnames,
                                     &stateptr, &task)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

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
static PyObject* py_write_pdb(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *stateptr;
    char *filename;
    FILE *fd;
    psfgen_data *data;
    int rc;

    char *kwnames[] = {(char*) "psfstate", (char*) "filename", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:write_pdb", kwnames,
                                     &stateptr, &filename)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

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

static PyObject* py_write_psf(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *stateptr;
    char *filename, *type;
    FILE *fd;
    psfgen_data *data;
    int rc, charmmfmt;

    char *kwnames[] = {(char*) "psfstate", (char*) "filename",
                       (char*) "type", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oss:write_psf", kwnames,
                                     &stateptr, &filename, &type)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

    fd = fopen(filename, "w");
    if (!fd) {
        PyErr_Format(PyExc_OSError, "cannot open psf file '%s' for writing",
                     filename);
        return NULL;
    }
    if (!strcmp(type, "charmm")) {
        charmmfmt = 1;
    } else if (!strcmp(type, "x-plor")) {
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

static PyObject* py_read_coords(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *stateptr;
    char *filename, *segid;
    int rc;
    psfgen_data *data;
    FILE *fd;
    char *kwnames[] = {(char*) "psfstate", (char*) "filename",
                       (char*) "segid", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oss:read_coords", kwnames,
                                     &stateptr, &filename, &segid)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

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
static PyObject* py_add_segment(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *stateptr;
    char *segname;
    char *filename = NULL;
    FILE *fd;
    psfgen_data *data;
    int rc;

    char *kwnames[] = {(char*) "psfstate", (char*) "segid", (char*) "pdbfile",
                        NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|s:add_segment", kwnames,
                                     &stateptr, &segname, &filename)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

    // Sanity check segment name
    segname = strtoupper(segname, data->all_caps);
    if (strlen(segname) > 7) {
        PyErr_Format(PyExc_ValueError, "segment name '%s' more than 7 characters",
                     segname);
        return NULL;
    }
    rc = topo_mol_segment(data->mol, segname);
    if (rc) { return NULL; }

    // If pdb file, do that before finishing the segment
    if (filename) {
        fd = fopen(filename, "r");
        if (!fd) {
            PyErr_Format(PyExc_FileNotFoundError,
                         "cannot open coordinate file '%s'", filename);
            return NULL;
        }
        rc = pdb_file_extract_residues(data->mol, fd, data->aliases, data->all_caps,
                                       data->outstream, python_msg);
        fclose(fd);
        if (rc) {
            PyErr_Format(PyExc_ValueError, "cannot read pdb file '%s'",
                    filename);
            return NULL;
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

static PyObject* py_get_segment(PyObject *self, PyObject *args, PyObject *kwargs) {
    // TODO: handle atoms arguments requests
    PyObject *stateptr, *objid;
    char *task;
    char *segid = NULL;
    char *resid = NULL;
    PyObject *result = NULL;
    psfgen_data *data;
    int rc;
    char *kwnames[] = {(char*) "psfstate", (char*) "task", (char*) "segid",
                       (char*) "resid", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|ss:get_segment", kwnames,
                                     &stateptr, &task, &segid, &resid)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) { return NULL; }

    // Ensure task argument is valid
    if (strcasecmp(task, "first") && strcasecmp(task, "last") \
     && strcasecmp(task, "resids") && strcasecmp(task, "residue") \
     && strcasecmp(task, "segids")) {
        PyErr_SetString(PyExc_ValueError,
                        "segment task must be in [first,last,resids,residue]");
        return NULL;
    }

    if (!strcasecmp(task, "segids")) {
        result = PyList_New(0);
        if (data->mol) {
           for (int i=0; i<hasharray_count(data->mol->segment_hash); ++i) {
#if PY_MAJOR_VERSION >=3
               objid = PyUnicode_FromString(data->mol->segment_array[i]->segid);
#else
               objid = PyString_FromString(data->mol->segment_array[i]->segid);
#endif
               rc = PyList_Append(result, objid);
               if (rc) { return NULL; }
           }
        }
        return result;
    }

    // Require segid argument for non-segid requests
    if (!segid) {
            PyErr_Format(PyExc_ValueError,
                        "segid argument must be passed for segment task '%s'",
                        task);
            return NULL;
    }

    // Identify the segment
    int segidx = (data->mol ? hasharray_index(data->mol->segment_hash, segid)
                            : HASHARRAY_FAIL);
    if (segidx == HASHARRAY_FAIL) {
        PyErr_Format(PyExc_ValueError,
                     "segid '%s' doesn't exist", segid);
        return NULL;
    }
    topo_mol_segment_t *seg = data->mol->segment_array[segidx];

    if (!strcasecmp(task, "first")) {
#if PY_MAJOR_VERSION >=3
        result = PyUnicode_FromString(seg->pfirst);
#else
        result = PyUnicode_FromString(seg->pfirst);
#endif
    } else if (!strcasecmp(task, "last")) {
#if PY_MAJOR_VERSION >=3
        result = PyUnicode_FromString(seg->plast);
#else
        result = PyString_FromString(seg->plast);
#endif
    } else if (!strcasecmp(task, "resids")) {
        result = PyList_New(0);
        for (int i=0; i<hasharray_count(seg->residue_hash); ++i) {
            if (hasharray_index(seg->residue_hash, seg->residue_array[i].resid)
                    != HASHARRAY_FAIL) {
#if PY_MAJOR_VERSION >=3
            objid = PyUnicode_FromString(seg->residue_array[i].resid);
#else
            objid = PyString_FromString(seg->residue_array[i].resid);
#endif
            rc = PyList_Append(result, objid);
            if (rc) { return NULL; }
            }
        }
    } else if (!strcasecmp(task, "residue")) {
        if (!resid) {
                PyErr_Format(PyExc_ValueError,
                            "resid argument must be passed for segment task '%s'",
                            task);
                return NULL;
        }
        int residx = hasharray_index(seg->residue_hash, resid);
        if (residx == HASHARRAY_FAIL) {
            PyErr_Format(PyExc_ValueError,
                         "invalid resid '%s' for segment '%s'", resid, segid);
            return NULL;
        }
#if PY_MAJOR_VERSION >=3
        result = PyUnicode_FromString(seg->residue_array[residx].name);
#else
        result = PyString_FromString(seg->residue_array[residx].name);
#endif
    }
    return result;
}

/* Module topology functions */
static PyObject* py_parse_topology(PyObject *self, PyObject *args, PyObject *kwargs) {
    char *filename;
    PyObject *stateptr;
    psfgen_data *data = NULL;
    FILE *fd = NULL;
    int rc;

    char *kwnames[] = {(char*) "psfstate", (char*) "filename", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os:parse_topology", kwnames,
                                     &stateptr, &filename)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) {
        return NULL;
    }

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

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* py_patch(PyObject *self, PyObject *args, PyObject *kwargs) {
    char *patchname;
    PyObject *stateptr, *targlist, *target;
    psfgen_data* data = NULL;
    topo_mol_ident_t *targets;
    int rc, ntargets;
    char *kwnames[] = {(char*) "psfstate", (char*) "patchname",
                       (char*) "targets", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OsO:patch", kwnames,
                                     &stateptr, &patchname, &targlist)) {
        return NULL;
    }
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) {
        return NULL;
    }
    // If targets are lists, just turn them into tuples here
    if (PyList_Check(targlist)) {
        targlist = PyList_AsTuple(targlist);
    }

    if (!PyTuple_Check(targlist)) {
        PyErr_SetString(PyExc_ValueError,
                        "patch targets must be a list or tuple!");
        return NULL;
    }
    ntargets = (int)PyTuple_Size(targlist);

    // Construct targets
        // It's a list/tuple of list/tuples
    targets = malloc(ntargets*sizeof(topo_mol_ident_t));
    for (int i=0; i<ntargets; ++i) {
        target = PyTuple_GetItem(targlist, i);
        if (PyList_Check(target)) {
            target = PyList_AsTuple(target);
        }
        if (!PyTuple_Check(target) || (int)PyTuple_Size(target) != 2) {
            PyErr_SetString(PyExc_ValueError,
                            "patch target must be a list or tuple (segid, resid)");
            return NULL;
        }
#if PY_MAJOR_VERSION >=3
        targets[i].segid = PyUnicode_AsUTF8(PyTuple_GetItem(target, 0));
        targets[i].resid = PyUnicode_AsUTF8(PyTuple_GetItem(target, 1));
#else
        targets[i].segid = PyString_AsString(PyTuple_GetItem(target, 0));
        targets[i].resid = PyString_AsString(PyTuple_GetItem(target, 1));
#endif
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

static PyObject* py_set_coords(PyObject *self, PyObject *args, PyObject *kwargs) {
    // Arguments will be parsed into here
    char *segid, *aname, *resid;
    PyObject *position, *stateptr;

    // Will be unpacked into these things
    topo_mol_ident_t target;
    psfgen_data* data = NULL;
    double x, y, z;
    int rc;

    // segid, resid, aname, position (char*, int, char*, tuple)
    char *kwnames[] = {(char*) "psfstate", (char*) "segid", (char*) "resid",
                       (char*) "aname", (char*) "position", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OsssO:set_coords", kwnames,
                                     &stateptr, &segid, &resid, &aname,
                                     &position)) {
        return NULL;
    }
    // Unpack molecule capsule
    data = PyCapsule_GetPointer(stateptr, NULL);
    if (!data || PyErr_Occurred()) {
        return NULL;
    }

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
    rc = topo_mol_set_xyz(data->mol, &target, x, y, z);
    if (rc) {
        PyErr_SetString(PyExc_ValueError, "failed to set coordinates");
        return NULL;
    }

    // TODO: python 2
    return PyLong_FromLong(rc);
}

/* Method definitions */
static PyMethodDef methods[] = {
    {(char *) "set_coords", (PyCFunction)py_set_coords, METH_VARARGS | METH_KEYWORDS},
    {(char *) "init_mol", (PyCFunction)py_init_mol, METH_NOARGS},
    {(char *) "del_mol", (PyCFunction)py_del_mol, METH_O},
    {(char *) "alias_residue", (PyCFunction)py_alias, METH_VARARGS | METH_KEYWORDS},
    {(char *) "get_segment", (PyCFunction)py_get_segment, METH_VARARGS | METH_KEYWORDS},
    {(char *) "parse_topology", (PyCFunction)py_parse_topology, METH_VARARGS | METH_KEYWORDS},
    {(char *) "read_coords", (PyCFunction)py_read_coords, METH_VARARGS | METH_KEYWORDS},
    {(char *) "patch", (PyCFunction)py_patch, METH_VARARGS | METH_KEYWORDS},
    {(char *) "regenerate", (PyCFunction)py_regenerate, METH_VARARGS | METH_KEYWORDS},
    {(char *) "write_psf", (PyCFunction)py_write_psf, METH_VARARGS | METH_KEYWORDS},
    {(char *) "write_pdb", (PyCFunction)py_write_pdb, METH_VARARGS | METH_KEYWORDS},
    {(char *) "add_segment", (PyCFunction)py_add_segment, METH_VARARGS | METH_KEYWORDS},
    {NULL, NULL, 0, NULL}
};

/* Module initialization functions */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef psfgendef = {
    PyModuleDef_HEAD_INIT,
    "_psfgen",
    NULL,
    -1,
    methods,
};
#endif

#if PY_MAJOR_VERSION >=3
PyMODINIT_FUNC PyInit__psfgen(void) {
    PyObject *m = PyModule_Create(&psfgendef);
#else
PyObject* init_psfgen(void) {
    PyObject *m = Py_InitModule((char *)"_psfgen", methods);
#endif

    return m;
}
