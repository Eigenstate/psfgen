
import _psfgen
import sys

# Definitions for psf file types
CHARMM="charmm"
XPLOR="x-plor"

# Definitions for auto arguments
AUTO_ANGLES, AUTO_DIHEDRALS, AUTO_NONE = "angles", "dihedrals", "none"

class PsfGen:

    """
    A Psf generator object. Represents the state of a single molecular system,
    with its own set of loaded topologies, residue and/or atom aliases, segments,
    and coordinates.
    """

    def __init__(self, output=sys.stdout):
        """
        Creates a PsfGen object.

        Args:
            outstream (str or stream object): Where to write output. Defaults
                to stdout. Must have fileno() attribute, or if a string, file
                will be opened.
        """
        if isinstance(output, str):
            self.output = open(output, 'wb')
        elif hasattr(output, "fileno"):
            self.output = output
        else:
            raise ValueError("output argument must be a str or open file")

        self._fileno = self.output.fileno()
        self._data = _psfgen.init_mol(outfd=self._fileno)

    def __del__(self):
        self.output.close()
        _psfgen.del_mol(self._data);

    #===========================================================================

    def read_topology(self, filename):
        """
        Parses a charmm format topology file into current library.

        Args:
            filename (str): File to parse

        Raises:
            FileNotFoundError: If the file cannot be opened for reading
            ValueError: If an error occurs during parsing
        """
        _psfgen.parse_topology(psfstate=self._data, filename=filename)

    #===========================================================================

    def alias_residue(self, top_resname, pdb_resname):
        """
        Provide an alternate name for a residue in the topology file compared
        to the PDB file.

        Args:
            top_resname (str): Resname in topology file
            pdb_resname (str): Equivalent resname in PDB files
        """
        _psfgen.alias(psfstate=self._data, type="residue",
                      name=top_resname, newname=pdb_resname)

    #===========================================================================

    def alias_atom(self, resname, top_atomname, pdb_atomname):
        """
        Provide an alternate name for a atom in a residue in the topology file
        compared to the PDB file.

        Args:
            resname (str): Residue name in which to make alias
            top_atomname (str): Atom name in the topology file
            pdb_atomname (str): Equivalent atom name in PDB file
        """
        _psfgen.alias(psfstate=self._data, type="atom",
                      resname=resname, name=top_atomname, newname=pdb_atomname)

    #===========================================================================

    def get_segids(self):
        """
        Obtain all the currently defined segment IDs.

        Returns:
            (list of str): All defined segids in current molecule
        """
        return _psfgen.query_segment(psfstate=self._data, task="segids")

    #===========================================================================

    def get_resids(self, segid):
        """
        Obtain all currently defined resids in a given segment

        Args:
            segid (str): Segment ID to query

        Returns:
            (list of str): All defined resids in given segment
        """
        return _psfgen.query_segment(psfstate=self._data, task="resids",
                                     segid=segid)

    #===========================================================================

    def get_patches(self, list_defaults=True, list_all=False):
        """
        Obtain information about available or applied patches in the current
        system state.

        Args:
            list_defaults (bool): If True, default patches will be listed too.
                Otherwise only explicitly applied patches will be shown.
                Defaults to True
            list_all (bool): List all available patches, not just applied ones.
                For psfgen internal reasons, "NONE", "None", and "none" will
                appear in this list.

        Returns:
            (list of 3-tuple): (patchname, segid, resid) of all applied patches
        """
        if list_all:
            return _psfgen.query_system(self._data, task="patches")

        return _psfgen.get_patches(self._data, listall=list_defaults)

    #===========================================================================

    def get_resname(self, segid, resid):
        """
        Obtains the residue name given a resid and segment

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (str): Residue name
        """
        # Handle integer resid, just turn it into a string
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.query_segment(psfstate=self._data, task="residue",
                                     segid=segid, resid=resid)

    #===========================================================================

    def get_atoms(self, segid, resid):
        """
        Obtains atom names in a given residue

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (list of str): Atom names in residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.get_atoms(psfstate=self._data, segid=segid, resid=resid)
    #===========================================================================

    def get_first(self, segid):
        """
        Get the name of the patch applied to the beginning of a given segment

        Args:
            segid (str): Segment ID to query

        Returns:
            (str): Patch name, or None
        """
        return _psfgen.query_segment(psfstate=self._data, task="first",
                                     segid=segid)

    #===========================================================================

    def get_last(self, segid):
        """
        Get the name of the patch applied to the end of a given segment

        Args:
            segid (str): Segment ID to query

        Returns:
            (str): Patch name, or None
        """
        return _psfgen.query_segment(psfstate=self._data, task="last",
                                    segid=segid)

    #===========================================================================

    def get_topologies(self):
        """
        Get all loaded topology files

        Returns:
            (list of str): Filenames that have been loaded
        """
        return _psfgen.query_system(psfstate=self._data, task="topologies")

    #===========================================================================

    def get_residue_types(self):
        """
        Get all defined residues

        Returns:
            (list of str): Defined residue types from current topologies
        """
        return _psfgen.query_system(psfstate=self._data, task="residues")

    #===========================================================================

    def add_segment(self, segid, first="none", last="none", pdbfile=None,
                    auto_angles=True, auto_dihedrals=True,
                    residues=None, mutate=None):
        """
        Adds a new segment to the internal molecule state.

        Args:
            segid (str): Name/ID of the new segment
            first (str): Patch to apply to first residue in segment, or
                None for topology file default setting
            last (str): Patch to apply to last residue in segment, or
                None for topology file in default setting
            pdbfile (str): PDB file to read residue information from, or
                None to create an empty segment
            auto_angles (bool): If angles should be autogenerated.
            auto_dihedrals (bool): If dihedrals should be autogenerated.
            residues (list of 2 or 3-tuple): (resid, residue, [chain]) of
                residues to append to the end of the segment. Chain can be None
                or unset in this tuple to use the current chain.
            mutate (list of 2 tuple): (resid, resname) of residues to alter.
                The given residue IDs will be set to the given residue name.
        """
        if self.get_segids() and segid in self.get_segids():
            raise ValueError("Duplicate segID '%s'" % segid)

        _psfgen.add_segment(psfstate=self._data, segid=segid, pdbfile=pdbfile,
                            first=first, last=last, auto_angles=auto_angles,
                            auto_dihedrals=auto_dihedrals, residues=residues,
                            mutate=mutate)

    #===========================================================================

    def read_coords(self, filename, segid):
        """
        Reads in coordinates from a PDB file, matching segment, residue, and
        atom names to the current segment.

        Args:
            filename (str): Filename of PDB file to read
            segid (str): Segment ID to assign coordinates to
        """
        if self.get_segids() and segid not in self.get_segids():
            raise ValueError("Can't read coordinates for segment '%s' as "
                             "it is undefined." % segid)
        _psfgen.read_coords(psfstate=self._data,
                            filename=filename,
                            segid=segid)

    #===========================================================================

    def set_coord(self, segid, resid, atomname, position):
        """
        Sets the coordinates of a given atom to new values.

        Args:
            segid (str): Segment ID of atom
            resid (str): Residue ID of atom
            atomname (str): Atom name
            position (3-tuple of double): New x, y, and z coordinates for atom

        Returns:
            (int): Number of coordinates set
        """
        return _psfgen.set_coord(psfstate=self._data, segid=segid,
                                 resid=resid, aname=atomname, position=position)

    #===========================================================================

    def guess_coords(self):
        """
        Sets unset coordinates using geometric assumptions
        """
        _psfgen.guess_coords(self._data)

    #===========================================================================

    def patch(self, patchname, targets):
        """
        Applies a patch to the molecule.

        Args:
            patchname (str): Name of the patch to apply
            targets (list of 2 tuple): (segid, resid) to apply patch to
        """

        _psfgen.patch(psfstate=self._data, patchname=patchname,
                      targets=targets)

    #===========================================================================

    #def multiply(self, ncopies, groups):
    #    identlist = []
    #    for group in groups:
    #        if len(group) == 2:
    #            identlist.append({"segid": str(group[0]),
    #                              "resid": str(group[1])})
    #        elif len(group) == 3:
    #            identlist.append({"segid": str(group[0]),
    #                              "resid": str(group[1]),
    #                              "aname":str(group[2])})
    #        else:
    #            raise ValueError("groups must contain lists of 2 or 3 elements")
    #    rc = topo_mol_multiply_atoms(mol=self.mol, targets=identlist,
    #                                 ncopies=ncopies)
    #    if rc:
    #        raise ValueError("failed to multiply atoms (error=%d)" % rc)

    #===========================================================================

    #def delAtom(self, segid, resid=None, atomName=None):
    #    ident={"segid": str(segid).encode()}
    #    if resID is not None:
    #        ident["resid"] = str(resid).encode()
    #    if atomName is not None:
    #        ident["aname"] = str(atomName).encode()
    #    topo_mol_delete_atom(mol=self.mol, target=ident)

    #===========================================================================

    def regenerate_angles(self):
        """
        Removes angles and regenerates them from bonds. Can be used after
        patching.
        """
        _psfgen.regenerate(self._data, task="angles")

    #===========================================================================

    def regenerate_dihedrals(self):
        """
        Removes dihedrals and regenerates them from angles. Can be used after
        patching. Usually, you should call `regenerate_angles` first.
        """
        _psfgen.regenerate(self._data, task="dihedrals")

    #===========================================================================

    def regenerate_resids(self):
        """
        Regenerates residue IDs by removing insertion codes and minimially
        modifying them for uniqueness.
        """
        _psfgen.regenerate(self._data, task="resids")

    #===========================================================================

    #def readPSF(self, filename):
    #    fd = fopen(filename, 'r')
    #    retval = psf_file_extract(self.mol, fd, None, None)
    #    fclose(fd)
    #    if retval:
    #        raise PsfgenFormatError("Error reading psf file '%s'"
    #                                % filename)

    #===========================================================================

    def write_psf(self, filename, type=XPLOR):
        """
        Writes the current molecule state out as a psf file

        Args:
            filename (str): Filename to write to
            type (str): Type of psf file to write, in ["charmm", "x-plor"]
        """
        _psfgen.write_psf(psfstate=self._data, filename=filename, type=type)

    #===========================================================================

    def write_pdb(self, filename):
        """
        Writes the current molecule state out as a pdb file

        Args:
            filename (str): Filename to write to
        """
        _psfgen.write_pdb(psfstate=self._data, filename=filename)

    #===========================================================================

    def query_segment(self, task, segid=None, resid=None):
        """
        Ask for information about a segment.

        Args:
            task (str): In [first, last, residue, resids, segids] depending on
                what information is desired.
            segid (str): Segment ID to query
            resid (str): Residue ID to query, if task is "residue"

        Returns:
            (str or list of str): Requested information
        """
        return _psfgen.segment(psfstate=self._data, task=task,
                               segid=segid, resid=resid)

    #===========================================================================


