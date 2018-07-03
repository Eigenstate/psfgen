
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

    def __init__(self, output=sys.stdout, case_sensitive=False):
        """
        Creates a PsfGen object.

        Args:
            outstream (str or stream object): Where to write output. Defaults
                to stdout. Must have fileno() attribute, or if a string, file
                will be opened.
            case_sensitive (bool): Whether or not residue names and definitions
                are considered to be case sensitive
        """
        if isinstance(output, str):
            self.output = open(output, 'wb')
        elif hasattr(output, "fileno"):
            self.output = output
        else:
            raise ValueError("output argument must be a str or open file")

        self._fileno = self.output.fileno()
        self._data = _psfgen.init_mol(outfd=self._fileno)

        self._read_topos = False # Cannot change case sensitivity if true
        self._allcaps = case_sensitive
        _psfgen.set_allcaps(psfstate=self._data, allcaps=not case_sensitive)

    #===========================================================================

    def __del__(self):

        if not self.output.closed and self.output is not sys.stdout:
            self.output.close()

        _psfgen.del_mol(self._data);

    #===========================================================================

    # This property decorator lets the case sensitivity be a boolean attribute
    # of the PsfGen instance that calls the C code to update the internal
    # psfgen_data* object when set.

    @property
    def case_sensitive(self):
        """
        Determines whether or not residue and other names have case sensitivity.
        Defaults to False. Can only be changed before reading in any topology
        files.
        """
        return self._allcaps

    @case_sensitive.setter
    def case_sensitive(self, value):
        if not isinstance(value, bool):
            raise ValueError("case_sensitive must be a boolean value")

        if self._read_topos:
            raise ValueError("Cannot change case sensitivity value after "
                             "reading in topology files")

        _psfgen.set_allcaps(psfstate=self._data, allcaps=not value)
        self._allcaps = value

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
        self._read_topos = True

    #===========================================================================

    def alias_residue(self, top_resname, pdb_resname):
        """
        Set a correspondence between a residue name in a PDB file and in the
        defined topologies. The topology residue names should be treated as the
        canonical set here.

        Args:
            top_resname (str): Resname in topology file
            pdb_resname (str): Equivalent resname in PDB files
        """
        _psfgen.alias(psfstate=self._data, type="residue",
                      newname=top_resname, name=pdb_resname)

    #===========================================================================

    def alias_atom(self, resname, top_atomname, pdb_atomname):
        """
        Set a correspondence between an atom name in a specific residue in a
        PDB file and in the defined topologies.

        Args:
            resname (str): Residue name in which to make alias
            top_atomname (str): Atom name in the topology file
            pdb_atomname (str): Equivalent atom name in PDB file
        """
        _psfgen.alias(psfstate=self._data, type="atom",
                      resname=resname, newname=top_atomname, name=pdb_atomname)

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

    def get_atom_names(self, segid, resid):
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

        return _psfgen.query_atoms(psfstate=self._data, segid=segid,
                                   resid=resid, task="name")

    #===========================================================================

    def get_masses(self, segid, resid):
        """
        Obtains atom masses in a given residue

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (list of float): Atom masses in residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.query_atoms(psfstate=self._data, segid=segid,
                                   resid=resid, task="mass")

    #===========================================================================

    def get_charges(self, segid, resid):
        """
        Obtains atom charges in a given residue

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (list of float): Atom charges in residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.query_atoms(psfstate=self._data, segid=segid,
                                   resid=resid, task="charge")

    #===========================================================================

    def get_atom_indices(self, segid, resid):
        """
        Obtains atom indices/IDs in a given residue

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (list of int): Atom IDs in residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.query_atoms(psfstate=self._data, segid=segid,
                                   resid=resid, task="atomid")

    #===========================================================================

    def get_coordinates(self, segid, resid):
        """
        Obtains atom coordinates in a given residue

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (list of 3-tuple): (x,y,z) position of all atoms in residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.query_atoms(psfstate=self._data, segid=segid,
                                   resid=resid, task="coordinates")

    #===========================================================================

    def get_velocities(self, segid, resid):
        """
        Obtains atom velocities in a given residue, if set

        Args:
            segid (str): Segment ID to query
            resid (str or int): Residue ID to query

        Returns:
            (list of 3-tuple): (vx,vy,vz) velocities of all atoms in residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        return _psfgen.query_atoms(psfstate=self._data, segid=segid,
                                   resid=resid, task="velocities")

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
        """
        _psfgen.set_coord(psfstate=self._data, segid=segid,
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

    def delete_atoms(self, segid, resid=None, atomname=None):
        """
        Deletes atoms from the molecule, with options for increasing
        specificity. As many atoms as match the passed options will be
        deleted, so for example if only segid is present, an entire segment
        will be removed.

        Args:
            segid (str): Segment ID of atom(s) to delete
            resid (str or int): Residue ID of atom(s) to delete, or None to
                delete entire matching segment
            atomname (str): Atom name to delete, or None to delete entire
                matching segment or residue
        """
        if isinstance(resid, int):
            resid = str(resid)

        _psfgen.delete_atoms(psfstate=self._data, segid=segid, resid=resid,
                             aname=atomname)

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

    def read_psf(self, filename, pdbfile=None, namdbinfile=None,
                 velnamdbinfile=None):
        """
        Reads structure information from a PSF file and adds it to the internal
        molecule state. Can also read insertion codes and coordinates from a
        PDB file, and/or coordinates from a NAMD binary file, if atoms are in
        the same order. Can read velocities from a NAMD binary file as well.

        The PSF file is read in and also kind of interpreted: topology files
        in the REMARKS section are loaded, and a new segment is created with
        the segid specified in the PSF file. Note that if the topology files
        cannot be found, they will still be listed as "loaded" in the internal
        PsfGen state.

        Args:
            filename (str): PSF file to read
            pdbfile (str): PDB file to obtain coordinates, elements, and
                insertion codes from. Atoms must be in same order as the PSF.
            namdbinfile (str): Binary NAMD file to read coordinates from. Will
                take priority over coordinates in pdbfile, if specified. Atoms
                must be in the same order as the PSF.
            velnamdbinfile (str): Binary NAMD file to read velocities from.
                Atoms must be in the same order as the PSF.
        """
        _psfgen.read_psf(psfstate=self._data, filename=filename,
                         pdbfile=pdbfile, namdbinfile=namdbinfile,
                         velnamdbinfile=velnamdbinfile)

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

    def write_namdbin(self, filename, velocity_filename=None):
        """
        Writes the current molecule state as a NAMD binary file. This file
        format can only contain coordinate information for each atom. Optionally
        write velocities to another NAMD binary file.

        Args:
            filename (str): Filename to write to
            velocity_filename (str): If present, filename to write velocities to
        """
        _psfgen.write_namdbin(psfstate=self._data, filename=filename,
                              velocity_filename=velocity_filename)

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


