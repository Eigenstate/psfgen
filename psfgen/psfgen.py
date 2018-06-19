
from psfgen.psfgen_core import *
import _psfgen

class PsfgenFormatError(IOError):
    pass

# Definitions for psf file types
CHARMM="charmm"
XPLOR="x-plor"

# Definitions for auto arguments
AUTO_ANGLES, AUTO_DIHEDRALS, AUTO_NONE = "angles", "dihedrals", "none"

class Psfgen:
    def __init__(self):
        self.data = _psfgen.init_mol()

    def __del__(self):
        _psfgen.del_mol(self.data);

    def readCharmmTopology(self, filename):
        """
        Parses a charmm format topology file into current library.

        Args:
            filename (str): File to parse

        Raises:
            FileNotFoundError: If the file cannot be opened for reading
            ValueError: If an error occurs during parsing
        """
        _psfgen.parse_topology(psfstate=self.data, filename=filename)

    def aliasResidue(self, topResName, pdbResName):
        _psfgen.alias(psfstate=self.data, type="residue",
                      name=topResName, newname=pdbResName)

    def aliasAtom(self, resName, topAtomName, pdbAtomName):
        _psfgen.alias(psfstate=self.data, type="atom",
                      resname=resName, name=topAtomName, newname=pdbAtomName)

    def addSegment(self, segid, first=None, last=None, auto=None,
                   pdb=None, residues=None, mutate=None):
        if segid in _psfgen.segment(psfstate=self.data, task="segids"):
            raise ValueError("Duplicate segID")

        # Handle first
        if first is not None:
            if topo_mol_segment_first(mol=self.mol, rname=first):
                raise ValueError("Failed to set patch for first residue")

        # Handle last
        if last is not None:
            if topo_mol_segment_last(mol=self.mol, rname=last):
                raise ValueError("Failed to set patch for last residue")

        # Handle auto
        if auto is not None:
            angles, dihedrals = [0,0]

            if isinstance(auto, str):
                auto = [auto]
            for arg in auto:
                if arg is AUTO_ANGLES:
                    angles = 1
                if arg is AUTO_DIHEDRALS:
                    dihedrals = 1
                if arg is AUTO_NONE:
                    angles=0
                    dihedrals=0
            if topo_mol_segment_auto_angles(mol=self.mol, autogen=angles):
                raise ValueError("Failed setting angle autogen")
            if topo_mol_segment_auto_dihedrals(mol=self.mol, autogen=dihedrals):
                raise ValueError("Failed setting dihedral autogen")

        # Handle pdb
        if isinstance(pdb, str):
            pdb = [pdb]
        for file in pdb:
            fd = fopen(file, 'r')
            retval = pdb_file_extract_residues(mol=self.mol, file=fd,
                                               h=self.aliases, all_caps=1,
                                               arg5=None, print_msg=None)
            fclose(fd)
            if retval:
                raise ValueError("Error reading residues from pdb file %s"
                                 % file)

        # Handle residues
        if residues is not None:
            print("doing residues")
            for resid, resname in residues.items():
                # TODO: chain argument not handled here
                if topo_mol_residue(mol=self.mol, resid=str(resid),
                                    resname=resname, chain=None):
                    raise ValueError("failed on residue %s:%s"
                                     % (resid,resname))

        # Handle mutate
        if mutate is not None:
            for resid, resname in mutate.items():
                if topo_mol_mutate(mol=self.mol, resid=str(resid),
                                   rname=resname):
                    raise ValueError("failed on mutate %s:%s" % (resid,resname))

        if topo_mol_end(self.mol):
            raise ValueError("failed on end of segment %s" % segid)

    def readCoords(self, pdbfile, segid):
        if segid not in _psfgen.segment(psfstate=self.data, task="segids"):
            raise ValueError("Can't read coordinates for segment '%s' as "
                             "it is undefined." % segid)
        _psfgen.read_coords(psfstate=self.data,
                            filename=pdbfile,
                            segid=segid)

    def setCoords(self, segid, resid, atomname, position):
        """
        Sets the coordinates of a given atom to new values.

        Args:
            segid (str): Segment ID of atom
            resid (str): Residue ID of atom
            atomname (str): Atom name
            position (3-tuple of double): New x, y, and z coordinates for atom

        Returns:
            True if successful
        """
        return not _psefgen.set_coords(mol=self.mol, segid=segid, resid=resid,
                                       aname=atomname, position=position)

    def guessCoords(self):
        if topo_mol_guess_xyz(self.mol):
            raise ValueError("failed on guessing coordinates")

    def patch(self, patchName, targets):
        identlist = [{"segid": segid.encode(),
                      "resid": str(resid).encode()}
                     for segid, resid in targets]
        if topo_mol_patch(mol=self.mol, targets=identlist, rname=patchName,
                          prepend=False, warn_angles=False,
                          warn_dihedrals=False, deflt=0):
            raise ValueError("failed to apply patch %s to %s"
                             % (patchName, targets))

    def multiply(self, ncopies, groups):
        identlist = []
        for group in groups:
            if len(group) == 2:
                identlist.append({"segid": str(group[0]),
                                  "resid": str(group[1])})
            elif len(group) == 3:
                identlist.append({"segid": str(group[0]),
                                  "resid": str(group[1]),
                                  "aname":str(group[2])})
            else:
                raise ValueError("groups must contain lists of 2 or 3 elements")
        rc = topo_mol_multiply_atoms(mol=self.mol, targets=identlist,
                                     ncopies=ncopies)
        if rc:
            raise ValueError("failed to multiply atoms (error=%d)" % rc)

    def delAtom(self, segid, resid=None, atomName=None):
        ident={"segid": str(segid).encode()}
        if resID is not None:
            ident["resid"] = str(resid).encode()
        if atomName is not None:
            ident["aname"] = str(atomName).encode()
        topo_mol_delete_atom(mol=self.mol, target=ident)

    def regenerateAngles(self):
        topo_mol_regenerate_angles(self.mol)

    def regenerateDihedrals(self):
        topo_mol_regenerate_dihedrals(self.mol)

    def readPSF(self, filename):
        fd = fopen(filename, 'r')
        retval = psf_file_extract(self.mol, fd, None, None)
        fclose(fd)
        if retval:
            raise PsfgenFormatError("Error reading psf file '%s'"
                                    % filename)

    def writePSF(self, filename, type=XPLOR):
        if not type in [XPLOR, CHARMM]:
            raise ValueError("type must be either '%s' or '%s'"
                             % (CHARMM, XPLOR))
        fd = fopen(filename, 'w')
        ischarmm = 0
        if type == CHARMM:
            ischarmm = 1
        topo_mol_write_psf(mol=self.mol, file=fd, charmmfmt=ischarmm,
                           nocmap=False, nopatches=False,
                           arg6=None, print_msg=None)
        fclose(fd)

    def writePDB(self, filename):
        fd = fopen(filename, 'w')
        retval = topo_mol_write_pdb(mol=self.mol, file=fd,
                                    arg3=None, print_msg=None)
        fclose(fd)
        if retval:
            raise IOError("failed on writing coordinates to pdb file %s"
                          % filename)

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
        return _psfgen.segment(psfstate=self.data, task=task,
                               segid=segid, resid=resid)


