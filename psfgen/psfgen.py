
from psfgen.psfgen_core import *

class PsfgenFormatError(IOError):
    pass

# Definitions for psf file types
CHARMM="charmm"
XPLOR="x-plor"

# Definitions for auto arguments
AUTO_ANGLES, AUTO_DIHEDRALS, AUTO_NONE = "angles", "dihedrals", "none"

class Psfgen:
    def __init__(self):
        self.defs = topo_defs_create()
        self.aliases = stringhash_create()
        self.mol = topo_mol_create(self.defs)

    def __del__(self):
        topo_mol_destroy(self.mol)
        topo_defs_destroy(self.defs)
        stringhash_destroy(self.aliases)

    def readCharmmTopology(self, filename):
        with open(filename, 'r') as fd:
            if charmm_parse_topo_defs(self.defs, fd.fileno(), None, None):
                raise PsfgenFormatError("Error reading topology file")

    def aliasResidue(self, topResName, pdbResName):
        if extract_alias_residue_define(self.aliases,
                                        str(topResName),
                                        str(pdbResName)):
            raise ValueError("failed on residue alias")

    def aliasAtom(self, resName, topAtomName, pdbAtomName):
        if extract_alias_atom_define(self.aliases, resName,
                                     topAtomName, pdbAtomName):
            raise ValueError("failed on atom alias")

    def addSegment(self, segID, first=None, last=None, auto=None,
                   pdb=None, residues=None, mutate=None):
        if topo_mol_segment(self.mol, segID):
            raise ValueError("Duplicate segID")

        # Handle first
        if first is not None:
            if topo_mol_segment_first(self.mol, first):
                raise ValueError("Failed to set patch for first residue")

        # Handle last
        if last is not None:
            if topo_mol_segment_last(self.mol, last):
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
            if topo_mol_segment_auto_angles(self.mol, angles):
                raise ValueError("Failed setting angle autogen")
            if topo_mol_segment_auto_dihedrals(self.mol, dihedrals):
                raise ValueError("Failed setting dihedral autogen")

        # Handle pdb
        if isinstance(pdb, str):
            pdb = [pdb]
        for file in pdb:
            with  open(file, 'r') as fd:
                if pdb_file_extract_residues(self.mol, fd.fileno(), self.aliases,
                                             None, None):
                    raise ValueError("Error reading residues from pdb file %s"
                                     % file)

        # Handle residues
        if residues is not None:
            print("doing residues")
            for resid, resname in residues.items():
                if topo_mol_residue(self.mol, str(resid), str(resname)):
                    raise ValueError("failed on residue %s:%s"
                                     % (resid,resname))

        # Handle mutate
        if mutate is not None:
            for resid, resname in mutate.items():
                if topo_mol_mutate(self.mol, str(resid), str(resname)):
                    raise ValueError("failed on mutate %s:%s" % (resid,resname))

        if topo_mol_end(self.mol):
            raise ValueError("failed on end of segment %s" % segID)

    def readCoords(self, pdbfile, segID=None):
        with open(pdbfile, 'r') as fd:
            if pdb_file_extract_coordinates(self.mol, fd.fileno(), segID,
                                            self.aliases, None, None):
                raise ValueError("failed on reading coordinates from pdb file %s"
                                 % pdbfile)

    def setCoords(self, segID, resID, aName, pos):
        target = {"segid": segID, "resid": resID, "aname": aName }
        if topo_mol_set_xyz(self.mol, target, pos[0], pos[1], pos[2]):
            raise ValueError("failed on coord for segid %d resid %d"
                             % (segID, resID))

    def guessCoords(self):
        if topo_mol_guess_xyz(self.mol):
            raise ValueError("failed on guessing coordinates")

    def patch(self, patchName, targets):
        identlist = [{"segid": str(segid), "resid": str(resid)}
                     for segid, resid in targets]
        if topo_mol_patch(self.mol, identlist, patchName, 0,0,0):
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
        rc = topo_mol_multiply_atoms(self.mol, identlist, ncopies)
        if rc:
            raise ValueError("failed to multiply atoms (error=%d)" % rc)

    def delAtom(self, segID, resID=None, atomName=None):
        ident={"segid": str(segID)}
        if resID is not None:
            ident["resid"] = str(resID)
        if atomName is not None:
            ident["aname"] = str(atomName)
        topo_mol_delete_atom(self.mol, ident)

    def readPSF(self, filename):
        with open(filename, 'r') as fd:
            if psf_file_extract(self.mol, fd.fileno(), None, None):
                raise PsfgenFormatError("Error reading psf file '%s'"
                                        % filename)

    def writePSF(self, filename, type=XPLOR):
        if not type in [XPLOR, CHARMM]:
            raise ValueError("type must be either '%s' or '%s'"
                             % (CHARMM, XPLOR))
        with open(filename, 'w') as fd:
            ischarmm = 0
            if type == CHARMM:
                ischarmm = 1
            topo_mol_write_psf(self.mol, fd.fileno(), ischarmm, None, None)

    def writePDB(self, filename):
        with open(filename, 'w') as fd:
            if topo_mol_write_pdb(self.mol, fd.fileno(), None, None):
                raise IOError("failed on writing coordinates to pdb file %s"
                              % filename)

