import prody as pr


def read_pdb(pdb_file):
    """Reads a PDB file and returns a ProDy atom group."""
    return pr.parsePDB(pdb_file)



