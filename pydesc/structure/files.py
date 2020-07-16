class PDBWriter:
    @staticmethod
    def create_pdb_string(substructure, transformed=True):
        """Returns string formatted as PDB file representing given (sub)structure.

        Arguments:
        substructure -- instance of pydesc AbstractStructure.
        transformed -- initially set to True, if so - creates PyMOL object
        with respect for all previous movements; otherwise uses coordinates
        from pdb file.
        """
        line_n = 0
        lines = []
        for mer in substructure:
            pdb_id = substructure.converter.get_pdb_id(mer.ind)
            for atom_name, atom in mer.atoms.items():
                if transformed:
                    coord = atom.get_coord(substructure.trt_matrix)
                else:
                    coord = atom.get_coord()
                icode = " " if pdb_id.icode is None else pdb_id.icode
                pdb_line = (
                    f"ATOM  {line_n:5d} {atom_name:>4s} {mer.name:3s}"
                    f"{pdb_id.chain:>2s}{pdb_id.ind:4d}{icode:1s}"
                    f"   {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                    f"{atom.occupancy:6.2f} {atom.b_factor:5.2f}         "
                    f"{atom.element:2s}"
                )
                line_n += 1
                lines.append(pdb_line)
        lines.append("END")
        return "\n".join(lines)
