#!/usr/bin/python
import re
from Bio.PDB import *


def remove_nonRNA_atoms_from_model(model):
    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if str(residue.resname).strip() not in ['A', 'C', 'G', 'U' ]:
                    residue_to_remove.append((chain.id, residue.id))
            if len(chain) == 0:
                chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return model

def remove_HOH_from_model(structure):
    list_of_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    residue_to_remove= []
    chain_to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if (residue.id[0] in ('W', 'H_SO4')) or residue.resname in list_of_aa:
                    residue_to_remove.append((chain.id, residue.id))

            if len(chain) == 0:
                chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return structure

def standardize_model(structure):
    list_of_nucleotides = ['A', 'C', 'G', 'U', '1MA', '2MG', '5MC', '5MU', '7MG', 'H2U', 'M2G', 'OMC', 'OMG', 'PSU', 'YG', '4SU', 'MIA', 'I', 'DA', 'DT',
                             'DC', 'DU', 'DI', 'DA', 'DG']
    non_standard_residues = {'MIA':'A', '1MA':'A', '2MG':'G', '5MC': 'C', '5MU': 'U', '7MG': 'G', 'H2U': 'U', 'M2G':'G', 'OMC':'C', 'OMG':'G', 'PSU':'U', 'YG':'G', '4SU':'U', 'DA': 'A', 'DT': 'T',
                             'DC': 'C', 'DU': 'U', 'DI': 'I', 'DA': 'A', 'DG:':'G'}
    acceptable_atoms = ['C2', 'C4', 'C6', 'C8', 'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6', 'C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'', 'O2\'', 'O3\'', 'O4\'', 'O5\'', 'OP1', 'OP2', 'P' ]
    residue_to_remove= []
    chain_to_remove = []
    atom_to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname not in list_of_nucleotides:
                    residue_to_remove.append((chain.id, residue.id))
                if residue.resname in non_standard_residues:
                    residue.resname = non_standard_residues[residue.resname]
                    residue_id = list(residue.id)
                    residue_id[0] = ' '
                    residue.id = tuple(residue_id)
                    for atom in residue:
                        if atom.id not in acceptable_atoms:
                            atom_to_remove.append((chain.id, residue.id, atom.id))
                            atom_id = list(atom.full_id[3])
                            atom_id[0] = ' '
                            atom_tuple = tuple(atom_id)
                            atom.full_id = [atom.full_id[0], atom.full_id[1], atom.full_id[2], atom_tuple, atom.full_id[4]]


            if len(chain) == 0:
                chain_to_remove.append(chain.id)
    for atom in atom_to_remove:
        model[atom[0]][atom[1]].detach_child(atom[2])

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return structure


def points_to_vector_converter(center_of_junction, atom):
    x = atom[0] - center_of_junction[0]
    y = atom[1] - center_of_junction[1]
    z = atom[2] - center_of_junction[2]
    vector = [x, y, z]
    return vector


def generate_fragments(list_of_pairs):
    fragments = []
    for elements in range(1, len(list_of_pairs)):
        if elements == 1:
            fragments.append([list_of_pairs[0][0], list_of_pairs[elements][0]])
        elif elements == len(list_of_pairs)-1:
            fragments.append([list_of_pairs[elements - 1][1], list_of_pairs[elements][0]])
            fragments.append([list_of_pairs[elements][1], list_of_pairs[0][1]])
        else:
            fragments.append([list_of_pairs[elements - 1][1], list_of_pairs[elements][0]])
    return fragments


class StructureSelection(Select):
    def __init__(self, list_of_residues):
        self.list_of_residues = list_of_residues

    def accept_residue(self, residue):
        fragment_list = generate_fragments(self.list_of_residues)
        resudue_number_to_include = []
        list_res = ['A', 'U', 'C', 'G']
        for pair in fragment_list:
            for numbers in range(pair[0], pair[1]+1):
                resudue_number_to_include.append(numbers)
        if (residue.id[1] in resudue_number_to_include) and (re.sub(r'\W+', '',residue.get_resname()) in list_res):
            return 1
        else:
            return 0

class StructureSelection2(Select):
    def __init__(self, list_of_residues):
        self.list_of_residues = list_of_residues

    def accept_residue(self, residue):

        if (residue.full_id in self.list_of_residues):
            return 1
        else:
            return 0
