#!/usr/bin/python
import re
from Bio.PDB import *


def remove_nonRNA_atoms_from_model(model):
    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if str(residue.resname).strip() not in ['A', 'C', 'G', 'U']:
                    residue_to_remove.append((chain.id, residue.id))
            if len(chain) == 0:
                chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return model


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
