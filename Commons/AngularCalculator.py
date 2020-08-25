#!/usr/bin/python
import vg
import numpy as np
from Bio.PDB import *
from Bio.PDB import PDBList
from pathlib import Path
from Bio.PDB.mmcifio import MMCIFIO

from Commons.TetriaryStructuresTools import StructureSelection2, points_to_vector_converter, generate_fragments
from Commons.Utilities import is_non_zero_file, pairs_generator


def centroid_calculator(list_of_atoms):
    """
    Function to find center for a set of atoms
    :param list_of_atoms: list of atom coordinates
    :return: center point
    """
    length = len(list_of_atoms)
    sum1 = map(sum, zip(*list_of_atoms))
    return [val/length for val in sum1]

def planar_angle_calculator(first_stem,second_stem):
    """
    Function to calculate planar angle between 2 vectors
    :param first_stem: stem coordinates
    :param second_stem: stem coordinates
    :return: value of planar angle between two vectors
    """
    return np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]])))

def euler_angle_calculator(first_stem,second_stem):

    """
    Function to calculate Euler angles between 2 vectors
    :param first_stem: stem coordinates
    :param second_stem: stem coordinates
    :return: value of euler angles between two vectors
    """
    angle_x = np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.x))
    angle_y = np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.y))
    angle_z = np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.z))
    return [format(angle_x, '.3f'), format(angle_y, '.3f'), format(angle_z, '.3f')]

def calculate_euler_angles_pairwise(list_of_stem_pairs, structure_name, structure_chains, save_substructure = True):
    print(structure_name)

    list_of_points = []
    list_of_euler_angles = []
    list_of_planar_angles = []
    list_of_residues = []
    base_list_of_residues=[]
    base_list_of_points = []

    pdb_data_folder = Path("PDB_files/")
    structure_name = structure_name.lower()

    try:
        structure = get_structure(pdb_data_folder, structure_name)
    except:
        return [], [], '', False

    model = structure[0]
    chain_test = []
    chain = None

    for chain_to_include in structure_chains:
        if not chain:
            chain = model[chain_to_include]
            chain_test.append(model[chain_to_include])
        else:
            chain_test.append(model[chain_to_include])

    for pair in list_of_stem_pairs[1]:
        residue_1 = get_residue(chain_test, pair[0])
        residue_2 = get_residue(chain_test, pair[1])

        list_of_residues = [*list_of_residues, [residue_1.full_id, residue_2.full_id]]

        atoms_to_include = [*Selection.unfold_entities(residue_1, 'A'), *Selection.unfold_entities(residue_2, 'A')]
        list_of_atoms = [[atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]] for atom in atoms_to_include]

        list_of_points = [*list_of_points, centroid_calculator(list_of_atoms)]

    for base_pair in list_of_stem_pairs[0]:
        base_residue_1 = get_residue(chain_test, base_pair[0])
        base_residue_2 = get_residue(chain_test, base_pair[1])

        base_list_of_residues = [*base_list_of_residues, [base_residue_1.full_id, base_residue_2.full_id]]

        base_atoms_to_include = [*Selection.unfold_entities(base_residue_1, 'A'), *Selection.unfold_entities(base_residue_2, 'A')]
        base_list_of_atoms = [[atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]] for atom in base_atoms_to_include]

        base_list_of_points = [*base_list_of_points, centroid_calculator(base_list_of_atoms)]
    center_of_junction = centroid_calculator(base_list_of_points)

    #substructure save
    if save_substructure:
        all_included_residues = []
        pairs = generate_fragments(list_of_residues)
        for pair in pairs:
            if pair[0][2] == pair[1][2]:
                for x in range(pair[0][3][1], pair[1][3][1]+1):
                    all_included_residues.append(get_residue(chain_test, x, True))
            else:
                for x in range(pair[0][3][1], len(model[pair[0][2]])):
                    all_included_residues.append(get_residue(chain_test, x, True))
                for x in range(1, pair[1][3][1]):
                    all_included_residues.append(get_residue(chain_test, x, True))
        print(all_included_residues)
        name_of_file = save_structure(structure, list_of_stem_pairs[1], structure_name, all_included_residues)

    else:
        name_of_file = ''

    #calculation of euler and planar angles
    for pair in pairs_generator(range(len(list_of_points))):
        first_stem = points_to_vector_converter(center_of_junction, list_of_points[pair[0]])
        second_stem = points_to_vector_converter(center_of_junction, list_of_points[pair[1]])

        planar_angle = planar_angle_calculator(first_stem,second_stem)
        list_of_planar_angles.append(format(planar_angle, '.3f'))

        list_of_euler_angles.append(euler_angle_calculator(first_stem,second_stem))

    return list_of_euler_angles, list_of_planar_angles, name_of_file, True

def save_structure(structure, list_of_stem_pairs, structure_PDB_ID, list_of_residues):
    """
    Function creates the file with junction 3D model (*.cif format) and returns its name
    :param structure: the whole PDB structure
    :param list_of_stem_pairs: list of all pairs that junction consist of
    :param structure_PDB_ID: PDB ID of whole structure
    :param list_of_residues: list of residues full IDs that junction consist of
    :return: name of the created file
    """
    list_of_fragments = generate_fragments(list_of_stem_pairs)
    stems_location = '_'.join(str(item) for sublist in list_of_fragments for item in sublist)
    name_of_file = structure_PDB_ID + '_' + str(len(list_of_stem_pairs)) + '-way_junction' + '_' + stems_location + '.cif'
    if not is_non_zero_file('./output/structures/' + (structure_PDB_ID + ".cif")):
        io=MMCIFIO()
        io.set_structure(structure)
        name_of_file = structure_PDB_ID + '_' + str(len(list_of_stem_pairs)) + '-way_junction' + '_' + stems_location + '.cif'
        io.save('./output/structures/' + name_of_file, StructureSelection2(list_of_residues))
    return name_of_file

def get_structure(pdb_data_folder, structure_PDB_ID):
    """
    Function to retrieve information about structure of specified molecule
    :param pdb_data_folder: path to folder contaiting all 3D structures
    :param structure_PDB_ID: structures PDB ID
    :return: the structure object
    """
    if (not is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".pdb"))) and (not is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".cif"))):
        PDBList().retrieve_pdb_file(structure_PDB_ID, pdir=pdb_data_folder)

    if is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".pdb")):
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".pdb"))
    if is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".cif")):
        try:
            parser_cif = MMCIFParser(QUIET=True)
            structure = parser_cif.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".cif"))
        except:
            PDBList().retrieve_pdb_file(structure_PDB_ID, pdir=pdb_data_folder, file_format='pdb')
            parser_pdb = PDBParser()
            structure = parser_pdb.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".pdb"))
    return structure

def get_residue(list_of_chains, res_num, get_full_id = False):
    """
    Function to retrieve the resudie object or residue full id
    :param list_of_chains: list of all avaiable chains
    :param res_num: number of residue to retrieve
    :param get_full_id: bolean if True the function returns residue_full_id, else returns residue object
    :return: residue_full_id/residue object
    """
    sum_of_chain_length = len(list_of_chains[0].child_list)
    chain_id = 0
    while res_num>sum_of_chain_length:
        chain_id += 1
        sum_of_chain_length += len(list_of_chains[chain_id].child_list)
    if get_full_id:
        return list_of_chains[chain_id].child_list[res_num-1-(sum_of_chain_length - len(list_of_chains[chain_id].child_list))].full_id
    else:
        return list_of_chains[chain_id].child_list[res_num - 1 - (sum_of_chain_length - len(list_of_chains[chain_id].child_list))]


