#!/usr/bin/python
import vg
import numpy as np
from Bio.PDB import *
from Bio.PDB import PDBList
from pathlib import Path
from Bio.PDB.mmcifio import MMCIFIO

from Commons.TetriaryStructuresTools import StructureSelection, points_to_vector_converter, generate_fragments
from Commons.Utilities import is_non_zero_file, pairs_generator


def centroid_calculator(lis):
    length = len(lis)
    sum1 = map(sum, zip(*lis))
#    return sum1[0]/length, sum1[1]/length, sum1[2]/length
    return [val/length for val in sum1]

def planar_angle_calculator(first_stem,second_stem):
    return vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]))

def euler_angle_calculator(first_stem,second_stem):
    angle_x = vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.x)
    angle_y = vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.y)
    angle_z = vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.z)
    return [format(angle_x, '.3f'), format(angle_y, '.3f'), format(angle_z, '.3f')]


def calculate_euler_angles_pairwise(list_of_stem_pairs, structure_name, structure_chains, method_used):
    print(structure_name)
    io = PDBIO()
    list_of_points = []
    list_of_euler_angles = []
    list_of_planar_angles = []
    list_of_residues = []
    pdb_data_folder = Path("PDB_files/")
    structure_name = structure_name.lower()
    # Download file if does not exist
    if (not is_non_zero_file(pdb_data_folder / (structure_name + ".pdb"))) and (not is_non_zero_file(pdb_data_folder / (structure_name + ".cif"))):
        PDBList().retrieve_pdb_file(structure_name,pdir=pdb_data_folder)

    if is_non_zero_file(pdb_data_folder / (structure_name + ".pdb")):
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure(structure_name, pdb_data_folder / (structure_name + ".pdb"))
    if is_non_zero_file(pdb_data_folder / (structure_name + ".cif")):
        try:
            parser_cif = MMCIFParser(QUIET=True)
            structure = parser_cif.get_structure(structure_name, pdb_data_folder / (structure_name + ".cif"))
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

    for pair in list_of_stem_pairs:
        sum_of_chain_lenght = len(chain.child_list)
        chain_id = 0
        while pair[0]>sum_of_chain_lenght:
            chain_id += 1
            sum_of_chain_lenght += len(chain_test[chain_id].child_list)


        residue_1 = chain_test[chain_id].child_list[pair[0]-1-(sum_of_chain_lenght - len(chain_test[chain_id].child_list))]


        sum_of_chain_lenght = len(chain.child_list)
        chain_id = 0
        while pair[1]>sum_of_chain_lenght:
            chain_id += 1
            sum_of_chain_lenght += len(chain_test[chain_id].child_list)


        residue_2 = chain_test[chain_id].child_list[pair[1]-1-(sum_of_chain_lenght - len(chain_test[chain_id].child_list))]


        list_of_residues.append([residue_1.id[1], residue_2.id[1]])

        atoms_to_include = []
        atoms_to_include += (Selection.unfold_entities(residue_1, 'A'))
        atoms_to_include += (Selection.unfold_entities(residue_2, 'A'))
        list_of_atoms = []
        for atom in atoms_to_include:
            list_of_atoms.append([atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]])
        center_of_pair = centroid_calculator(list_of_atoms)
        list_of_points.append(center_of_pair)

    center_of_junction = centroid_calculator(list_of_points)


    name_of_file = save_structure(structure, list_of_stem_pairs, structure_name, list_of_residues)

    for pair in pairs_generator(range(len(list_of_points))):
        first_stem = points_to_vector_converter(center_of_junction, list_of_points[pair[0]])
        second_stem = points_to_vector_converter(center_of_junction, list_of_points[pair[1]])

        planar_angle = planar_angle_calculator(first_stem,second_stem)
        list_of_planar_angles.append(format(planar_angle, '.3f'))

        list_of_euler_angles.append(euler_angle_calculator(first_stem,second_stem))

    return list_of_euler_angles, list_of_planar_angles, name_of_file, True

def save_structure(structure, list_of_stem_pairs, structure_name, list_of_residues):
    list_of_fragments = generate_fragments(list_of_stem_pairs)
    stems_location = '_'.join(str(item) for sublist in list_of_fragments for item in sublist)
    name_of_file = structure_name + '_' + str(len(list_of_stem_pairs)) + '-way_junction' + '_' + stems_location + '.cif'
    if not is_non_zero_file('./output/structures/' + (structure_name + ".cif")):
        io=MMCIFIO()
        io.set_structure(structure)

        name_of_file = structure_name + '_' + str(len(list_of_stem_pairs)) + '-way_junction' + '_' + stems_location + '.cif'
        io.save('./output/structures/' + name_of_file, StructureSelection(list_of_residues))
    return name_of_file
