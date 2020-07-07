from unittest import TestCase
import numpy as np
from pathlib import Path
from Bio.PDB import *
from shutil import rmtree
import io
import os

from Commons.AngularCalculator import planar_angle_calculator
from Commons.AngularCalculator import centroid_calculator, save_structure


class TestAngularCalculator(TestCase):
    def test_euler_angle_calculator(self):

        from Commons.AngularCalculator import euler_angle_calculator
        first_stem = [0,0,1]
        second_stem =[0,1,0]

        self.assertEqual(euler_angle_calculator(first_stem,second_stem)[0], 90)
        self.assertTrue(np.math.isnan(euler_angle_calculator(first_stem, second_stem)[1]))
        self.assertTrue(np.math.isnan(euler_angle_calculator(first_stem, second_stem)[2]))

    def test_euler_angle_calculator_no_change(self):

        from Commons.AngularCalculator import euler_angle_calculator
        first_stem = [0,0,0]
        second_stem = [0,0,0]

        self.assertTrue(np.math.isnan(euler_angle_calculator(first_stem, second_stem)[0]))
        self.assertTrue(np.math.isnan(euler_angle_calculator(first_stem, second_stem)[1]))
        self.assertTrue(np.math.isnan(euler_angle_calculator(first_stem, second_stem)[2]))

    def test_planar_angle_calculator(self):


        first_stem = [0,0,1]
        second_stem =[0,1,0]

        self.assertEqual(planar_angle_calculator(first_stem,second_stem), 90.0)

    def test_planar_angle_calculator_no_change(self):

        first_stem = [6,6,6]
        second_stem =[1,1,1]

        self.assertEqual(planar_angle_calculator(first_stem,second_stem), 0.0)


    def test_centroid_calculator_cube(self):

        #defines points for cube
        set_of_points_cube = [[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]]
        self.assertEqual(centroid_calculator(set_of_points_cube), [0.5, 0.5, 0.5])


    def test_centroid_calculator_triangle(self):
        #triangle on the same plane
        set_of_points_rectangle = [[0,0,0], [0,3,0], [3,0,0]]
        self.assertEqual(centroid_calculator(set_of_points_rectangle), [1, 1, 0])

    def test_save_structure(self):
        if not os.path.exists(Path('./output/structures')):
            os.makedirs(Path('./output/structures'))
        parser_cif = MMCIFParser()
        pdb_data_folder = Path('../test_files/')
        structure = parser_cif.get_structure('1b23', pdb_data_folder / '1b23.cif')
        list_of_residues = [[7, 66], [10, 25], [27, 43], [49, 65]]
        structure_name = '1b23'
        list_of_stem_pairs= [[7, 64], [10, 24], [26, 42], [47, 63]]
        save_structure(structure, list_of_stem_pairs, structure_name, list_of_residues)
        tst_path = Path('./output/structures/1b23_4-way_junction_7_10_24_26_42_47_63_64.pdb')
        ref_path = Path('../test_files/1b23_4-way_junction_7_10_24_26_42_47_63_64.pdb')
        self.assertListEqual(
            list(io.open(tst_path)),
            list(io.open(ref_path)))
        rmtree(Path("./output"))

