#!/usr/bin/python
from Commons.AngularCalculator import calculate_euler_angles_pairwise
from Commons.DataStructures import Junction
from Commons.SecondaryStructureTools import find_junction
from Commons.TetriaryStructuresTools import generate_fragments


class JunctionFinder:
    def __init__(self, text, sequence, name_of_structure, chains, method):
        self.text = text
        self.sequence = sequence
        self.name_of_structure = name_of_structure
        self.chains = chains
        self.method = method

    @staticmethod
    def find_junctions(text, sequence, name_of_structure, chains, method):

        types_of_junction = []
        list_of_junctions = []

        for start in range(len(text)):
            current_junction = Junction()
            bool_found, number_of_stems, position_of_connectors = find_junction(text, start)

            if bool_found:
                list_of_euler_angles, planar_angle, junction_substructure_file= calculate_euler_angles_pairwise(position_of_connectors, name_of_structure, chains, method)
                current_junction.type = str(number_of_stems) + '-way junction'
                current_junction.list_of_segments_ranges = []
                current_junction.list_of_segments_ranges = generate_fragments(position_of_connectors)
                current_junction.lengths_of_segments = []
                current_junction.list_of_segment_db = []
                current_junction.list_of_segment_seq = []
                current_junction.list_of_angles = []
                current_junction.planar_angle = []
                current_junction.name_of_file = junction_substructure_file

                for i,pair in enumerate(current_junction.list_of_segments_ranges):
                    current_junction.lengths_of_segments.append(pair[1] - pair[0] - 1)
                    current_junction.list_of_segment_db.append(text[pair[0]:pair[1] + 1])
                    current_junction.list_of_segment_seq.append(sequence[pair[0]:pair[1] + 1])
                    current_junction.list_of_angles.append(list_of_euler_angles[i])
                    current_junction.planar_angle.append(planar_angle[i])

                list_of_junctions.append(current_junction)
                types_of_junction.append(number_of_stems)

        return types_of_junction, list_of_junctions
