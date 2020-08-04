#!/usr/bin/python
from Commons.AngularCalculator import calculate_euler_angles_pairwise
from Commons.DataStructures import Junction, Stem, Connector
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
            current_junction.list_of_stems = []
            current_stem = Stem()
            current_connector = Connector()

            bool_found, number_of_stems, position_of_connectors = find_junction(text, start)

            if bool_found:
                list_of_euler_angles, planar_angle, junction_substructure_file= calculate_euler_angles_pairwise(position_of_connectors, name_of_structure, chains, method)
                current_junction.type = str(number_of_stems) + '-way junction'
                current_junction.name_of_file = junction_substructure_file

                current_stem.segment_length = 1
                current_stem.list_of_connectors = []

                current_connector.list_of_segments_ranges = []
                current_connector.list_of_segments_ranges = generate_fragments(position_of_connectors)
                current_connector.lengths_of_segments = []
                current_connector.list_of_segment_db = []
                current_connector.list_of_segment_seq = []
                current_connector.list_of_angles = []
                current_connector.planar_angle = []


                for i,pair in enumerate(current_connector.list_of_segments_ranges):
                    current_connector.lengths_of_segments.append(pair[1] -1 - pair[0])
                    current_connector.list_of_segment_db.append(text[pair[0] - 1:pair[1]])
                    current_connector.list_of_segment_seq.append(sequence[pair[0] - 1:pair[1]])
                    current_connector.list_of_angles.append(list_of_euler_angles[i])
                    current_connector.planar_angle.append(planar_angle[i])

                current_stem.list_of_connectors.append(current_connector)
                current_junction.list_of_stems.append(current_stem)
                list_of_junctions.append(current_junction)
                types_of_junction.append(number_of_stems)

        return types_of_junction, list_of_junctions
