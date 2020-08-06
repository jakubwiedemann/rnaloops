#!/usr/bin/python
from lxml import etree as ET
from Commons.Utilities import extract


def xml_generate(list_of_records):

    root = ET.Element("File")
    for record in list_of_records:
        pdb_structure = ET.SubElement(root, "PDB_Structure")
        ET.SubElement(pdb_structure, "field1", name="PDB_ID").text = record.name

        for junction_record in record.list_of_junctions:
            junction = ET.SubElement(pdb_structure, "Junction")
            ET.SubElement(junction, "field1", name="Type_of_junction").text = junction_record.type
            ET.SubElement(junction, "field2", name="Junction_PDB_File").text = junction_record.name_of_file

            stems = ET.SubElement(junction, "Stems")
            for stem_record in junction_record.list_of_stems:
                stem = ET.SubElement(stems, "Stem")
                ET.SubElement(stem, "field1", name="Stem length").text = str(stem_record.segment_length)
                ET.SubElement(stem, "field1", name="First strand sequence").text = str(stem_record.first_strand_sequence)
                ET.SubElement(stem, "field1", name="Second strand squence").text = str(stem_record.second_strand_sequence)
                ET.SubElement(stem, "field1", name="First strand DB_notation").text = str(stem_record.first_strand_db)
                ET.SubElement(stem, "field1", name="Second strand DB_notation").text = str(stem_record.second_strand_db)

                connectors = ET.SubElement(stem, "Connectors")
                for connector_records in stem_record.list_of_connectors:

                        connector = ET.SubElement(connectors, "Connector")
                        ET.SubElement(connector, "field1", name="ID").text = str(connector_records.connector_id)
                        ET.SubElement(connector, "field2", name="Range").text = str(connector_records.list_of_segments_ranges)
                        ET.SubElement(connector, "field3", name="Lenght").text = str(connector_records.lengths_of_segments)
                        ET.SubElement(connector, "field4", name="Sequence").text = str(connector_records.list_of_segment_seq)
                        ET.SubElement(connector, "field5", name="DB_Notation").text = str(connector_records.list_of_segment_db)
                        ET.SubElement(connector, "field6", name="Planar_angle").text = str(connector_records.planar_angle)

                        euler_angles = ET.SubElement(connector, "Euler_Angles")
                        ET.SubElement(euler_angles, "field1", name="Angle_X").text = str(extract(0, connector_records.list_of_angles))
                        ET.SubElement(euler_angles, "field2", name="Angle_Y").text = str(extract(1, connector_records.list_of_angles))
                        ET.SubElement(euler_angles, "field3", name="Angle_Z").text = str(extract(2, connector_records.list_of_angles))

    tree = ET.ElementTree(root)

    tree.write("./output/RESULTS.xml", pretty_print=True)
