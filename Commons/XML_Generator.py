#!/usr/bin/python
import xml.etree.cElementTree as ET


def xml_generate(list_of_records):

    root = ET.Element("File")
    for record in list_of_records:
        pdb_structure = ET.SubElement(root, "PDB_Structure")
        ET.SubElement(pdb_structure, "field1", name="PDB_ID").text = record.name

        for junction_record in record.list_of_junctions:
            junction = ET.SubElement(pdb_structure, "Junction")
            ET.SubElement(junction, "field1", name="Type_of_junction").text = junction_record.type
            ET.SubElement(junction, "field2", name="Junction_PDB_File").text = junction_record.name_of_file

            for connector_id, (list_of_segments_ranges, lenghts_of_segments, list_of_segment_seq,  list_of_segment_db, list_of_angles, planar_angle) \
                    in enumerate(zip(junction_record.list_of_segments_ranges, junction_record.lengths_of_segments, junction_record.list_of_segment_seq,  junction_record.list_of_segment_db, junction_record.list_of_angles, junction_record.planar_angle)):
                connectors = ET.SubElement(junction, "Connectors")

                connector = ET.SubElement(connectors, "Connector")
                ET.SubElement(connector, "field1", name="ID").text = str(connector_id)
                ET.SubElement(connector, "field2", name="Range").text = str(list_of_segments_ranges)
                ET.SubElement(connector, "field3", name="Lenght").text = str(lenghts_of_segments)
                ET.SubElement(connector, "field4", name="Sequence").text = str(list_of_segment_seq)
                ET.SubElement(connector, "field5", name="DB_Notation").text = str(list_of_segment_db)
                ET.SubElement(connector, "field6", name="Planar_angle").text = str(planar_angle)

                euler_angles = ET.SubElement(connector, "Euler_Angles")
                ET.SubElement(euler_angles, "field1", name="Angle_X").text = str(list_of_angles[0])
                ET.SubElement(euler_angles, "field2", name="Angle_Y").text = str(list_of_angles[1])
                ET.SubElement(euler_angles, "field3", name="Angle_Z").text = str(list_of_angles[2])

    tree = ET.ElementTree(root)

    tree.write("./output/RESULTS.xml")
