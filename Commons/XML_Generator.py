#!/usr/bin/python
from lxml import etree as ET
from Commons.Utilities import extract


def xml_generate(list_of_records):

    separator = ', '
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
                ET.SubElement(stem, "field1", name="Stem length").text = str(separator.join(map(str,stem_record.segment_length)))
                ET.SubElement(stem, "field2", name="Duplex sequence").text = str(separator.join(map(str,stem_record.sequence)))
                ET.SubElement(stem, "field3", name="Connectors").text = str(separator.join(map(str,stem_record.list_of_connectors)))


            connectors = ET.SubElement(junction, "Connectors")
            for connector_records in junction_record.list_of_connectors:
                for connector_record in connector_records:
                    connector = ET.SubElement(connectors, "Connector")
                    ET.SubElement(connector, "field1", name="ID").text = str(connector_record.connector_id)
                    ET.SubElement(connector, "field2", name="Range").text = str(separator.join(map(str,connector_record.list_of_segments_ranges)))
                    ET.SubElement(connector, "field3", name="Lenght").text = str(connector_record.lengths_of_segments)
                    ET.SubElement(connector, "field4", name="Sequence").text = str(connector_record.list_of_segment_seq)
                    ET.SubElement(connector, "field5", name="DB_Notation").text = str(connector_record.list_of_segment_db)
                    ET.SubElement(connector, "field6", name="Planar_angle").text = str(separator.join(map(str,connector_record.planar_angle)))

                    euler_angles = ET.SubElement(connector, "Euler_Angles")
                    ET.SubElement(euler_angles, "field1", name="Angle_X").text = str(separator.join(map(str,extract(0, connector_record.list_of_angles))))
                    ET.SubElement(euler_angles, "field2", name="Angle_Y").text = str(separator.join(map(str,extract(1, connector_record.list_of_angles))))
                    ET.SubElement(euler_angles, "field3", name="Angle_Z").text = str(separator.join(map(str,extract(2, connector_record.list_of_angles))))

    tree = ET.ElementTree(root)

    tree.write("./output/RESULTS.xml", pretty_print=True)
