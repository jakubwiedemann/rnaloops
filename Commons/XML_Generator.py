#!/usr/bin/python
from lxml import etree as ET
from Commons.Utilities import extract, extract_XML
import glob


def xml_generate(list_of_records):

    separator = ', '
    root = ET.Element("File")
    for record in list_of_records:
        pdb_structure = ET.SubElement(root, "PDB_Structure")
        ET.SubElement(pdb_structure, "field1", name="PDB_ID").text = record.name

        for junction_record in record.list_of_junctions:
            junction = ET.SubElement(pdb_structure, "Junction")
            ET.SubElement(junction, "field1", name="Type_of_junction").text = junction_record.type
            ET.SubElement(junction, "field2", name="Junction_PDB_File").text = str(separator.join(map(str,junction_record.name_of_file)))
            #ET.SubElement(junction, "field3", name="Stem lengths").text = str(separator.join(map(str,junction_record.segment_length)))
            ET.SubElement(junction, "field3", name="DB_Notation").text = str(junction_record.sequence).rstrip()
            ET.SubElement(junction, "field4", name="Sequence").text = str(junction_record.list_of_segment_db).rstrip()

            stems = ET.SubElement(junction, "Stems")
            for stem_record in junction_record.list_of_stems:
                stem = ET.SubElement(stems, "Stem")
                ET.SubElement(stem, "field1", name="Duplex sequence").text = str(separator.join(map(str,stem_record.sequence))).rstrip()
                ET.SubElement(stem, "field2", name="Connectors").text = str(separator.join(map(str,sorted(stem_record.list_of_connectors))))
                ET.SubElement(stem, "field3", name="Duplex lengths").text = str(stem_record.segment_length)

            connectors = ET.SubElement(junction, "Connectors")
            for connector_records in junction_record.list_of_connectors:
                for connector_record in connector_records:
                    connector = ET.SubElement(connectors, "Connector")
                    ET.SubElement(connector, "field1", name="ID").text = str(connector_record.connector_id)
                    ET.SubElement(connector, "field2", name="Range").text = str(separator.join(map(str,connector_record.list_of_segments_ranges)))
                    ET.SubElement(connector, "field3", name="Length").text = str(connector_record.lengths_of_segments)
                    ET.SubElement(connector, "field4", name="Sequence").text = str(connector_record.list_of_segment_seq[1:-1]).rstrip()
                    #ET.SubElement(connector, "field5", name="DB_Notation").text = str(connector_record.list_of_segment_db)
                    ET.SubElement(connector, "field5", name="Planar_angle").text = str(separator.join(map(str,connector_record.planar_angle[0])))
                    #ET.SubElement(connector, "field5", name="Range_id").text = str(separator.join(map(str,connector_record.list_of_segments_ranges_id)))

                    euler_angles = ET.SubElement(connector, "Euler_Angles")
                    ET.SubElement(euler_angles, "field1", name="Angle_X").text = str(separator.join(map(str,extract_XML(0, connector_record.list_of_angles))))
                    ET.SubElement(euler_angles, "field2", name="Angle_Y").text = str(separator.join(map(str,extract_XML(1, connector_record.list_of_angles))))
                    ET.SubElement(euler_angles, "field3", name="Angle_Z").text = str(separator.join(map(str,extract_XML(2, connector_record.list_of_angles))))

    tree = ET.ElementTree(root)

    tree.write("./output/RESULTS.xml", pretty_print=True)

def generate_secondary_structure_summary(secondary_structure):
    separator = ', '
    root = ET.Element("File")
    for record in secondary_structure:
        pdb_structure = ET.SubElement(root, "PDB_Structure")
        ET.SubElement(pdb_structure, "field1", name="PDB_ID").text = record.name

        for junction_record in record.list_of_junctions:
            junction = ET.SubElement(pdb_structure, "Junction")
            ET.SubElement(junction, "field1", name="Type_of_junction").text = junction_record.type

            connectors = ET.SubElement(junction, "Connectors")
            for connector_records in junction_record.list_of_connectors:
                for connector_record in connector_records:
                    connector = ET.SubElement(connectors, "Connector")
                    ET.SubElement(connector, "field1", name="ID").text = str(connector_record.connector_id)
                    ET.SubElement(connector, "field2", name="Range").text = str(separator.join(map(str,connector_record.list_of_segments_ranges)))
                    ET.SubElement(connector, "field3", name="Length").text = str(connector_record.lengths_of_segments)
                    ET.SubElement(connector, "field4", name="Sequence").text = str(connector_record.list_of_segment_seq)
                    ET.SubElement(connector, "field5", name="DB_Notation").text = str(connector_record.list_of_segment_db)

def xml_generate_single_rec(record):

    separator = ', '
    root = ET.Element("File")

    pdb_structure = ET.SubElement(root, "PDB_Structure")
    ET.SubElement(pdb_structure, "field1", name="PDB_ID").text = record.name

    for junction_record in record.list_of_junctions:
        junction = ET.SubElement(pdb_structure, "Junction")
        ET.SubElement(junction, "field1", name="Type_of_junction").text = junction_record.type
        ET.SubElement(junction, "field2", name="Junction_PDB_File").text = str(separator.join(map(str,junction_record.name_of_file)))
        #ET.SubElement(junction, "field3", name="Stem lengths").text = str(separator.join(map(str,junction_record.segment_length)))
        ET.SubElement(junction, "field3", name="DB_Notation").text = str(junction_record.sequence).rstrip()
        ET.SubElement(junction, "field4", name="Sequence").text = str(junction_record.list_of_segment_db).rstrip()

        stems = ET.SubElement(junction, "Stems")
        for stem_record in junction_record.list_of_stems:
            stem = ET.SubElement(stems, "Stem")
            ET.SubElement(stem, "field1", name="Duplex sequence").text = str(separator.join(map(str,stem_record.sequence))).rstrip()
            ET.SubElement(stem, "field2", name="Connectors").text = str(separator.join(map(str,sorted(stem_record.list_of_connectors))))
            ET.SubElement(stem, "field3", name="Duplex lengths").text = str(stem_record.segment_length)

        connectors = ET.SubElement(junction, "Connectors")
        for connector_records in junction_record.list_of_connectors:
            for connector_record in connector_records:
                connector = ET.SubElement(connectors, "Connector")
                ET.SubElement(connector, "field1", name="ID").text = str(connector_record.connector_id)
                ET.SubElement(connector, "field2", name="Range").text = str(separator.join(map(str,connector_record.list_of_segments_ranges)))
                ET.SubElement(connector, "field3", name="Length").text = str(connector_record.lengths_of_segments)
                ET.SubElement(connector, "field4", name="Sequence").text = str(connector_record.list_of_segment_seq[1:-1]).rstrip()
                #ET.SubElement(connector, "field5", name="DB_Notation").text = str(connector_record.list_of_segment_db)
                ET.SubElement(connector, "field5", name="Planar_angle").text = str(separator.join(map(str,connector_record.planar_angle[0])))
                #ET.SubElement(connector, "field5", name="Range_id").text = str(separator.join(map(str,connector_record.list_of_segments_ranges_id)))

                euler_angles = ET.SubElement(connector, "Euler_Angles")
                ET.SubElement(euler_angles, "field1", name="Angle_X").text = str(separator.join(map(str,extract_XML(0, connector_record.list_of_angles))))
                ET.SubElement(euler_angles, "field2", name="Angle_Y").text = str(separator.join(map(str,extract_XML(1, connector_record.list_of_angles))))
                ET.SubElement(euler_angles, "field3", name="Angle_Z").text = str(separator.join(map(str,extract_XML(2, connector_record.list_of_angles))))

    tree = ET.ElementTree(root)

    tree.write("./output/single_records/"+record.name +".xml", pretty_print=True)


def run():
    xml_files = glob.glob("./output/single_records/*.xml")
    xml_element_tree = None
    for xml_file in xml_files:
        data = ET.parse(xml_file).getroot()
        # print ET.tostring(data)
        for result in data.iter('File'):
            if xml_element_tree is None:
                xml_element_tree = data
                insertion_point = xml_element_tree.findall("./File")[0]
            else:
                insertion_point.extend(result)
    if xml_element_tree is not None:

        tree = ET.ElementTree(xml_element_tree)

        tree.write("./output/merged.xml", pretty_print=True)

def newRunRun():

    xml_files = glob.glob("./output/single_records/*.xml")
    f = open("merged1.xml", "w")
    f.write("<File>\n")
    for xmlFile in xml_files:
        f1 = open(xmlFile,"r")
        lines = f1.readlines()
        lines = lines[1:-1]
        f.writelines(lines)
    f.write("</File>")
    f.close()
