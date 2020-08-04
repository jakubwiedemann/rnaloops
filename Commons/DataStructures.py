#!/usr/bin/python


class Record:
    def __init__(self, name, method='', chains=[], segment_length=[]):
        self.name = name
        self.method = method
        self.chains = chains
        self.segment_length = segment_length


class Stem:
    def __init__(self, segment_length=1, list_of_junctions=[]):
        self.segment_length = 1
        self.list_of_junctions = list_of_junctions


class Junction:
    def __init__(self):
        self.type = ''
        self.list_of_segments_ranges = None
        self.lengths_of_segments = None
        self.list_of_segment_seq = None
        self.list_of_segment_db = None
        self.list_of_angles = None
        self.planar_angle = None
        self.name_of_file = None
