#!/usr/bin/python
import os


def find_matching_parenthesis(text, opening_position):
    closing_position = opening_position
    counter = 1
    while counter > 0:
        closing_position += 1
        char = text[closing_position]
        if char == "(":
            counter += 1
        if char == ")":
            counter -= 1
    return closing_position


def find_junction(text, junction_start):
    stems_identified = 0
    list_of_pairs = []
    if text[junction_start] == "(":
        last = find_matching_parenthesis(text, junction_start)
        list_of_pairs.append([junction_start+1, last+1])
        junction_start += 1
        stems_identified += 1

        while junction_start < last:
            if text[junction_start] == "(":
                end = find_matching_parenthesis(text, junction_start)
                list_of_pairs.append([junction_start+1, end+1])

                junction_start = end + 1
                stems_identified += 1

            else:
                junction_start += 1
    else:
        junction_start += 1
    if stems_identified > 2:
        return True, stems_identified, list_of_pairs
    else:
        return False, stems_identified, list_of_pairs
