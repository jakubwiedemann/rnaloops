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


def find_junction(db_sequence, junction_start, common_stem_length_calc = True):
    stems_identified = 0
    list_of_pairs = []
    if db_sequence[junction_start] == "(":
        last = find_matching_parenthesis(db_sequence, junction_start)
        list_of_pairs.append([junction_start+1, last+1])
        junction_start += 1
        stems_identified += 1

        while junction_start < last:
            if db_sequence[junction_start] == "(":
                end = find_matching_parenthesis(db_sequence, junction_start)
                list_of_pairs.append([junction_start+1, end+1])

                junction_start = end + 1
                stems_identified += 1

            else:
                junction_start += 1
    else:
        junction_start += 1
    if stems_identified > 2:
        if common_stem_length_calc:
            return True, stems_identified, list_of_pairs, common_stem_length(db_sequence, list_of_pairs)
        else:
            return True, stems_identified, list_of_pairs
    else:
        return False, stems_identified, list_of_pairs


def common_stem_length(db_sequence, list_of_pairs):
    length_table = []
    pair_length_table = []
    for i,pair in enumerate(list_of_pairs):
        start = pair[0]-1
        open_count = 0
        end = pair[1]-1
        end_count = 0
        if i == 0:
            while db_sequence[start] == '(' and start >= 0:
                start -=1
                open_count += 1
        else:
            while db_sequence[start] == '(':
                start +=1
                open_count += 1
        length_table.append(open_count)
        if i == 0:
            while db_sequence[end] == ')' and end < len(db_sequence)-1:
                end +=1
                end_count += 1
        else:
            while db_sequence[end] == ')':
                end -=1
                end_count += 1
        length_table.append(end_count)
        pair_length_table.append([open_count, end_count])
    pair_common_length_table = [min(pair) for pair in pair_length_table]
    return pair_common_length_table

