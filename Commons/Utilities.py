#!/usr/bin/python
import os


def pairs_generator(lst):
    n = len(lst)
    for i in range(n):
        yield lst[i],lst[(i+1)%n]


def is_non_zero_file(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) > 0
