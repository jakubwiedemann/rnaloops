#!/usr/bin/python
import config
import os
from Modes.MultipleCustom import multiple_file_custom_mode


if __name__ == "__main__":

    # Setting directories
    if not os.path.exists('./output'):
        os.makedirs('output')
    if not os.path.exists('./output/structures'):
        os.makedirs('./output/structures')
    if not os.path.exists('./PDB_files'):
        os.makedirs('PDB_files')

    # Run main application
    multiple_file_custom_mode(config.path_to_dotbracket_files)
