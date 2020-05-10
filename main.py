#!/usr/bin/python
import config
import os
from Modes.MultipleCustom import multiple_file_custom_mode, multiple_files_from_pdbee


if __name__ == "__main__":

    # Setting directories
    if not os.path.exists('./output'):
        os.makedirs('output')
    if not os.path.exists('./output/structures'):
        os.makedirs('./output/structures')
    if not os.path.exists('./PDB_files'):
        os.makedirs('PDB_files')

    # Run main application
    if config.mode.upper() == 'MULTIPLE_CUSTOM':
        multiple_file_custom_mode(config.path_to_dotbracket_files)
    elif config.mode.upper() == 'MULTIPLE':
        multiple_files_from_pdbee(config.path_to_dotbracket_files)
    else:
        print('Select valid mode')
