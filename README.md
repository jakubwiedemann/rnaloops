# nWayJunction_release

Usage example

Before running, user have to set up config file. As a default we have 2 params: 
*mode (set to multiple),
*path_to_dotbracket_files (set to ./dotbracket_files/test_multiple/). 
In ./dotbracket_files/test_multiple/ folder is placed a dot-bracket structure generated for 1E8O by RNApdbee (http://rnapdbee.cs.put.poznan.pl/) that can be use as an example. 1E80 contains a 3-way junction that will be found by the nWayJunction tool.

After running, in ./output folder program generates RESULT.xml file and folder /structures with 3D structures of the found n-way junctions.

Workflow description

As an input tool uses a file with the dot-bracket representation of the structure (to get this file from the 3D structure you can use RNApdbee).

Based on the dot-bracket notation algorithm finds an n-way junction. In the next step, the algorithm looks for PDB file in folder, if the file does not exist the structure is retrieved automatically and saved in the mentioned folder. For identified junctions based on 3D structure algorithm represents outgoing stems as a vector in 3D space and for pairs of stems calculates the value of planar angle and Euler angles. The final step is the generation of RESULT.xml file that contains all results.
