# Parsing PDB File

## Author:
Fadlilah Nur Hasanah (fadlilahhsnh@gmail.com)

## Introduction

The writing of this program is aimed to learn the pdb file format and how to extract information from it. 
This program gives summary of the amino acid compositions. By parsing the pdb file, one can obtain the residue positions
to do a lot of further analysis such as modelling. This program only uses those positions to calculate the most 
distance residues.
Additionally, one can see the summary of amino acid composition in the output. Hydrophobicity and amino acid charge 
information could help determine the structure, the function of the protein, and predict the protein's interaction with another proteins.

Commentary about pdb file format:
- PDB file can't be parsed using space/tab separation. This may work on some small protein, but not generally. Therefore line coordinate-based separation.

## Usage:

parsing_pdb.py (inputPDBfile)

1fcn.pdb included in this project is one pdb file example. One can use any other PDB Files
downloaded from pdb databank. (https://www.rcsb.org/)

## Output:
The output will be printed out on the terminal.
The output shows:
 - Amino acid composition of the given protein in percentage of each amino acid type
 - Hyprophobicity composition of the given protein
 - Atomic composition (count of atom N, C(alpha), O, and S)
 - Charged residue composition
 - Heteroatoms composition
 - Most distance residues