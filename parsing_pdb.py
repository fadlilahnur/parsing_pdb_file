import sys
from math import sqrt


# class to save information of each amino acid from pdb file
class AminoAcid:
    def __init__(self, atom, amino_acid_type, chain_alphabet, residue_number, location):
        self.atom = atom
        self.type = amino_acid_type
        self.chain_alphabet = chain_alphabet
        self.residue_number = residue_number
        self.location = location

    def get_atom(self):
        return self.atom

    # saving location
    def get_x(self):
        return self.location.get('x')

    def get_y(self):
        return self.location.get('y')

    def get_z(self):
        return self.location.get('z')

    # amino acid type
    def get_type(self):
        return self.type

    # chain information
    def get_chain(self):
        return self.chain_alphabet

    def get_residue_number(self):
        return self.residue_number


# function to calculate distance of two atoms
# param atom1 and atom2 are two atoms, their distance are to be computed
# return the distance of two given atoms
def atom_distance(atom1, atom2):
    # calculate 3D distance with their coordinates
    return sqrt(pow(atom1.get_x() - atom2.get_x(), 2)
                + pow(atom1.get_y() - atom2.get_y(), 2)
                + pow(atom1.get_y() - atom2.get_z(), 2))


# calculate the amino acids compositions
# param pdb_object: opened pdb file
# returns: amino_acid_quantities, hydrophobic_quantity, hydrophilic_quantity, atomic_quantities, charged_quantities,
#         het_atom_composition, chain_map_amino_acids (map chain alphabet and the list of amino acids in that chain)
def amino_acids_composition(pdb_object):
    # initial dictionary for 20 Proteins and their corresponding quantity
    amino_acid_quantities = dict()
    # hydrophobicity need to be cross-checked
    # list of hydrophobic and hydrophilic residues as data base
    hydrophobic = ({'ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'VAL', 'MET'})
    hydrophilic = ({'ASN', 'GLN', 'GLY', 'LYS', 'HIS', 'ASP', 'PRO', 'GLU', 'ARG', 'SER', 'THR', 'TRP', 'TYR'})
    # counter for the hydrophobicity categories
    hydrophilic_quantity = 0
    hydrophobic_quantity = 0
    # container for the atomic quantities, maps the atom name and its corresponding quantity
    atomic_quantities = dict()
    # charged amino acids data base
    positive_charged = ({'LYS', 'HIS', 'ARG'})
    negative_charged = ({'ASP', 'GLU'})
    # container to save the charges quantity
    charged_quantities = dict({'positive': 0,
                               'negative': 0})
    # container that maps the hetero atom with its corresponding quantity
    het_atom_composition = dict()
    # container that store each amino acids of all chains (key : chain alphabet, value : list of amino acids)
    chain_map_amino_acids = dict({'A': []})
    chains = ['A']

    last_residue = 0
    last_chain = "A"
    location_getter = dict({'x': 6, 'y': 7, 'z': 8})
    # open file and read it per line
    pdb_content = pdb_object.readlines()
    for line in pdb_content:
        # find if the line contains ATOM
        current_line = line.split()
        if current_line[0] == 'ATOM':
            atom_line = current_line
            # differentiate if there is a 4 digits amino acid
            if len(atom_line[2]) == 7:
                current_residue = atom_line[2][4] + atom_line[2][5] + atom_line[2][6]
                atom_type = atom_line[2][0] + atom_line[2][1] + atom_line[2][2]
                location_getter = dict({'x': 5, 'y': 6, 'z': 7})
                atom = str(atom_line[10])
                # differentiate the case, where the residue number is above 999
                if len(atom_line[3]) == 1:
                    residue_number = atom_line[4]
                    current_chain = atom_line[3]
                else:
                    residue_number = atom_line[3].removeprefix()
                    current_chain = atom_line[3][0]
                    location_getter = dict({'x': 4, 'y': 5, 'z': 6})
            else:
                if len(atom_line[3]) == 4:
                    current_residue = atom_line[3][1] + atom_line[3][2] + atom_line[3][3]
                    atom_type = str(atom_line[2])
                    atom = str(atom_line[11])
                else:
                    current_residue = str(atom_line[3])
                    atom_type = str(atom_line[2])
                    atom = str(atom_line[11])
                    # differentiate the case, where the residue number is above 999
                if len(atom_line[4]) == 1:
                    residue_number = atom_line[5]
                    current_chain = atom_line[4]
                else:
                    residue_number = atom_line[4].removeprefix()
                    current_chain = atom_line[4][0]
                    location_getter = dict({'x': 5, 'y': 6, 'z': 7})

            # counting atoms
            if atom is not None:
                if atom not in atomic_quantities:
                    atomic_quantities[atom] = 0
                atomic_quantities[atom] += 1

            # update value for new chain
            if last_chain != current_chain:
                last_residue = 0
                last_chain = current_chain
                chains.append(current_chain)
                chain_map_amino_acids[current_chain] = []

            # storing CA locations

            if atom_type == 'CA':
                chain_map_amino_acids.get(current_chain).append(
                    AminoAcid(atom_type, current_residue, current_chain,
                              residue_number,
                              {'x': float(current_line[location_getter.get('x')]),
                               'y': float(current_line[location_getter.get('y')]),
                               'z': float(current_line[location_getter.get('z')])}))
                # output_object.writelines("current res" + current_residue + "\n"
                #                          + residue_number + "\n"
                #                          + current_line[location_getter.get('x')] + "\n")

            # if the current chain is different from the previous chain it should update the residue number
            if residue_number != last_residue:
                last_residue = residue_number
                # counting amino acids
                if current_residue not in amino_acid_quantities:
                    amino_acid_quantities[current_residue] = 0
                amino_acid_quantities[current_residue] += 1
                # counting hydrophobicity-based categories
                if current_residue in hydrophobic:
                    hydrophobic_quantity += 1
                elif current_residue in hydrophilic:
                    hydrophilic_quantity += 1
                # counting charges
                if current_residue in positive_charged:
                    charged_quantities['positive'] += 1
                elif current_residue in negative_charged:
                    charged_quantities['negative'] += 1
        # counting hetero atoms
        if current_line[0] == 'HETATM':
            # check if the chain alphabet and residue number became one (in case of residue number >999)
            if len(current_line[4]) == 1:
                residue_number = current_line[5]
            else:
                prefix = current_line[4][0]
                residue_number = current_line[4].removeprefix(prefix)
            if residue_number != last_residue:
                last_residue = residue_number
                current_residue = str(current_line[3])
                # exclude water (HOH)
                if current_residue != 'HOH':
                    if current_residue not in het_atom_composition:
                        het_atom_composition[current_residue] = 0
                    het_atom_composition[current_residue] += 1

    return amino_acid_quantities, hydrophobic_quantity, hydrophilic_quantity, atomic_quantities, charged_quantities, \
        het_atom_composition, chain_map_amino_acids


# function to calculate most distance atom in an atom list
# param atom_list is the list where the atoms are stored
# returns the pair of most distance atom and the distance
def most_distance_atom(atom_list):
    distance = 0
    pair = ()
    for atom1 in atom_list:
        for atom2 in atom_list:
            tmp = atom_distance(atom1, atom2)
            if tmp > distance:
                distance = tmp
                pair = (atom1, atom2)
    return pair, round(distance, 2)


# main
if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: parsing_pdb.py <inputPDBfile>")
    # allowing only pdb file as an input
    if not sys.argv[1].endswith('.pdb'):
        sys.exit("usage: parsing_pdb.py <inputPDBfile>")

    # get the file from the cmd line argument
    pdb_file = sys.argv[1]

    # printing the result to output file
    with open(pdb_file, "r") as pdb_file_object:
        amino_acids_composition = amino_acids_composition(pdb_file_object)
        print("Amino Acid Composition: \n")
        total_amino_acids = 0
        aa_composition = amino_acids_composition[0]
        # summing all amino acids quantities
        for key in aa_composition:
            total_amino_acids += aa_composition.get(key)
        # printing out the amino acid composition quantities and percentages
        for key in aa_composition:
            quantity = aa_composition.get(key)
            percentage = round(quantity / total_amino_acids * 100, 1)
            print(key + " " + str(quantity) + " " + str(percentage) + "% \n")

        # printing out hydrophobicity based composition
        hydrophobic_percentage = round(amino_acids_composition[1] / total_amino_acids * 100, 1)
        hydrophilic_percentage = round(amino_acids_composition[2] / total_amino_acids * 100, 1)
        print("\nAmino acid composition categorized on the basis of hydrophobicity:\n")
        print("Hydrophobic " + str(hydrophobic_percentage) + "%\n")
        print("Hydrophilic " + str(hydrophilic_percentage) + "%\n \n")

        # printing out atomic composition
        atomic_dictionary = amino_acids_composition[3]
        print("Atomic composition:\n")
        total_atom = 0
        # for key in atomic_dictionary:
        #     total_atom += atomic_dictionary.get(key)
        for key in atomic_dictionary:
            quantity = atomic_dictionary.get(key)
            # percentage = round(quantity / total_atom * 100)
            print(str(key) + " " + str(quantity) + "\n")

        # printing out charged amino acids quantities
        charged_amino_acid_quantities = amino_acids_composition[4]

        print("\nTotal number of positively charged residues: "
              + str(charged_amino_acid_quantities['positive']) + "\n")
        print("Total number of negatively charged residues: "
              + str(charged_amino_acid_quantities['negative']) + "\n \n")
        # printing out hetero atom composition
        print("Number of heteroatoms: "
              + str(len(amino_acids_composition[5])) + "\n")
        for key in amino_acids_composition[5]:
            print(key + "\n")

        # looking for the most distance pair
        amino_acid_lists = amino_acids_composition[6]
        most_distance = 0
        most_distance_pair = ()
        for aa_list in amino_acid_lists:
            current_pair = most_distance_atom(amino_acid_lists.get(aa_list))
            if current_pair[1] > most_distance:
                most_distance = current_pair[1]
                most_distance_pair = current_pair[0]
        # print out the most distance pair and its distance in Angstrom
        print("Distance between most distant residues " +
              str(most_distance_pair[0].get_residue_number()) +
              " " + most_distance_pair[0].get_chain() + " " +
              "and " +
              str(most_distance_pair[1].get_residue_number()) + " " +
              most_distance_pair[1].get_chain() + " :" + str(most_distance)
              + " Angstrom\n")

    pdb_file_object.close()
