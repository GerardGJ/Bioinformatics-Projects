#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os
import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset

NACCESS_BIN = '/home/gerard/Desktop/Biophys/Exercise/BioPhysics/soft/NACCESS/naccess'


parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.dirname(__file__) + '/data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('pdbqt_file',help='Input PDBQT', type=open)

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing PDB', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Getting Charges and Atom type from PDBQT
print('Parsing PDBQT', args.pdbqt_file.name)
params=[{}]

for line in args.pdbqt_file:
    line = line.rstrip()
    params.append({'charge': line[69:76], 'type': line[77:].replace(' ','')})

total_charge = 0.
for at in st.get_atoms():
    at.xtra['atom_type'] = params[at.serial_number]['type']
    at.xtra['charge'] = float(params[at.serial_number]['charge'])
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    total_charge += at.xtra['charge']
print('Total Charge: {:8.2f}'.format(total_charge))

# Calculating surfaces
# Srf goes to .xtra['EXP_NACCESS'] field
srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BIN)


# Simple Test for atom fields
print(vars(st[0]['A'][2]['O']))
print(st[0]['A'][2]['O'].xtra['atom_type'])
vdw_p = st[0]['A'][2]['O'].xtra['vdw'].eps
print(vdw_p)
#print(atom.get_name(),atom.xtra['charge'])


#Exercise 3 for just one Residue:
electrostatic_acumulation = [0,0,0]
total_vdw = 0
for atom1 in st[0]['A'][2].get_atoms(): #We get one atom
    for atom2 in st[0]['A'][2].get_atoms(): #We get the second atom
        r = atom2-atom1
        if r != 0: #Here me make sure we dont compute the interactions between the same atom
            er_list = [1,80,(86.9525/(1-7.7839*np.exp(-0.3153*r)))-8.5525] #Calculate the Aproximation
            e_electrostatic_list = []
            for er in er_list: #Calculate the electrostatic interaction for each pair of atoms and for each constant
                elec = 332.16*((atom1.xtra['charge']*atom2.xtra['charge'])/(er*r))
                e_electrostatic_list.append(elec)
            #print('Electrostati interactions:',atom1.get_name(),atom2.get_name())
            #print(e_electrostatic_list)
            for x in range(3):
                electrostatic_acumulation[x] += e_electrostatic_list[x]
            #print('Vdw:')
            eij = (atom1.xtra['vdw'].eps * atom2.xtra['vdw'].eps)**0.5
            sij = (atom1.xtra['vdw'].sig * atom2.xtra['vdw'].sig)**0.5
            evdw = 4 * eij * (((sij/r)**12)-((sij/r)**6)) #Here we calculate the Vdw
            total_vdw += evdw
            #print(evdw)
            #print('')
print('')
print('Electrostatic interaction:',electrostatic_acumulation)
print('vdw:',evdw)


#3 Electrostatic b/w atoms:
electrostatic_acumulation = [0,0,0]
total_vdw = 0
explored_res = []
for res1 in st.get_residues(): #We get one residue
    res_electro_acu = [0,0,0]
    res_vdw = 0
    for res2 in st.get_residues(): #We get the other residue
        if res1 != res2 and res2 not in explored_res: #As we only calculate the energies for atoms not in the same residue we make sure they are not the same
            for atom1 in res1.get_atoms():#We get an atom from the first res
                for atom2 in res2.get_atoms():#We get an atom from the second res
                    r = atom2-atom1
                    er_list = [1,80,(86.9525/(1-7.7839*np.exp(-0.3153*r)))-8.5525] #Calculate the Mehler-Solmajer 
                    e_electrostatic_list = []
                    for er in er_list:
                        elec = 332.16*((atom1.xtra['charge']*atom2.xtra['charge'])/(er*r)) #Calculate the energy for each of the constants
                        e_electrostatic_list.append(elec)
                    for x in range(3):
                        electrostatic_acumulation[x] += e_electrostatic_list[x] #add them
                        res_electro_acu[x] += e_electrostatic_list[x]
                    eij = np.sqrt(atom1.xtra['vdw'].eps * atom2.xtra['vdw'].eps) #get epsilon ij
                    sij = np.sqrt(atom1.xtra['vdw'].sig * atom2.xtra['vdw'].sig)#get sigma ij
                    evdw = 4 * eij * (((sij/r)**12)-((sij/r)**6)) #Calculate the vdw
                    total_vdw += evdw #add them to the total of vdw energies
                    res_vdw += evdw
    print(res1.get_id())
    print(res_electro_acu)
    print(res_vdw)
    explored_res.append(res1)
print('')
print('')
print(electrostatic_acumulation)
print(total_vdw)


#Exercise 4 for all res:
total_solvation = 0
for res in st.get_residues():
    for atom in res.get_atoms():
        if atom.element != 'H':
            sig = atom.xtra['vdw'].sig
            asa = atom.xtra['EXP_NACCESS']
            total_solvation += sig*float(asa)
            #print(atom.xtra['EXP_NACCESS'], atom.element)
print('Solvation_energie:',total_solvation)
