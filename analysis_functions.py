#1. Imports

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

#2. Functions

    #2.1. Minimum heavy atom distance between two structural features: takes an MDA Universe and two residue specifiers as input. Outputs np.array([time, distance])

def calc_dist(u, feature_A_specifier, feature_B_specifier):
    """
    u (MDAnalysis.Universe object)
    feature_A_specifier (MDAnalysis selection string)
    feature_B_specifier (MDAnalysis selection string)
    """
    
    feature_A, feature_B = u.select_atoms(feature_A_specifier + " and not type H"), u.select_atoms(feature_B_specifier + " and not type H")   
    
    print("Calculating SB distance between " + str(feature_A[0].residue) + " and " + str(feature_B[0].residue) + ".")
    
    time = []
    ha_min_dist = []
    for ts in u.trajectory:
        ha_dist_ts = []
        for atom_A in feature_A:
            for atom_B in feature_B:
                ha_dist = np.linalg.norm(atom_A.position-atom_B.position)
                ha_dist_ts.append(ha_dist)
        ha_min_dist.append(min(ha_dist_ts))
        
        time.append(ts.time)
    
    time_np = np.transpose(np.array([time]))
    ha_min_dist_np = np.transpose(np.array([ha_min_dist]))
    
    fct_output = np.concatenate((time_np, ha_min_dist_np), axis=1)
    
    return fct_output

    #2.2. Calculates backbone H-bonds between two sets of residues.  Outputs np.array([time, number_of_h-bonds])

def calc_h_bonds(u, sel1str, sel2str,  dssp_cutoff=-0.5):
    """
    u (MDAnalysis.Universe object)
    sel1str (MDAnalysis selection string)
    sel2str (MDAnalysis selection string)
    dssp_cutoff (float, in kcal/mol)
    """
    
    sel1=u.select_atoms(sel1str)
    sel2=u.select_atoms(sel2str)
    h_b_trj=[]
    for ts in u.trajectory:
        n_hbonds=0
        for residue1 in sel1.residues:
            for residue2 in sel2.residues:
                try:
                    donor_atomgroup=residue1.atoms
                    acceptor_atomgroup=residue2.atoms
                    #Use DSSP definition of an H-bond (https://doi.org/10.1002/bip.360221211)
                    #Select donor and acceptor atoms
                    donor_N = donor_atomgroup.select_atoms("name N")
                    donor_H = donor_atomgroup.select_atoms("name HN H H01 H02 H03")
                    acceptor_O = acceptor_atomgroup.select_atoms("name O")
                    acceptor_C = acceptor_atomgroup.select_atoms("name C")
                    #Calculate N-C, N-O, H-O, and H-C vector
                    NC_vec = donor_N.positions[0]-acceptor_C.positions[0]
                    NO_vec = donor_N.positions[0]-acceptor_O.positions[0]
                    HC_vec = donor_H.positions[0]-acceptor_C.positions[0]
                    HO_vec = donor_H.positions[0]-acceptor_O.positions[0]
                    #Compute DSSP energy
                    dssp_E = 0.42*0.20*332.0*(1/np.linalg.norm(NO_vec)+1/np.linalg.norm(HC_vec)-1/np.linalg.norm(HO_vec)-1/np.linalg.norm(NC_vec)) #kcal/mol
                    #Check criterium
                    if dssp_E <= dssp_cutoff:
                        n_hbonds+=1
                    else:
                        pass
                except: #A common exception would be a proline in the list of potential donor residues, which does not have an backbone amid H
                    pass
                try:
                    donor_atomgroup=residue2.atoms
                    acceptor_atomgroup=residue1.atoms
                    #Use DSSP definition of an H-bond (https://doi.org/10.1002/bip.360221211)
                    #Select donor and acceptor atoms
                    donor_N = donor_atomgroup.select_atoms("name N")
                    donor_H = donor_atomgroup.select_atoms("name HN H H01 H02 H03")
                    acceptor_O = acceptor_atomgroup.select_atoms("name O")
                    acceptor_C = acceptor_atomgroup.select_atoms("name C")
                    #Calculate N-C, N-O, H-O, and H-C vector
                    NC_vec = donor_N.positions[0]-acceptor_C.positions[0]
                    NO_vec = donor_N.positions[0]-acceptor_O.positions[0]
                    HC_vec = donor_H.positions[0]-acceptor_C.positions[0]
                    HO_vec = donor_H.positions[0]-acceptor_O.positions[0]
                    #Compute DSSP energy
                    dssp_E = 0.42*0.20*332.0*(1/np.linalg.norm(NO_vec)+1/np.linalg.norm(HC_vec)-1/np.linalg.norm(HO_vec)-1/np.linalg.norm(NC_vec)) #kcal/mol
                    #Check criterium
                    if dssp_E <= dssp_cutoff:
                        n_hbonds+=1
                    else:
                        pass
                except:
                    pass
        h_b_trj.append([ts.time, n_hbonds])
    return h_b_trj