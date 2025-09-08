# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 19:08:02 2020

@author: kevin

Copied from CPyT for the course CIEM 2301: Offshore Geotechnics
github.com/triaduct/cpyt
github.com/triaduct/offshoregeotechnics
"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt 

"""
NOTE: 
    
    -ve implies downwards direction
"""
#%%
def koppejan(cpt, pile_dia, target_depth):
    """
    4D/8D method according to Koppejan
    """
    ## Check if the required CPT data exists
    if len(cpt[cpt.z < target_depth - 4*pile_dia]) == 0:
        raise ValueError("ERROR: CPT doesn't cover depth required for averaging technique. Choose a more suitable CPT in the .txt file\n")
    
    # Find the minimum qcavg (dependent on endpoint of qc1)
    qc_avg = 9999            # Initialise qcavg
    
    for endpoint in cpt[(cpt.z <= target_depth-0.7*pile_dia) & (cpt.z >= target_depth - 4*pile_dia)].z:    # Identify endpoint where qcavg is minimal
        ## qc1
        """Mean of values until endpoint"""
        cpt1 = cpt[(cpt.z <= target_depth) & (cpt.z >= endpoint)]
        qc1 = cpt1.qc.mean()
    
        ## qc2
        """Min. path rule from 4D to 0.7D"""
        cpt2 = cpt[(cpt.z <= target_depth) & (cpt.z >= endpoint)]
        cpt2 = cpt2.iloc[::-1]                                        # Flip cpt for use in expanding minimum
        cpt2["min_path"] = cpt2.qc.expanding().min()                  # Expanding minimum
        qc2 = cpt2.min_path.mean()
        
        ## qc3 
        """Min. path rule from pile tip to 8D"""
        start = cpt2.min_path.iloc[-1]                               # Start min path from the end of min path for qc2
        cpt3 = cpt[(cpt.z >= target_depth) & (cpt.z <= target_depth + 8*pile_dia)]   # Search for min within this cpt
        cpt3 = cpt3.iloc[::-1]                                        # Flip cpt for use in expanding minimum
        cpt3["min_path"] = np.nan
        cpt3.qc.iloc[0] = start                                      # Continue on from path from qc2
        cpt3.min_path = cpt3.qc.expanding().min()                    # Expanding minimum
        qc3 = cpt3.min_path.mean()
        
        ## qc_avg
        qc_avg_temp = 0.5*(0.5*(qc1+qc2) + qc3)
        if qc_avg_temp < qc_avg:
            qc_avg = qc_avg_temp
            qc1_min = qc1
            qc2_min = qc2
            qc3_min = qc3
            qc2_min_path = np.array([cpt2.z,cpt2.min_path])
            qc3_min_path = np.array([cpt3.z,cpt3.min_path])
            endpoint_final = endpoint
            
    return qc_avg


def lcpc(cpt, pile_dia, target_depth):

    D15 = cpt[(cpt.z >= target_depth - 1.5*pile_dia) & (cpt.z <= target_depth + 1.5*pile_dia)]   # Region across which the averaging will be applied
    qc_avg = D15.qc.mean()
    D15.loc[D15.qc >= 1.3*qc_avg,"qc"] = 1.3*qc_avg        # Truncate qc values
    D15.loc[D15.qc <= 0.7*qc_avg,"qc"] = 0.7*qc_avg        # Truncate qc values
    qc_avg = D15.qc.mean()
    
    return qc_avg

def de_boorder(cpt: pd.DataFrame, pile_dia: float, target_depth: float) -> float:
    HD_a = 6.5    # Distance over which the cosine function is applied above the pile
    HD_b = 10.5   # Distance over which the cosine function is applied below the pile
    f = 13.5      # Damping factor
    s_a = 0.56    # Reshapes the weight related to the stiffness ratio
    s_b = 0.79    # Reshapes the weight related to the stiffness ratio

    cpt = cpt.loc[
        (cpt["z"] <= target_depth + HD_a * pile_dia) &
        (cpt["z"] >= target_depth - HD_b * pile_dia)
    ].copy()
  
    cpt.loc[:, "HD"] = np.where(cpt["z"] >= target_depth, HD_a, HD_b)       # HD: above vs below pile tip
    cpt.loc[:, "x"] = np.abs((target_depth - cpt["z"]) / (pile_dia * cpt["HD"]))    # Distance factor
    cpt.loc[:, "w1"] = np.exp(-f * cpt["x"]) * np.cos(0.5 * np.pi * cpt["x"])   # cosine damping
    near_z_ix = cpt["z"].sub(target_depth).abs().idxmin()
    qc_tip = cpt.loc[near_z_ix, "qc"]                                           # qc at depth closest to pile tip
    cpt.loc[:, "s"] = np.where(cpt["z"] >= target_depth, s_a, s_b)
    cpt.loc[:, "w2"] = (qc_tip / cpt["qc"]) ** cpt["s"]         # Stiffness ratio weighting
    cpt.loc[:, "w3"] = cpt["w1"] * cpt["w2"]                    # Total weightingghting

    qc_w = cpt["qc"] * cpt["w3"] / cpt["w3"].sum()
    qc_avg = qc_w.sum()

    return qc_avg
