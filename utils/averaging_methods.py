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


def de_boorder(cpt, pile_dia, target_depth):
    HD_a = 8.3              # Distance over which the cosine function is applied above the pile
    HD_b = 15.5             # "..." below the pile
    f = 13.5                # Damping factor
    s = 0.9                 # Reshapes the weight related to the stiffness ratio
    
    cpt = cpt.loc[(cpt.z <= target_depth + HD_a*pile_dia) & (cpt.z >= target_depth - HD_b*pile_dia)]
    cpt["HD"] = np.nan
    cpt.loc[cpt.z >= target_depth,"HD"] = HD_a
    cpt.loc[cpt.z <= target_depth, "HD"] = HD_b
    cpt["x"] = abs((target_depth-cpt.z)/(pile_dia*cpt.HD))
    cpt["w1"] = np.exp(-f*cpt.x)*np.cos(0.5*np.pi*cpt.x)   # First weight relating to the cosine dampening function and distance to pile tip
    near_z_ix = cpt.z.sub(target_depth).abs().idxmin()       # Index of row with depth closest to pile depth
    qc_tip = cpt.loc[near_z_ix, "qc"]
    cpt["w2"] = (qc_tip/cpt.qc)**s                        # Wegiht of one point related to stiffness ratio
    cpt["w3"] = cpt.w1*cpt.w2                              # Total weight of qc at one point
    
    qc_w = cpt.qc*cpt.w3/(cpt.w3.sum())
    qc_avg = qc_w.sum()
    
    return qc_avg


def boulanger_dejong(cpt, pile_dia, target_depth):
    """
    Method prescribed in Boulanger & DeJong (2018). Recommend averaging method
    in the Unified pile design method (Lehane et al., 2022)
    
    :dia:           Diameter of penetrometer.
                    NOTE: The Boulanger & DeJong method was not explicitly formulated for 
                    piles (namely scale effects). This was the further research of de 
                    Lange (2017), de Boorder (2021) and so on.
    :target_depth:  depth for calcualtion, with respect to cpt.z (-ve => downwards)       

    """
    z50_ref=4.0
    mz=3.0
    m50=0.5
    mq=2
            
    qt_tip = cpt.qc.iloc[(cpt.z-target_depth).abs().argsort()[:2]]      # Cone resistance at :target_depth:
    cpt[["C1","C2","z50","w1","w2","w1w2","wc"]] = np.nan
    cpt["z_norm"] = -1*(cpt.z-target_depth)/pile_dia    # Normalised depth. -ve values are above the pile tip, +ve are below
    
    # Equation 6
    cpt.C1.loc[cpt.z_norm >= 0] = 1
    cpt.C1.loc[(cpt.z_norm >= -4) & (cpt.z_norm < 0)] = 1 + cpt.z_norm/8
    cpt.C1.loc[cpt.z_norm < -4] = 0.5
    
    cpt.C2.loc[cpt.z_norm >= 0] = 1         # Equal to unity below the pile tip (+ve numbers imply downwards)
    cpt.C2.loc[cpt.z_norm <= 0] = 0.8       # Equal to 0.8 above the pile tip

    cpt.z50 = 1+2*(cpt.C2*z50_ref - 1)*(1-(1/(1+(qt_tip/cpt.qc)**m50)))     # Equation 7
    
    cpt.w1 = cpt.C1/(1+(cpt.z_norm/z50_ref)**mz)        # Equation 5
    cpt.w2 = np.sqrt(2/(1+(cpt.qc/qt_tip)**mq))         # Equation 8
    cpt.w1w2 = cpt.w1 * cpt.w2                          
    cpt.wc = (cpt.w1*cpt.w2)/cpt.w1w2.sum()             # Equation 3
    qc_w = (cpt.qc*cpt.wc)/(cpt.wc.sum()) 
    qc_avg = qc_w.sum()
    
    return qc_avg
    
