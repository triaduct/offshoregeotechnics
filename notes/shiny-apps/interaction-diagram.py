# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 13:59:07 2025

@author: kduffy
"""


from shiny import App, render, ui
import numpy as np
import matplotlib.pyplot as plt

def vulpe_2014_hm_plane(v, D, d, su_mudline, k):
    """
    Returns the failure surface in terms of h/h_star and m/m_star
    """
    # === Derived parameters ===
    su0 = su_mudline + k*d
    embedment_ratio = d / D

    # === Normalised load space ===
    n = 1000  # higher resolution for better contours
    h = np.linspace(-2, 2, n)
    m = np.linspace(-2, 2, n)
    h_grid, m_grid = np.meshgrid(h, m, indexing='ij')

    # === Fitting functions ===
    q = 4.69
    p = 2.12
    h_star = 1 - v**q
    m_star = 1 - v**p

    alpha = 2.28 - 1.03 * embedment_ratio if v > 0.5 else 2.55 - 1.43 * embedment_ratio
    beta  = 0.05 - 1.15 * embedment_ratio if v > 0.5 else -0.09 - 0.88 * embedment_ratio

    # === Failure function ===
    vals = (np.abs(h_grid / h_star) ** alpha
            + (m_grid / m_star) ** alpha
            + 2 * beta * ((h_grid * m_grid) / (h_star * m_star)))

    # === Extract points near failure surface ===
    tolerance = 0.001
    mask = np.abs(vals - 1) < tolerance
    h_vals = h_grid[mask]/h_star
    m_vals = m_grid[mask]/m_star

    return h_vals, m_vals
+