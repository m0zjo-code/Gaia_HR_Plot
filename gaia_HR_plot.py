# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 23:33:33 2021
M0ZJO
"""

# Import processing libs
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.gaia import Gaia
import logging as log
from matplotlib import cm
import numpy as np

# Attempt to load data from disk
# If no data file is present, download it (May take a long time!)
# Note that this request is awaiting a job to be carried out on the remote server, don't worry if it looks like nothing is happening
try:
    log.info("Getting the DR2 results from file...")
    df = pd.read_csv("hr_data.csv")
except OSError:
    log.info("File not found, downloading...")
    job = Gaia.launch_job_async("select top 200000000"
        " bp_rp, phot_g_mean_mag+5*log10(parallax)-10 as mg, ra, dec" 
        " from gaiadr2.gaia_source"
        " where parallax_over_error > 10"
        " and visibility_periods_used > 8"
        " and phot_g_mean_flux_over_error > 50"
        " and phot_bp_mean_flux_over_error > 20"
        " and phot_rp_mean_flux_over_error > 20"
        " and phot_bp_rp_excess_factor <"
        " 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)"
        " and phot_bp_rp_excess_factor >"
        " 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)"
        " and astrometric_chi2_al/(astrometric_n_good_obs_al-5)<"
        "1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))"
        +" and 1000/parallax <= 200",
        dump_to_file=True,
	verbose=True,
        output_format="csv",
        output_file="hr_data.csv")
    print(job)
    r = job.get_results()
    df = r.to_pandas()

# Set up figure
fig = plt.figure(
    figsize=(14, 14),
    facecolor='black'
    )
ax = fig.add_axes([.1, .1, .85, .8])

ax.set_facecolor('black')
ax.set_title('Gaia Hertzsprung-Russell Diagram', color='white', fontsize=20)
ax.set_xlabel('Colour Index (B-V)', color='white', fontsize=14)
ax.set_ylabel('Magnitude (Mv)', color='white', fontsize=14)
ax.set_xlim(min(df['bp_rp']) - .1, max(df['bp_rp']) + .1)
ax.set_ylim(max(df['mg']) + 1.0, min(df['mg']) - 1.0)
ax.tick_params(top='off', right='off', direction='in', colors='white')

im = ax.scatter( df['bp_rp'], df['mg'],
    marker='.',
    s=[1] * len(df),
    c=np.array(df['bp_rp']),
    cmap = cm.jet,
    linewidth=0)
# Set scaled colourbar - TODO need to rescale this to match the visible spectrum
fig.colorbar(im, ax=ax)

# Save a higher resolution version as file
plt.savefig("Hertzsprung-Russell.png", facecolor='black', edgecolor='white', dpi=240)
# Display in GUI 
plt.show()
