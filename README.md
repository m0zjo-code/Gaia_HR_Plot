# Gaia HR Plot
Example on how to create a hertzsprung russell diagram from a sample of the Gaia dataset 

This example will download a chunk of the Gaia dataset, and plot on a HR diagram
If the dataset (approx 90 MByte) has already been downloaded, the data will be loaded from disk

The following libs are required (install via pip):
- matplotlib
- pandas
- astroquery

![Example HR Plot](Hertzsprung-Russell.png "Example Plot")

This work has made use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.

TODO - need to properly scale the colourmap
