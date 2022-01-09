This code implements the TSMOR algorithm by Nair and Balajewicz, 2019 on 2D steady flow over forward facing step.

To use TSMOR,
1. cd into src directory
2. run main.m file 

Snapshots and transport fields are precomputed in snapshots and transport_fields directories respectively.

You can either load these transport fields or you can choose to have them computed. Set the following in line 24 of main.m:
load_saved_transport = 1 : to load saved transport fields (default)
load_saved_transport = 0 : to evaluate the tranports

Transport fields corresponding to 1 to 9 Fourier coefficients are available. Set the following in line 26 of main.m:
nfc = 9 (default) : Number of Fourier coefficients
