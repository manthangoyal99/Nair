1.  Make a 'run' folder that is distinct from the code folder, say at 
    ~/Documents/Research/RUNS/SteadyROM/Euler1DSteady
2.  'cd' to this folder
3.  Generate the training database snapshots with the command:
        "python path_to_Euler1DSteady_PY/Quasi1D_Steady_Snapshots_Driver.py".
    This should create 'training_snapshots.npz' in the current folder. You will
    also get to visualize the results.
4.  Generate the testing (i.e., validation) database snapshots with the command:
        "python path_to_Euler1DSteady_PY/Quasi1D_Steady_Snapshots_Driver.py -v".
    This should create 'testing_snapshots.npz' in the current folder. You will
    also get to visualize the results.
5.  Run the TSMOR with the command:
        "python path_to_Euler1DSteady_PY/Quasi1D_Steady_TSMOR_Driver.py -n 4 -t".
    This instructs the program to use 4 basis functions for the grid distortion,
    and to do offline training in addition to the online ROM run. The offline
    training should create 'transport_fields_nf4.npz' in the current folder. You
    will also get to visualize the results.