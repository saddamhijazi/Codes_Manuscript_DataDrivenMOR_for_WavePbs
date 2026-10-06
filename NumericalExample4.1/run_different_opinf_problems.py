import subprocess
import numpy as np

r_vals = np.int32(np.linspace(10,100,10))

filename = "opinf_optimize_symMat_EarlyStopping.py"

for val in r_vals:

    with open(filename, "r") as file:
        lines = file.readlines()
    
    lines[58] = f"r = {val}\n"  # Line 59 is index 58

    with open(filename, "w") as file:
        file.writelines(lines)
    print(f'Inferring the reduced operators for the reduced dimension {val} \n')
    # Run the script
    subprocess.run(["python3", filename])
