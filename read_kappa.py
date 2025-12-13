import h5py
import numpy as np

mesh = [15, 15, 1]  # q-point mesh

kappa_file = f"kappa-m{mesh[0]}{mesh[1]}{mesh[2]}.hdf5"
print(f"==> Reading {kappa_file} ...")

with h5py.File(kappa_file, "r") as f:
    print("Data:")
    for key in f:
        print(f"  {key:25} : array-shape {f[key].shape}")

    temperatures = f["temperature"][:]
    kappa = f["kappa"][:]
    gamma = f["gamma"][:]
    freqs = f["frequency"][:]

print(freqs[:,0:3])