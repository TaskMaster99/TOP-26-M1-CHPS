#!/usr/bin/env python3

import matplotlib.pyplot as plt

processes = [1, 2, 3, 4]
mlups = [13.71, 28.26, 54.89, 28.28]

ideal_mlups = [mlups[0] * p for p in processes]

plt.figure(figsize=(10, 6))

plt.plot(processes, mlups, marker='o', linestyle='-', linewidth=2, color='#007acc', label='Actual Performance (MLUPS)')
plt.plot(processes, ideal_mlups, linestyle='--', color='gray', alpha=0.7, label='Ideal Linear Scaling')

plt.title('LBM Simulation Scalability', fontsize=14)
plt.xlabel('Number of MPI Processes', fontsize=12)
plt.ylabel('FOM (Million Lattice Updates per Second)', fontsize=12)
plt.xticks(processes)
plt.grid(True, which='both', linestyle=':', alpha=0.5)
plt.legend()

plt.tight_layout()
plt.savefig("performance/plot/perf.png")