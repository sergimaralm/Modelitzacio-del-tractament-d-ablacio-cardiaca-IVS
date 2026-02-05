# Cardiac Ablation Treatment Modeling (IVS)

This repository contains the implementation of various numerical methods to simulate temperature evolution within the Interventricular Septum (IVS) during cardiac ablation treatment.

The main objective is to solve the one-dimensional **Heat Equation** to model how thermal energy propagates through the tissue, comparing different temporal discretization schemes.

## Repository Contents

The project is structured into Python modules, each implementing a specific numerical algorithm:

| File | Description |
| :--- | :--- |
| `eulerexplicito.py` | Implementation of the **Explicit Euler** method. |
| `eulerimplicit.py` | Implementation of the **Implicit Euler** method (unconditionally stable). |
| `crank-nicholson.py` | Implementation of the **Crank-Nicolson** method (2nd order, unconditionally stable). |
| `soluciotemps.py` | Analysis of temporal evolution and comparisons. |
| `funcions.py` | Auxiliary library with shared functions (analytical solution, Jacobi method, etc.). |
| `figures/` | Folder where automatically generated plots are saved. |
| `Pràctica 1.pdf` | Theoretical documentation and project description. |

## Requirements and Installation

To run this project, you will need **Python 3.x** and standard scientific computing libraries.

1. Clone the repository:
   ```bash
   git clone [https://github.com/sergimaralm/Modelitzacio-del-tractament-d-ablacio-cardiaca-IVS.git](https://github.com/sergimaralm/Modelitzacio-del-tractament-d-ablacio-cardiaca-IVS.git)
   cd Modelitzacio-del-tractament-d-ablacio-cardiaca-IVS
