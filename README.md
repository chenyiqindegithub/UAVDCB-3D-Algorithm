# Optimizing 3D Charging Direction Set for Directional UAV Chargers with Dual-Conical Charging Beams in 3D-WRSNs

This repository contains the implementation code for the paper:

> **Optimizing 3D Charging Direction Set for Directional UAV Chargers with Dual-Conical Charging Beams in 3D-WRSNs**  
> Zhenguo Gao* (Senior Member, IEEE), Yiqin Chen, Hsiao-Chun Wu (Fellow, IEEE), Zhufang Kuang, Shih-Hau Fang (Senior Member, IEEE)

---

## Overview

This project implements a comprehensive framework for optimizing charging direction sets in 3D Wireless Rechargeable Sensor Networks (3D-WRSNs) using UAV-mounted directional chargers with dual-conical charging beams. The proposed **UAVDCB-3D** algorithm efficiently computes optimal charging directions while ensuring full network coverage.


---

## Repository Structure

```
.
├── README.md                           # This file
├── RuncMFEDS.m                         # Main MFEDS algorithm implementation
├── AngleSetManager.m                   # Helper class for angle range management
├── AngleRangeContainer.m               # Data structure for angle range queries
├── wrsn_config.m                       # Main simulation framework and configuration
│
├── Dir_dual/                           # Direction computation algorithms
│   ├── DirNode_DualCone.m              # Node-based direction method
│   ├── DirACC_DualCone.m               # ACC method
│   ├── DirGCC_DualCone_Sampled.m       # GCC method
│   └── DirPloyhedron_DualCone.m        # Fixed polyhedron-based method
│
├── sysLKH/                             # TSP solver implementations
│   ├── lkhTSP.m                        # LKH heuristic solver interface
│   ├── greedyTSP.m                     # Greedy TSP algorithm
│   └── acoTSP.m                        # Ant Colony Optimization TSP
│
└── LKH/                                # LKH executable and configuration files
    ├── LKH.exe                         # LKH solver executable
    └── TSPLIB/                         # Problem instances directory
```

---

## Requirements

- MATLAB R2022b or later
- Optimization Toolbox (for `linprog` solver)
- Statistics and Machine Learning Toolbox (for `pdist`, `squareform`)
- LKH-3 solver (included in `LKH/` directory)

---

