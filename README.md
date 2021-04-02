Globally Optimal Solution to Inverse Kinematics of 7DOF Serial Manipulator
=====
Solver for inverse kinematics tasks for general serial manipulators with seven revolute joints.
It combines symbolical and numerical approaches to find a global solution to a polynomial objective function.
A toolbox for Matlab.

Usage
-----
Solver for each of the manipulators is stored in separate folders.
DH parameters of the manipulator are saved in the MATLAB file `manipulator.mat`.
For direct solving by the POP solver, please, refer to the script file `gloptipoly.m`.
For the approach with the symbolic reduction step, please, refer to the script file `GB.m`.
These files are self-explanatory.

For a detailed description of this software, please refer to [[1]](#1).

Several MATLAB toolboxes are required to be installed: GloptiPoly [[2]](#2), YALMIP [[3]](#3), and MOSEK [[4]](#4).
For symbolic preprocessing, Maple installation is needed.

Citing
-----
If you are using the software for (scientific) publications, please cite the following source:
```bibtex
@misc{softwareTrutmanGlobal,
  title = {{Global 7DOF IKT -- Globally Optimal Solution to Inverse Kinematics of 7DOF Serial Manipulator}},
  author = {Trutman, Pavel},
  howpublished = {\url{https://github.com/PavelTrutman/Global-7DOF-IKT}},
  year = {2021}
}
```
Please also cite the relevant publication:
```bibtex
@article{trutman2020globally,
  title = {{Globally Optimal Solution to Inverse Kinematics of 7DOF Serial Manipulator}},
  author = {Trutman, Pavel and Mohab, Safey El Din and Henrion, Didier and Pajdla, Tomas},
  journal = {arXiv preprint arXiv:2007.12550},
  year = {2020}
}
```

References
-----
<a id="1">[1]</a>
Pavel Trutman, Mohab Safey El Din, Didier Henrion, and Tomas Pajdla.
Globally Optimal Solution to Inverse Kinematics of 7DOF Serial Manipulator.
arXiv preprint arXiv:2007.12550, 2020.

<a id="2">[2]</a>
Didier Henrion, and Jean-Bernard Lasserre.
GloptiPoly: Global optimization over polynomials with Matlab and SeDuMi.
ACM Transactions on Mathematical Software (TOMS) 29.2: 165-194, 2003.

<a id="3">[3]</a>
Johan Löfberg.
YALMIP: A toolbox for modeling and optimization in MATLAB.
In Proceedings of the CACSD Conference, Taipei, Taiwan, 2004.

<a id="4">[4]</a>
MOSEK ApS.
The MOSEK optimization toolbox for MATLAB manual.
Version 8.0, 2016.
http://docs.mosek.com/8.0/toolbox/index.html
