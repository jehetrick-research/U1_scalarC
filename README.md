## U(1) lattice gauge simulation in (scalar) C


* Update is Metropolis algorithm

* Usage:

```bash
$ pureU1 -h
pureU1 command line options (defaults):
        -N (8)  -Nx (8) -Nt (8)
         or: -L 8x8
        -beta (1.000000)
        -warms (0)
        -trajecs (1)
        -meas (1)
        -seed (1509432272)
        -init [cold, hot, or file filename] (cold)
```

Used this project to reproduce 
[DeGrand and Toussaint, Phys Rev. D22, 2478 (1980)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.22.2478).

Includes addition of flux via twisted BCs.

- [ ] **Need to add Jupyter Notebook summary and plots.**
