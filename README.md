This repository contains the code implementing the attacks against PROV described in the paper "Polynomial-Time Key-Recovery Attack on the PROV Specification", by River Moreira Ferreira and Ludovic Perret in 2024.
We recall than we takes as reference the implementation of Pierre Pébereau with the paper "One vector to rule them all: Key recovery from one vector in UOV schemes" in 2023. 


The code should be executed on a machine with sagemath 9.5 avalaible, with the following command:
```
$ sage Attack.sage
```

The attack can be performed on every security levels with parameter suggested in PROV specification. The security level can be choosed with the following command where n equals 1,3 or 5.
```
$ sage Attack.sage n
```
If the parameter n is ommited (or incorrect ), the attack is performed on the security level 1.