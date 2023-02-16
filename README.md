# Certifiable solvers for the relative pose problem between two calibrated cameras

Estimate the relative pose between two calibrated cameras. 
DIfferent formulations of the problem 
are used to certify the solution

Refer to our paper [HERE](https://mapir.isa.uma.es/mapirwebsite/?p=1718) 
for more information.



**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GNUv3](https://github.com/mergarsal/RedCertRPp/blob/main/LICENSE)


If you use this code for your research, please cite:

1. Iterative certifier

```
@ARTICLE{,
    author = {Garcia-Salguero, Mercedes and Gonzalez-Jimenez, Javier},
     title = {Fast certifiable relative pose estimation with gravity prior},
   journal = {Artificial Intelligence},
      year = {2023},
       url = {http://mapir.isa.uma.es/papersrepo/2023/2023_mercedes_AI_priorRPp_doc.pdf},
       doi = {https://doi.org/10.1016/j.artint.2023.103862}
}
```

2. Redundant formulations for the relative pose problem
```
@ARTICLE{,
    author = {Garcia-Salguero, Mercedes and Briales, Jesus and Gonzalez-Jimenez, Javier},
     month = {April},
     title = {A Tighter Relaxation for the Relative Pose Problem Between Cameras},
   journal = {Journal of Mathematical Imaging and Vision},
    volume = {},
    number = {},
      year = {2022},
       url = {https://link.springer.com/article/10.1007/s10851-022-01085-z},
       doi = {10.1007/s10851-022-01085-z},
     pages = {}
}
```


# Dependencies

The certifier requires : 
1. *Optimization* by D. Rosen. 
We use our fork

```
https://github.com/mergarsal/Optimization
```
2. *Iterative certifier* by us. 
```
https://github.com/mergarsal/QCQPIterCertifier
```

Use 
```
git submodule update --init --recursive
```
To download the dependency


## Build
```
git clone https://github.com/mergarsal/RedCertRPp.git
cd RedCertRPpPrior

mkdir build & cd build 

cmake .. 

make -jX

```

The compiled examples should be inside the `bin` directory. Run: 
```
        ./bin/example_certifier
```
 

