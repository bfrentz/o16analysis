# o16analysis

This is the analysis program to the 13C(a, n)16O -> 4a + n experiment that took place at the University of Notre Dame Nuclear Science Laboratory in Spring 2015.

## Running

This can only be run on the Notre Dame CRC (as it is the storage location for the experimental data). This also requires an installation of the root analysis software distributed by CERN. 

The commands to run the program (after starting root):
```
.L libExpEvent.so
.L analysis4.6.C
Run(int energy)
```

where energy is a user input for energy, corresponding to energy= 24, 25, 27, 28, or 29. It will aggregate all data from those energies according to the include files difining the different runs. 

It can also be run in batch mode using the scripts contained on the CRC. 

### Authorship

Frentz, Bryce<br/>
Gyurjinyan, Armen<br/>
Tan, Wanpeng<br/>
Sauer, Ethan

Updated January 2019
