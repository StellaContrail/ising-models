# Ising model / XY model
Fortran programs that simulate system consisted of spins using Ising model and XY model  
This repository includes three different types of models.  
* 2-dimensional Ising model  
Ising model is an ideal system which consists of only 1/2 type of spins.  
In real world, Fe(Iron) for instance, is the kind of system.  
Thus this simulation can be also defined as "Simulation of Paramagnetic Material"
![Energy of 2D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising2d/energy.png)   
The temperature dependence of the energy (2D Ising model)  
![Magnetization of 2D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising2d/mag_h0.png)  
The temperature dependance of the magnetization (2D Ising model)
* 3-dimensional Ising model  
The 2 dimensional version of ising model stated above is not suitable in a sense of realistic, because paramagnetic materials exist in the 3D world, not 2D world.  
This program runs same algorithm as 2D Ising model, but extended to 3 dimensions.
![Energy of 3D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising3d/energy.png)  
The temperature dependence of the energy (3D Ising model)  
![Magnetization of 3D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising3d/precise_mag.png)  
The temperature dependence of the magnetization (3D Ising model)
* XY model  
This model is an ideal system which consists of spins which can face in any directions.  
To realize the continous spin, all spins are represented as (cosx, sinx)  
![Energy of XY Model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising_xy/energy.png)  
The temperature dependence of energy (XY Model)  
![Magnetization of XY Model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising_xy/mag.png)  
The temperature dependence of magnetization (XY Model)
