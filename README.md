# Ising model / XY model
FORTRANを用いてイジング模型をシミュレーションします。  
Simulating of Ising model using FORTRAN.  
This repository includes three different types of models.  
* 2-dimensional Ising model  
イジング模型は1/2 スピンのみで構成された系です。  
1/2スピンを持つ粒子が格子状に配置されおり、常にup/downのどちらの状態にあります。  
これらの粒子は隣り合った他の粒子と相互作用し、Monte Carlo法 (Metropolis–Hastingsアルゴリズム) に従ってスピンの向きを変えます。  
このイジング模型と呼ばれるモデルでは非常によく相転移を再現し、強磁性体の性質をシミュレートするモデルとなります。  
そのため、以下に挙げる複数のシミュレーションの結果を表す画像は強磁性体の性質に合うものとなっています。  
Ising model is an ideal system which consists of only 1/2 spins.  
Suppose multiple spins are arranged in a lattice, and those spins are always in one of the states of up or down.  
Each of the spins can interact with its neighbors, also can be flipped according to Monte Carlo method (Metropolis–Hastings algorithm).  
This model, called Ising model, is able to simulate the phase transition very well.  
In this program, the interaction factor J is set to be 1 so as to be model of ferromagnetism.  
Therefore, result pictures given below shows properties of ferromagnetism very well.  
The temperature dependence of the energy (2D Ising model)  
![Energy of 2D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising2d/energy.png)   
The temperature dependance of the magnetization (2D Ising model)
![Magnetization of 2D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising2d/mag_h0.png)  
* 3-dimensional Ising model  
The 2 dimensional version of ising model stated above is not suitable in a sense of realistic, since paramagnetic materials exist in the 3D world, not 2D world.  
This program runs same algorithm as 2D Ising model, but extended to 3 dimensions.  
The temperature dependence of the energy (3D Ising model)  
![Energy of 3D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising3d/energy.png)  
The temperature dependence of the magnetization (3D Ising model)
![Magnetization of 3D Ising model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising3d/precise_mag.png)  
* XY model  
This model is an ideal system which consists of spins which can face in any directions.  
To realize the continous spin, all spins are represented as (cosx, sinx)  
The temperature dependence of energy (XY Model)  
![Energy of XY Model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising_xy/energy.png)  
The temperature dependence of magnetization (XY Model)
![Magnetization of XY Model](https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising_xy/mag.png)  
