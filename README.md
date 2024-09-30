<h1> Magnetization Simulation ( 2D/3D Ising model / XY model ) </h1>
<p>
    Simulation of magnetization in ferromagnetic metals.
</p>

<h2> Description and Simulation Results </h2>
<h3> 2D Ising model </h3>
<div class="description">
    <p>
        Ising Model is an ideal system that consists of only 1/2 spins.
        The spins are arranged in equidistance from each other, making themselves lattice field.
        Since each spins interact with each other due to the magnetic moments, spin pairs pointing opposite directions tend to align their directions.
        <br>
        In the computer simulation, the Monte Carlo method (Metropolis-Hastings algorithm) is used for example to reproduce this "flipping" phenomenon.
        ( Sometimes the Rejection sampling method, Molecular dynamics method are used instead, but we will skip that since we don't use the algorithms here. )
        More specifically, this algorithm is used to calculate the "probability of flipping" of exp(-beta E) since the system obeys Gibbs-Boltzmann distribution.
        <br>
        We will introduce other models later, but this simple model reproduces the phase transition of magnetization very well, simulating the basic properties of ferromagnetic materials.
    </p>
</div>

<div class="result">
    <div class="fig">
        <figure>
            <img src="https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising2d/mag_h0.png" alt="Magnetization of 2D Ising model">
            <br>
            <figcaption>
                <b>
                    Temperature dependance of the magnetization (2D Ising model)
                </b>
            </figcaption>
        </figure>
        <p>
            You can see that around the critical point (1/kBTc = -0.43), the magnetization curve jumps which indicates the phase transition at Tc.
            The phase transition can be obtained by simulating the model in more than 2 dimension.
        </p>
    </div>
    <div class="fig">
        <figure>
            <img src="https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising2d/energy.png" alt="Energy of 2D Ising model">
            <br>
            <figcaption>
                <b>
                    Temperature dependence of the energy (2D Ising model)
                </b>
            </figcaption>
        </figure>
        <p>
            We already have an exact solution for 2D Ising model without external magnetic field, but not for the case with external magnetic field.
            The computer simulation can also be useful for calculating the approximated solution for such a case.
        </p>
    </div>
</div>

<br>

<h3> 3D Ising model </h3>
<div class="description">
    <p>
        As described above, the Ising model is a lattice system with 1/2 spins with an exact solution when there's no external magnetic field.
        The lattice system can be extended to a 3D lattice system, with no exact solution found ever.
        It is also a good example to show the usefulness of computer simulations in a sense of providing the approximated solution even when the exact solution doesn't exist.
    </p>
</div>

<div class="result">
    <div class="fig">
        <figure>
            <img src="https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising3d/precise_mag.png" alt="Magnetization of 3D Ising model">
            <br>
            <figcaption>
                <b>
                    Temperature dependence of the magnetization (3D Ising model)
                </b>
            </figcaption>
        </figure>
    </div>
    <div class="fig">
        <figure>
            <img src="https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising3d/energy.png" alt="Energy of 3D Ising model">
            <br>
            <figcaption>
                <b>
                    Temperature dependence of the energy (3D Ising model)  
                </b>
            </figcaption>
        </figure>
    </div>
</div>

<br>

<h3> XY model </h3>
<div class="description">
    <p>
        XY model is different from the Ising model, in a sense of spins having continuous directions.
        Since spins can have continuous directions, it is suggested the system can have topological solitons, creating vortices with discrete circulation.
        Such a phase transition is already well studied and named "Berezinskii-Kosterlitz-Thouless transition".
        <br>
        The continuous spin direction can be expressed as cosine and sine function, and flipping as reflection.
        We don't show the spin distribution here, but spin information is stored in an array in the program, so you can write them out to a file to check the topological soliton.
    </p>
</div>

<div class="result">
    <div class="fig">
        <figure>
            <img src="https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising_xy/mag.png" alt="Magnetization of XY Model">
            <br>
            <figcaption>
                <b>
                    Temperature dependence of magnetization (XY Model)
                </b>
            </figcaption>
            <p>
                The BKT transition is successfully reproduced in the low-temperature region.
            </p>
        </figure>
    </div>
    <div class="fig">
        <figure>
            <img src="https://raw.githubusercontent.com/StellaContrail/IsingModel/master/ising_xy/energy.png" alt="Energy of XY Model">
            <br>
            <figcaption>
                <b>
                    The temperature dependence of energy (XY Model)   
                </b>
            </figcaption>
        </figure>
    </div>
</div>
