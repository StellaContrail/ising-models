<h1> Magnetization Simulation ( 2D/3D Ising model / XY model ) </h1>
<p>
    金属における磁性をコンピュータシミュレーションにより再現します。
</p>
<p>
    Simulation of magnetization in ferromagnetic metals.
</p>

<h2> Description and Simulation Results </h2>
<h3> 2D Ising model </h3>
<div class="description">
    <p>
        Ising模型とはSpinがLattice上に均等に分布した系によって磁性をシミュレーションするモデルです。
        Spinは半整数(-1/2, +1/2)をとり、隣り合ったSpin同士の相互作用はそれぞれのSpinベクトルの内積に相互作用の強さを記述するJ因子、および逆向きになったときにエネルギーが最大になるようにマイナス符号がかけられた形で表現されます。
        <br>
        シミュレーションでは全体の磁化がexp(-beta E)の重みに従って分布するように、Monte Carlo法 (Metropolis–Hastingsアルゴリズム) に従ってスピンの向きを変えます。
        <br>
        後にも様々な模型を紹介しますが、特にこのイジング模型と呼ばれるモデルでは簡単な割に非常によく温度による磁性の相転移を再現し、常磁性体の性質をうまくシミュレーションするモデルとなります。
    </p>
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
        先述したIsing模型は二次元のLatticeに分布したSpinをシミュレーションするモデルでしたが、このプログラムでは三次元に拡張したLattice上でシミュレーションを行っています。
        三次元のIsing模型では厳密解が得られないことが知られており、シミュレーションの有用性を発揮させられる良い例となっています。
    </p>
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
        XY模型はこれまで紹介したIsing模型とは異なり、連続的な向きを持つスピンでLatticeが構成されているような模型です。
        連続的な向きを持つため位相欠陥を起こしうることが考えられると思いますが、実際に低温領域においてXY模型は離散的な循環を持つ量子渦を発生させ、自発的対称性の破れを誘起させます。
        このようなXY模型に特異的な相転移はBKT転移として既に理論面、そして実験面からも発見されています。
        <br>
        ここではスピンの分布についてはプロットしていませんが、スピンの情報は配列に確保してあるため、逐次データを保存することで相転移を観測することが可能です。
        シュミレーションではスピンの向きを(cosX, sinX)として表現し、スピンをフリップさせることは鏡像操作を行うことで実現しています。
    </p>
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
