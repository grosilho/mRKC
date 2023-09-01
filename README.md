# mRKC
This code is a C++ implementation of several explicit stabilized methods for solving ordinary differential equations
$$
y'=f(y),\qquad y(0)=y_0
$$
and multirate ordinary differential equations
$$
y'=f_F(y)+f_S(y),\qquad y(0)=y_0,
$$
where $f_F$ is a stiff but cheap term and $f_S$ is a mildly stiff but expensive term.

Explicit stabilized methods use an increased number of stages to increase stability, in contrast to standard methods that use more stages to increase accuracy. Due to this different strategy, the stability domain grows quadratically along the negative real axis and the methods have no step size restriction (despite being explicit).

In this code we implement the following explicit stabilized methods: [RKC](), [ROCK], [RKL], RKU[^RKU] for the non multirate problem $y'=f(y)'$.

### References
This code has been used for the numerical experiments presented in:

>Abdulle, A., Grote, M. J., & Rosilho de Souza, G. (2022). _Explicit stabilized multirate method for stiff differential equations_. Mathematics of Computation, 91(338), 2681–2714.


## Install and Run
One possible way to run the code is by creating a Docker image from the Dockerfile provided here and running the code within a container. Otherwise, one could compile the code as usual.

### Docker
To create the Docker image, from the root directory of the repository, run:
> docker build -t multirate .  

For running the code, execute:
> docker run --rm -ti -v "$(pwd)/results":/multirate/results multirate ARGS_LIST

For the `ARGS_LIST`, see below.

### Compile
For compilation, [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) is needed. However, if the code has been downloaded with `git clone ...` then a compatible version of `eigen` is cloned automatically by the makefile[^eigen]. Hence, the following command will download the dependency and compile the code:

> mkdir build && cd build && cmake -S .. -B . && make

Run the code from the ``build`` directory as
> ./MultirateIntegrators ARGS_LIST

For the `ARGS_LIST`, see the next section.

## Command line arguments
For the full list of arguments use the `--help` option or check out the [`ARGS_LIST.md`](ARGS_LIST.md) file, here we just provide an example.

To run a simulation
- solving the PDE Brusselator benchmark (problem 5): `-test 5`,
- with the mRKC method and step size $\Delta t=0.01$: `-rk mRKC -dt 1e-2`,
- with output in `.bin` format every 10 time steps: `-bin true -ofreq 10`,
execute:

```
./MultirateIntegrators -test 5 -rk mRKC -dt 1e-2 -bin true -ofreq 10
```

The output is stored in the `results` folder.


<p align="center">
  <img src="./docs/img/logo-MICROCARD.png" height="80"/>
  <img src="./docs/img/EuroHPC.jpg" height="80" />
  <img src="./docs/img/WBF_SBFI_EU_Frameworkprogramme_E_RGB_pos_quer.png" height="80" />
</p>

[^eigen]: If the code has been downloaded as a zip, then you should manually place a compatible version of eigen (e.g. 3.4.1) in the `/external/eigen`folder.

