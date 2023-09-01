# mRKC
WORK IN PROGRESS

### The modified equation
bla

### Combination of RKC methods
bla 


### References
This code has been used for the numerical experiments presented in:

>Abdulle, A., Grote, M. J., & Rosilho de Souza, G. (2022). _Explicit stabilized multirate method for stiff differential equations_. Mathematics of Computation, 91(338), 2681–2714.


## Install and Run
The easiest way to run the code is by creating a Docker image from the Dockerfile provided here and running the code within a container. Otherwise, one could install the dependencies and compile the code as usual.

### Docker
To create the Docker image, from the root directory of the repository, run:
> docker build -t mRKC .  

For running the code, execute:
> docker run --rm -ti -v "$(pwd)/results":/mRKC/results mRKC ARGS_LIST

For the `ARGS_LIST`, see below.

### Compile
Before compiling the code install the dependencies: [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [gmsh](https://gmsh.info/), and [fftw3](https://www.fftw.org/). Then, compile[^comp_err] the code running:
> mkdir build && cd build && cmake -S .. -B . && make

Run the code from the ``build`` directory as
> ./main ARGS_LIST

For the `ARGS_LIST`, see the next section.

## Command line arguments
For the full list of arguments use the `--help` option or check out the [`ARGS_LIST.md`](ARGS_LIST.md) file, here we just provide an example.

To run a simulation

- with step size $\Delta t=0.01 ms$ and mesh size $\Delta x=10\mu m=1e\text{-}3cm$: `-dt 1e-2 -dG 1e-3`,
- for a duration of $10ms$: `-tend 10`,
- with output in `.vtk` format every 10 time steps: `-vtk true -ofreq 10`,
- on a block of 2x10 cells: `-nx 10 -ny 2`,
- with cells of length $c_l=100\mu m=0.01 cm$ and width $c_w=20\mu m=0.002 cm$: `-cl 1e-2 -cw 2e-3`,
- both leftmost cells are stimulated: `-nic 2`,
- with vertical gap junctions having a wave shape with amplitude $1\mu m=1e\text{-}4 cm$, frequency $2$ and being smooth: `-va 1e-4 -vf 2 -vs true`,
- with horizontal gap junctions having a wave shape with amplitude $2\mu m=2e\text{-}4 cm$, length $30\mu m=0.003 cm$ and being squared: `-ha 2e-4 -hl 3e-3 -hs false`,
- with horizontal gap junctions randomly placed along the horizontal cell side and their permeability being zero with probability 30% : `-halt 2 -hprob 0.3`,

execute:

```
./bemi -dt 1e-2 -dG 1e-3 -tend 10 -vtk true -ofreq 10 -nx 10 -ny 2 -cl 1e-2 
-cw 2e-3 -nic 2 -va 1e-4 -vf 2 -vs true -ha 2e-4 -hl 3e-3 -hs false -halt 2 
-hprob 0.3
```

The output is stored in the `results` folder.

# Acknowledgements
Assyr???


<p align="center">
  <img src="./docs/img/logo-MICROCARD.png" height="80"/>
  <img src="./docs/img/EuroHPC.jpg" height="80" />
  <img src="./docs/img/WBF_SBFI_EU_Frameworkprogramme_E_RGB_pos_quer.png" height="80" />
</p>


[^comp_err]: If ``cmake`` doesn't  find the ``fftw3`` library, we suggest to compile ``fftw3`` using ``cmake`` (instead of the traditional ``./configure``). Check out the `Dockerfile` to see the commands that we used. 

