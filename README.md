# PyCylinder Flow

![Static Badge](https://img.shields.io/badge/python-blue?logo=python&logoColor=yellow)
![GitHub](https://img.shields.io/github/license/SkirOwen/pycylinder_flow?color=green)

This code is a python translation of the MATLAB code of Jonas Latt; 
see section `Original Code` for details.

## Installation
Install the requirements from `requirements.txt`.  
It is recommended to use a virtual environment for that.

## Usage
`cylinder.py` runs the flow field over a grid of size `(lx, ly)` 
(preset as quite course in the file and you probably want to refine).  
To reach a stable regime `max_t > 2000`.

It spits out a bunch of data files containing the horizontal `vz`
and vertical `vy` velocity on the grid for each time step.

You can then use `plot_cylinder.m` to visualise the flow in increasing time steps.
Note, here you only want to take the flow from timesteps 
2000+ since the early time is literally 'turning on' the upstream flow
(from zero) so there is some funky behaviour we don't want to consider there.
Seems like you have to manually copy over maxT and tPlot from whatever 
you selected in the previous case.
I'd check the visualisation to make sure you're in the 'uniform flow' 
regime and not the setup time region.

You can play with the obstacle size and position (obst_{x,y,r}),
and the Reynolds number to generate more data.

```shell
python run.py -s 100 50 -i 400
```


```bash
$ ffmpeg -framerate 60 -pattern_type glob -i '*.png' -s:v 1920x1080 -c:v libx264 -pix_fmt yuv420p out.mp4
```

```bash
ffmpeg -i out.mp4 -vf "fps=10,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 output.gif
```

## Original code
The original (unmodified) code in `MATLAB` is located in `./matlab_original`.


## License
The code is under the `GNU General Public License`, to continue under the same license as the
original `MATLAB` code.
The original license has been preserved where it was present in the original files.