# voronoi-crystals
Python package to create polycrystals by voronoi tesselation

(The code actually performing the voronoi tesselation has been moved from ComputatonalPhysics/lammps-utilities to here 29th June 2019.)

## Installation 
The troubles of making python work with pyvoro and glumpy (which this package depends upon) are not addressed here, so the command below only works if all required packages are already installed. 
```
pip install git+https://github.com/henriasv/voronoi-crystals
```

## Usage example 
```
import numpy as np
from voronoi_crystals import Polycrystal, run_example, glumpyVisualize

L = 100
box = np.asarray([[0, L], [0, L], [0, L]])
myPolycrystal = Polycrystal(N=2, substance="methane_hydrate_SI_mw", box=box)
glumpyVisualize(myPolycrystal)
myPolycrystal.dumpPositionsLammpsData("hydrate_polycrystal.data")
```

This script will 
* create a polycrystal 
* visualize the polycrystal  
* dump the positions of tha particles in tthe polycrystal in a file called ```hydrate_polycrystal.data``` in the current working directory. 

### Example lammps script to minimize the polycrystal 
```
units		real
dimension 	3
boundary 	p p p
atom_style 	atomic
pair_style 	sw

timestep 10.0

read_data test.data

pair_coeff * * water_methane_hydrate.sw O C

group water 	type 	1
group methane 	type 	2


minimize 1e-4 1e-4 1000 1000

write_data minimized.data
```


### Example lammps script to relax and then deform the polycrystal 
```
units		real
dimension 	3
boundary 	p p p
atom_style 	atomic
pair_style 	sw

timestep 10.0

read_data minimized.data

variable confining_pressure equal 100
variable production_temperature equal 263.15
variable t_damp equal 2000
variable p_damp equal 20000


pair_coeff * * water_methane_hydrate.sw O C

group water 	type 	1
group methane 	type 	2

neigh_modify delay 0 every 1 check yes

change_box all triclinic

fix npt all npt temp ${production_temperature} ${production_temperature} ${t_damp} x ${confining_pressure} ${confining_pressure} ${p_damp} y ${confining_pressure} ${confining_pressure} ${p_damp} z ${confining_pressure} ${confining_pressure} ${p_damp}

velocity all create $(1.8*v_production_temperature) 35134623 loop geom

thermo 10
thermo_style custom step time temp press etotal pxx pyy pzz pxy pxz pyz lx ly lz xy xz yz 

run 100 # to inspect performance
thermo 100
run 1000 

fix deform all deform 1 xz delta 100
run 500000

```

