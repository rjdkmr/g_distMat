## g_distMat

### This tool is moved to <code>[gmx_clusterByFeatures](https://gmx-clusterbyfeatures.readthedocs.io/)</code> as sub-command <code>[distmat](https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/distmat.html)</code>.

**<code>[gmx_clusterByFeatures](https://gmx-clusterbyfeatures.readthedocs.io/)</code>** is easy to install by one command: **``sudo pip3 install gmx-clusterByFeatures``**. **<code>[distmat](https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/distmat.html)</code>** command is same as **``g_distMat``**. **HELP** and **USAGE** for **``distmat``** is provided in this **[link](https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/distmat.html)**. Moreover, **<code>[gmx_clusterByFeatures](https://gmx-clusterbyfeatures.readthedocs.io/en/latest/usage.html)</code>** include **[``matplot``](https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/matplot.html)**, which can be used to plot the data obtained from **<code>[distmat](https://gmx-clusterbyfeatures.readthedocs.io/en/latest/commands/distmat.html)</code>**.

***

### About

This tool is similar to gromacs tool g_mdmat. It calculates average minimum-distance matrix and other related matrices between two atom-groups from molecular dynamics trajectory (GROMACS, NAMD or AMBER).

This tool might be helpful to analyze binding interface in protein-protein, protein-ligand or other type of bio-molecular complexes.

##### Features

* Average minimum residues-distance matrix between two atom-groups that may be identical or different.
* Related Standard deviation and variance Matrix.
* Contact-map (fraction of contacts over entire trajectory) for the given distance cut-off.
* Parallel Calculation -- multi-threading to use more than one core of the multi-core processor.
* Compare fluctuations in two different trajectory


### Requirements
**Gromacs versions supported**: 4.5.x, 4.6.x, 5.0.x, 5.1.x and 2016.x

To compile and install, GROMACS libraries ``libgmx.a(.so)`` (v4.5.x/v4.6.x) or ``libgromacs.a(.so)`` (v5.0.x/v5.1.x/v2016.x) are required.


### Download
```bash
git clone https://github.com/rjdkmr/g_distMat
```


### Installation

```bash
cd g_distMat
mkdir build
cd build
cmake ..  -DGMX_PATH=/opt/gromacs -DCMAKE_INSTALL_PREFIX=/opt/g_distMat
make
make install
```

Directory <code>/opt/gromacs</code> should contains <code>include</code> and <code> lib </code> directories.

If fftw library <code> libfftw3f.so or libfftw3f.a </code> are not present in standard locations:
```bash
-DFFTW_LIB=/path/to/fftw3/lib
```

### Usage
```bash
g_distMat -h
```

### Description

It calculates average minimum-distance matrix of residues between two
atom-groups. Simultaneously, it calculates variance  and standard-deviation
matrices. Also, it calculates contact-frequency map over the trajectory for
the residues that are within a minimum distance given by "-ct" option value.

To speed up the calculation, it uses all available cores of the CPU using
multi-threading. Number of threads/cores could be change by "-nt" option.

To calculate distance variances or deviation (fluctuations) in a trajectory
with respect to average distances from another trajectory, use [-f traj_for_average.xtc]  
and [-f2 traj_for_variance.xtc]. The averages will be
calculated from "traj_for_average.xtc". Subsequently, variances and deviation
will be calculated for "traj_for_variance.xtc" with respect to previosly
calculated averages.

### Options

```

Option     Filename  Type         Description
------------------------------------------------------------
  -f       traj.xtc  Input        Trajectory: xtc trr trj gro g96 pdb cpt
  -s      topol.tpr  Input        Structure+mass(db): tpr tpb tpa gro g96 pdb
  -n      index.ndx  Input, Opt.  Index file
 -f2       traj.xtc  Input, Opt.  Trajectory: xtc trr trj gro g96 pdb cpt
-mean   average.dat  Output       Generic data file
-var   variance.dat  Output, Opt. Generic data file
-std stdeviation.dat  Output, Opt. Generic data file
-cmap contact_map.dat  Output, Opt. Generic data file

Option       Type   Value   Description
------------------------------------------------------
-[no]h       bool   yes     Print help info and quit
-[no]version bool   no      Print version info and quit
-nice        int    19      Set the nicelevel
-b           time   0       First frame (ps) to read from trajectory
-e           time   0       Last frame (ps) to read from trajectory
-dt          time   0       Only use frame when t MOD dt = first time (ps)
-ct          real   0.4     cut-off distance (nm) for contact map
-nt          int    8       number of threads for multi-threading

```

### Ouput
Output files contain values in two-dimensional matrix format. This type of file could be visualized with gnuplot as shown in following example:

```
gnuplot> plot 'contact_map.dat' matrix with image
```
Color of this plot could be easily changed as per the gnuplot manual.
