# Installation

falcon is a Finite Element Analysis software built on top of the [Jem-Jive](https://software.dynaflow.com/jive/) libraries. As such, one needs to obtain and compile them. There are two ways to go about:
- **1. Local install** (jem-jive dependencies installed on root, user priviledges required)
- **2. Containerization** (jem-jive and its depedencies within container, can be used on clusters)

## 1. Local install

### Install pre-requisites
On terminal, execute:
```sh
sudo apt-get install git build-essential g++ zlib1g-dev libreadline-dev freeglut3-dev
```

### Obtain Jem-Jive
Obtain Jem and Jive either from the [Dynaflow website](https://software.dynaflow.com/jive/) or from the [Unofficial Github repo](https://github.com/ritukeshbharali/jemjive-3.0). Following the second approach,  on terminal, execute:
```sh
cd /home/user
git clone https://github.com/ritukeshbharali/jemjive-3.0.git
```

### Build Jem
- Enter the jem directory
- For sequential version, on terminal, execute:
```sh
./configure
make -j 4
export JEMDIR='/home/user/libraries/jemjive/jem-3.0’>> ~/.bashrc
```

- For MPI version, on terminal, execute:
```sh
./configure --cxx=mpic++
make -j 4
export JEMDIR='/home/user/libraries/jemjive/jem-3.0’>> ~/.bashrc
```

### Build Jive
- Enter the jive directory
- On terminal, run:
```sh
./configure
make -j 4
export JIVEDIR='/home/user/libraries/jemjive/jive-3.0’>> ~/.bashrc
```

### Build falcon
Assuming falcon is downloaded and extracted, go to 'build' directory. On the terminal, execute:
```sh
make opt -j5
```

### Tests
Go to 'tests' directory, which contain some input files to check the installation. On the terminal , execute:
```sh
./do_tests.sh
```

### Demos
The directory 'demo' contains some demo problems (larger than those in tests). Some of them are benchmark problems. They can be run as:
```sh
/path-to-falcon/build/falcon-opt problem.pro
```

## 2. Containerization

Through containerization, all dependencies of falcon (including Jem and Jive builds) are packaged into a single container(.sif file). This way, one can use same environment both, on a local machine as well as on an HPC clusters. Here are the steps:

### Install apptainer
Add software-properties-common to be able to use to add-apt-repository command
```sh
sudo apt update
sudo apt install -y software-properties-common
```

Add apptainer repository and install
```sh
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

### Build the container
Enter 'apptainer' directory. On terminal, execute:
```sh
apptainer build jive30.sif jive30.def
```

### Build falcon
Once apptainer build the container image 'jive30.sif', it can be used to build falcon, both on a local machine as well as on HPC clusters. Here is how it is done:
- Enter the 'build' directory. On terminal, execute:
```sh
cp ALL_MAKEFILES/makefile.apptainer ./makefile
apptainer exec /path-to-apptainer/jive30.sif make opt -j4
```

### Tests
Go to 'tests' directory, which contain some input files to check the installation. On the terminal , execute:
```sh
apptainer exec /path-to-apptainer/jive30.sif ./do_tests.sh
```

### Running falcon simulations on local machine
- Enter the directory that contains a problem to be run (i.e., contains a .pro file, say problem.pro).
- For a serial (multithreaded) run, on terminal, execute:
```sh
apptainer exec /path-to-apptainer/jive30.sif /path-to-falcon/build/falcon-opt problem.pro
```

- For a MPI run with, say 4 procs, on terminal, execute:
```sh
apptainer exec /path-to-apptainer/jive30.sif mpirun -np 4 /path-to-falcon/build/falcon-opt problem.pro
```

### Running falcon simulations on HPC clusters
Most HPC clusters typically has apptainer or singularity installed. singularity is the old name of apptainer. The 'jive30.sif' must work with both apptainer and singularity! To build and run falcon simulations on the cluster, here is the way to go:
- **Build the container** jive30.sif on the local machine. Steps shown above.
- Copy to container from local machine to the cluster.
- Copy falcon directory to the cluster
- To build falcon on cluster, execute on terminal:
```sh
cp ALL_MAKEFILES/makefile.apptainer ./makefile
apptainer exec /path-to-apptainer/jive30.sif make opt -j4
```
- To run simulations, add the following lines to the job script. For a serial (multithreaded):
```sh
apptainer exec /path-to-apptainer/jive30.sif /path-to-falcon/build/falcon-opt problem.pro
```

- For a MPI run with, say 4 procs, on terminal, execute:
```sh
apptainer exec /path-to-apptainer/jive30.sif mpirun -np 4 /path-to-falcon/build/falcon-opt problem.pro
```

**Note**: If the local machine/cluster does not have apptainer but singularity, replace 'apptainer exec' with 'singularity exec'