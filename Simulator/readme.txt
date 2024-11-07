This folder contains the simulator used to run the experiments.

Before running any experiment, it is mandatory to create the makefile in the <build> sub-folder. To do this, the user can type the following command
on the terminal:

qmake-qt4 cforaging.pro

Since the Qt4 libraries are deprecated and have been removed from Ubuntu 20.04, the user is referred to the following link:

https://ubuntuhandbook.org/index.php/2020/07/install-qt4-ubuntu-20-04/

Alternatively, the user can install a virtual machine with Ubuntu18.04 (or older versions).
Once the makefile has been created, the user must type the command "make" on the terminal. The user can ignore all the warning messages. If no errors
occur, the user should find an executable file, called "cforaging", in the <build> sub-folder. Instead, if errors are raised, please check all the paths
in the cforaging.pro file and repeat the steps described above.

To simplify the possibility to replicate the experiments, we provide the script launcher.sh, which requires two mandatory parameters:

1) configuration file;
2) number of replications.

There are three optional parameters:

3) seed to initialize the random number generator (default is 1);
4) directory containing the *.sample (default is current directory);
5) directory where the output files will be stored (default is current directory).

A typical usage is to copy the launcher.sh in the <cforaging> sub-folder and run the command from that sub-folder.
To run a simulation, you can type the following command on the terminal:

./launcher.sh cforaging.ini 10 1 ./ files/

The command will perform 10 replications of the experiment (starting with the seed 1). The experimental parameters are loaded by the cforaging.ini
configuration file. The *.sample files are supposed to be in the current directory. If this is not the case, an error is raised. Finally, all the
output files will be stored in the files sub-folder (which must exist within the current directory).
