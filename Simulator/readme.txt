This folder contains the simulator used to run the experiments.

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
