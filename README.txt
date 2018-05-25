Philip Tully 2011

--First, to run NEST with MPI configured:

export LD_LIBRARY_PATH=/lib/:/usr/local/lib/nest/
lamboot

--Second, to add this module ('pt_module') to nest.Models():
sudo rm -r bootstrap-module-100725/
sudo rm -r build-module-100725/
sudo ./compile-module.sh module-100725/

--then, 
ipython
import nest
nest.Models()                     // 'bcpnn_synapse', 'iaf_cond_alpha_bias' should NOT show up
nest.Install('pt_module')
nest.Models()                     // now, 'bcpnn_synapse' is available

--use the new synapse model!




Notes by Florian Fiebig 2017

--First, install NEST 2.2.2 and Python 2.7 (or check its installation on the Supercomputer)
Check the start of the compile script for correct loading of these modules and compiler swap(might not work without manual typing sometime). 
Check the bottom of the compile script for the correct nest install path

--Second, to add this module ('pt_module') to nest.Models():
sudo rm -r bootstrap-module-100725/
sudo rm -r build-module-100725/
sudo ./compile-module.sh module-100725/

--then, in python
import nest
nest.sr('(/cfs/klemming/nobackup/f/fiebig/170501_Beskow_BCPNN/share/nest/sli) addpath')
nest.Install('/cfs/klemming/nobackup/f/fiebig/170501_Beskow_BCPNN/lib/nest/pt_module')
print nest.Models('synapses')                   // now, 'bcpnn_synapse' is available

--use the new synapse model!
