import numpy as np
import nest

t_sim = 1000.
f_in = 100.

on_beskow = True
if (not 'bcpnn_synapse' in nest.Models('synapses')):
    if on_beskow:
	nest.sr('(/cfs/klemming/nobackup/f/fiebig/170501_Beskow_BCPNN/share/nest/sli) addpath')
	nest.Install('/cfs/klemming/nobackup/f/fiebig/170501_Beskow_BCPNN/lib/nest/pt_module')
    else:
	nest.Install('pt_module')

initial_weight = np.log(nest.GetDefaults('bcpnn_synapse')['p_ij']/(nest.GetDefaults('bcpnn_synapse')['p_i']*nest.GetDefaults('bcpnn_synapse')['p_j']))
initial_bias = np.log(nest.GetDefaults('bcpnn_synapse')['p_j'])
syn_param = {"weight": initial_weight, "bias": initial_bias,"K":1.0,"delay":1.0,"tau_i":10.0,"tau_j":10.0,"tau_e":100.0,"tau_p":1000.0}

neuron1_spike_gen = nest.Create('poisson_generator', params={'rate': f_in})
neuron2_spike_gen = nest.Create('poisson_generator', params={'rate': f_in})

neuron1 = nest.Create("iaf_neuron")
neuron2 = nest.Create("iaf_neuron")

nest.Connect(neuron1_spike_gen, neuron1, params={'weight':  1000.0})
nest.Connect(neuron2_spike_gen, neuron2, params={'weight':  1000.0})
nest.Connect(neuron1,neuron2,params=syn_param, model="bcpnn_synapse")
nest.Simulate(t_sim)
