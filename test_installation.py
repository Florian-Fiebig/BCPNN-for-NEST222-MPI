import nest
print nest.Models()                     # 'bcpnn_synapse' should NOT show up

on_beskow = True
if (not 'bcpnn_synapse' in nest.Models('synapses')):
    if on_beskow:
	nest.sr('(/cfs/klemming/nobackup/f/fiebig/170501_Beskow_BCPNN/share/nest/sli) addpath')
	nest.Install('/cfs/klemming/nobackup/f/fiebig/170501_Beskow_BCPNN/lib/nest/pt_module')
    else:
	try:
		nest.Install('pt_module')
	except:
		nest.Install('pt_module')

status=('bcpnn_synapse' in nest.Models('synapses'))
print 'bcpnn now available', status
nest.GetDefaults('bcpnn_synapse')
