import numpy as np
from math import floor, ceil

from utils import dotdict, split_in_groups


# Modifiable params

is_verbose = [False, True][0]

is_msn_active = [False, True][1]
is_snr_active = [False, True][1]
is_gpe_active = [False, True, 30.0][0]
is_gpe_a_active = [False, True][0]
is_gpe_b_active = [False, True][0]
is_gpe_c_active = [False, True][0]


# Definition of parameters

params = dotdict()


# SIMULATION params ############################################################

params.sim = dotdict()
params.sim.step = 0.01 #Humphries 2009
params.sim.n_threads = 6
params.sim.verbose = is_verbose


# NETWORK params ###############################################################

params.net = dotdict()
params.net.p_channels = [0.4, 0.4, 0.2] #Tomkins 2014, 3 canales de 40%, 40% y 20% de neuronas, pag 5
params.net.ctx_afferents_per_neuron = 250.0 #Humphries 2010

params.net.n_msn_per_mm3 = 84900 #Humphries 2010, pag 2
params.net.box_side_mm = 0.3 #Tomkins 2014, pag 4
params.net.box_side_um = params.net.box_side_mm * 1000.0
params.net.box_volume_mm3 = params.net.box_side_mm ** 3

params.net.msn = dotdict()
params.net.msn.n = int(params.net.box_volume_mm3 * params.net.n_msn_per_mm3)
params.net.msn.d1_n = int(params.net.msn.n * 0.5)
params.net.msn.d2_n = int(params.net.msn.n * 0.5)
params.net.msn.is_active = is_msn_active
params.net.msn.reduction_factor = 1.0 / 10.0 #Cuando se usa una población poisson en lugar del msn, reducimos el número de neuronas para acelerar los calculos.

params.net.gpe = dotdict()
params.net.gpe.is_active = is_gpe_active
params.net.gpe.n = ceil(46000.0/300.0) #Fountas
params.net.gpe.has_local_collaterals = True

params.net.gpe_a = dotdict()
params.net.gpe_a.is_active = is_gpe_a_active
params.net.gpe_a.n = ceil(params.net.gpe.n * 0.0405)
params.net.gpe_b = dotdict()
params.net.gpe_b.is_active = is_gpe_b_active
params.net.gpe_b.n = ceil(params.net.gpe.n * 0.85)
params.net.gpe_c = dotdict()
params.net.gpe_c.is_active = is_gpe_c_active
params.net.gpe_c.n = ceil(params.net.gpe.n * 0.1095)

params.net.snr = dotdict()
params.net.snr.n = 3000 #int(26300/300) #Fountas, Oorschot 1996
params.net.snr.is_active = is_snr_active

# NEURONS params ###############################################################

params.neuron = dotdict()
params.neuron.msn_d1 = dotdict()
params.neuron.msn_d1.C_m = 15.0
params.neuron.msn_d1.E_rev_AMPA = 0.0 
params.neuron.msn_d1.E_rev_NMDA = 0.0
params.neuron.msn_d1.E_rev_GABA_A = -60.0
params.neuron.msn_d1.tau_decay_AMPA = 6.0 
params.neuron.msn_d1.tau_decay_NMDA = 160
params.neuron.msn_d1.tau_decay_GABA_A = 4.0
params.neuron.msn_d1.Kplus = 1.0
params.neuron.msn_d1.V_T = -30.0
params.neuron.msn_d1.vp = -80.0 #-76.5
params.neuron.msn_d1.V_th = 40.0
params.neuron.msn_d1.a = 0.01
params.neuron.msn_d1.b = -20.0
params.neuron.msn_d1.c = -55.0
params.neuron.msn_d1.d = 91.0
params.neuron.msn_d1.y1 = 0.3
params.neuron.msn_d1.y2 = 0.3
params.neuron.msn_d1.alpha_1 = 6.3 #Humphries
params.neuron.msn_d1.alpha_2 = 0.215 #Humphries
params.neuron.msn_d1.alpha = 0.0 #D1 type
params.neuron.msn_d1.c_1 = 0.0289 #D1 type
params.neuron.msn_d1.c_2 = 0.331 #D1 type

params.neuron.msn_d2 = dotdict()
params.neuron.msn_d2.C_m = 15.0
params.neuron.msn_d2.E_rev_AMPA = 0.0 
params.neuron.msn_d2.E_rev_NMDA = 0.0
params.neuron.msn_d2.E_rev_GABA_A = -60.0
params.neuron.msn_d2.tau_decay_AMPA = 6.0 
params.neuron.msn_d2.tau_decay_NMDA = 160.0
params.neuron.msn_d2.tau_decay_GABA_A = 11.0
params.neuron.msn_d2.Kplus = 1.0
params.neuron.msn_d2.V_T = -30.0
params.neuron.msn_d2.vp = -80.0 #-76.5
params.neuron.msn_d2.V_th = 40.0
params.neuron.msn_d2.a = 0.01
params.neuron.msn_d2.b = -20.0
params.neuron.msn_d2.c = -55.0
params.neuron.msn_d2.d = 91.0
params.neuron.msn_d2.y1 = 0.3 #phi1 in the repo, the proportion of active D1 receptors
params.neuron.msn_d2.y2 = 0.3 #phi2 in the repo, the proportion of active D2 receptors
params.neuron.msn_d2.alpha_1 = 6.3 #Humphries beta1 in the repo, synaptic dopamine, scaling coefficients determining the relationship between dopamine receptor occupancy and the effect magnitude
params.neuron.msn_d2.alpha_2 = 0.215 #Humphries beta2 in the repo
params.neuron.msn_d2.alpha = 0.032 #D2 type # alpha in the repo, lowers current values for D2
params.neuron.msn_d2.c_1 = 0.0 #D2 type K in the repo, potassium, enhancing the inward-rectifying potassium current
params.neuron.msn_d2.c_2 = 0.0 #D2 type L in the repo, L-type calcium current, lowering the activation threshold of the calcium current 

params.neuron.snr = dotdict()
params.neuron.snr.a = 0.113 #Fountas 
params.neuron.snr.b = 11.057 #Fountas
params.neuron.snr.c = -62.7 #Fountas
params.neuron.snr.d = 138.4 #Fountas
params.neuron.snr.V_r = -64.58 # -55.8, #Fountas
params.neuron.snr.V_T = -51.8 # -55.2, #Fountas
params.neuron.snr.V_th = 9.8 #20.0 # 9.8 #Fountas, peak
params.neuron.snr.C_m = 172.1 # 80.0, #Fountas #200.0 C_sim
params.neuron.snr.Kplus = 0.7836 # 1.731, #Fountas
params.neuron.snr.tau_decay_msn = 5.2
params.neuron.snr.tau_decay_snr = 3.0
params.neuron.snr.E_rev_msn = -80.0 
params.neuron.snr.E_rev_snr = -80.0
params.neuron.snr.I_e = 690.4 


# SYNAPSES params ###########################################################

params.syn = dotdict()

params.syn.ctx_msn = dotdict()
params.syn.ctx_msn.conn_spec = dotdict()
params.syn.ctx_msn.conn_spec.rule = 'one_to_one'
params.syn.ctx_msn.syn_spec_ampa = dotdict()
params.syn.ctx_msn.syn_spec_ampa.synapse_model = 'static_synapse'
params.syn.ctx_msn.syn_spec_ampa.delay = 10.0 #A mano, no importa
params.syn.ctx_msn.syn_spec_ampa.weight = 6.1 / params.neuron.msn_d1.tau_decay_AMPA #Humphries 2009b, pag 1178, tabla 2
params.syn.ctx_msn.syn_spec_ampa.receptor_type = 1
params.syn.ctx_msn.syn_spec_nmda = dotdict()
params.syn.ctx_msn.syn_spec_nmda.synapse_model = 'static_synapse'
params.syn.ctx_msn.syn_spec_nmda.delay = 10.0 #A mano, no importa
params.syn.ctx_msn.syn_spec_nmda.weight = 3.05 / params.neuron.msn_d1.tau_decay_NMDA #Humphries 2009b, pag 1178, tabla 2
params.syn.ctx_msn.syn_spec_nmda.receptor_type = 2


params.syn.msn_msn = dotdict()
params.syn.msn_msn.conn_spec_random = dotdict()
params.syn.msn_msn.conn_spec_random.rule = 'pairwise_bernoulli'
params.syn.msn_msn.conn_spec_random.p = 728.0 / params.net.msn.n #Tomkins 2014, pag 5
params.syn.msn_msn.conn_spec_random.allow_autapses = False
params.syn.msn_msn.syn_spec_random = dotdict()
params.syn.msn_msn.syn_spec_random.synapse_model = 'static_synapse'
params.syn.msn_msn.syn_spec_random.delay = {"distribution": "uniform", "low": 1.0, "high": 2.0} #Humphries, código
params.syn.msn_msn.syn_spec_random.weight = 0.25 #A mano
params.syn.msn_msn.syn_spec_random.receptor_type = 3

params.syn.msn_snr = dotdict()
params.syn.msn_snr.conn_spec = dotdict()
params.syn.msn_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.msn_snr.conn_spec.p = 0.033 #Fountas
params.syn.msn_snr.syn_spec_static = dotdict()
params.syn.msn_snr.syn_spec_static.synapse_model = 'static_synapse'
params.syn.msn_snr.syn_spec_static.delay = 1.0 #Fountas
params.syn.msn_snr.syn_spec_static.weight = 57.654700925145136 # Búsqueda local siguiendo procedimiento de Fountas #4.5 #Fountas
params.syn.msn_snr.syn_spec_static.receptor_type = 1

params.syn.msn_snr.syn_spec_plastic = dotdict()
params.syn.msn_snr.syn_spec_plastic.synapse_model = 'tsodyks2_synapse'
params.syn.msn_snr.syn_spec_plastic.U = 0.0192 #Fountas
params.syn.msn_snr.syn_spec_plastic.u = 1.0#0.0192 #Fountas
params.syn.msn_snr.syn_spec_plastic.x = 0.0192#1.0
params.syn.msn_snr.syn_spec_plastic.tau_rec = 623.0 #Fountas
params.syn.msn_snr.syn_spec_plastic.tau_fac = 559.0 #Fountas
params.syn.msn_snr.syn_spec_plastic.delay = 1.0 #Fountas
params.syn.msn_snr.syn_spec_plastic.weight = 156.3 #Fountas
params.syn.msn_snr.syn_spec_plastic.receptor_type = 1

params.syn.snr_snr = dotdict()
params.syn.snr_snr.conn_spec = dotdict()
params.syn.snr_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.snr_snr.conn_spec.p = 0.1 #Fountas
params.syn.snr_snr.syn_spec = dotdict()
params.syn.snr_snr.syn_spec.synapse_model = 'static_synapse'
params.syn.snr_snr.syn_spec.delay = 1.0 #Fountas
params.syn.snr_snr.syn_spec.weight = 0.32539679447089986 # Búsqueda local siguiendo procedimiento de Fountas #0.2 #Fountas
params.syn.snr_snr.syn_spec.receptor_type = 2

params.syn.gpe_gpe = dotdict()
params.syn.gpe_gpe.conn_spec = dotdict()
params.syn.gpe_gpe.conn_spec.rule = 'pairwise_bernoulli'
params.syn.gpe_gpe.conn_spec.p = 0.1

params.syn.gpe_gpe.syn_spec = dotdict()
params.syn.gpe_gpe.syn_spec.model = 'static_synapse'
params.syn.gpe_gpe.syn_spec.delay = 1.0
params.syn.gpe_gpe.syn_spec.weight = 0.765
params.syn.gpe_gpe.syn_spec.receptor_type = 4


params.syn.msn_gpe = dotdict()
params.syn.msn_gpe.conn_spec = dotdict()
params.syn.msn_gpe.conn_spec.rule = 'pairwise_bernoulli'
params.syn.msn_gpe.conn_spec.p = 0.033 #Fountas

params.syn.msn_gpe.syn_spec_plastic = dotdict()
params.syn.msn_gpe.syn_spec_plastic.model = 'tsodyks2_synapse'
params.syn.msn_gpe.syn_spec_plastic.U = 0.24 #Fountas
params.syn.msn_gpe.syn_spec_plastic.u = 1.0#0.24 #Fountas
params.syn.msn_gpe.syn_spec_plastic.x = 0.24#1.0
params.syn.msn_gpe.syn_spec_plastic.tau_rec = 11.0 #Fountas
params.syn.msn_gpe.syn_spec_plastic.tau_fac = 73.0 #Fountas
params.syn.msn_gpe.syn_spec_plastic.delay = 5.0 #Fountas
params.syn.msn_gpe.syn_spec_plastic.weight = 21.6 #Fountas
params.syn.msn_gpe.syn_spec_plastic.receptor_type = 3

params.syn.msn_gpe.syn_spec_static = dotdict()
params.syn.msn_gpe.syn_spec_static.model = 'static_synapse'
params.syn.msn_gpe.syn_spec_static.delay = 5.0 #Fountas
params.syn.msn_gpe.syn_spec_static.weight = 10.0 # Ajuste manual siguiendo procedimiento de Fountas #5.435 #Fountas
params.syn.msn_gpe.syn_spec_static.receptor_type = 3

params.syn.gpe_snr = dotdict()
params.syn.gpe_snr.conn_spec = dotdict()
params.syn.gpe_snr.conn_spec.rule = 'pairwise_bernoulli'
params.syn.gpe_snr.conn_spec.p = 0.1066

params.syn.gpe_snr.syn_spec_static = dotdict()
params.syn.gpe_snr.syn_spec_static.model = 'static_synapse'
params.syn.gpe_snr.syn_spec_static.delay = 3.0
params.syn.gpe_snr.syn_spec_static.weight = 59.67154181984671 #Búsqueda local siguiendo procedimiento de Fountas #93.0 #Fountas
params.syn.gpe_snr.syn_spec_static.receptor_type = 4

params.syn.gpe_snr.syn_spec_plastic = dotdict()
params.syn.gpe_snr.syn_spec_plastic.model = 'tsodyks2_synapse'
params.syn.gpe_snr.syn_spec_plastic.U = 0.196
params.syn.gpe_snr.syn_spec_plastic.u = 1.0#0.196
params.syn.gpe_snr.syn_spec_plastic.x = 0.196#1.0
params.syn.gpe_snr.syn_spec_plastic.tau_rec = 969.0
params.syn.gpe_snr.syn_spec_plastic.tau_fac = 0.0
params.syn.gpe_snr.syn_spec_plastic.delay = 3.0
params.syn.gpe_snr.syn_spec_plastic.weight = 603.9
params.syn.gpe_snr.syn_spec_plastic.receptor_type = 4


