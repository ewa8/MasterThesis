import nest
import numpy as np
from bg_params import params
from utils import verbose, split_in_groups, raster_plot
from copy import deepcopy
import matplotlib.pyplot as plt
# import utils


class BasalGanglia:

    def __init__(s, curr, _params=None):
        s.params = _params if _params is not None else params.copy()
        s.curr = curr #if external current needed

    def configure(s):
        #creating the msn, poisson generators and setting the number of striatal parrot neurons to be equal to one for each neuron
        s._create_nodes()
        #connecting each poisson generator to a parrot neuron, connecting parrot neuron to msn
        s._create_connections()
        s._create_detectors()

    def _create_nodes(s):
        # creating 1 poisson generator for each channel
        s.poisson = [nest.Create('poisson_generator', 1, {'rate': 0.0}) for _ in range(len(s.params.net.p_channels))]

        # adding Poisson generators for the STN and the GPE
        s.poisson_stn = nest.Create('poisson_generator', 1, {'rate': 0.0})
        s.poisson_gpe = [nest.Create('poisson_generator', 1, {'rate': 0.0}) for _ in range(len(s.params.net.p_channels))]
        

        s.str_channels = [[] for _ in s.params.net.p_channels]
        if s.params.net.msn.is_active:
            s._create_msn()

        if s.params.net.gpe.is_active:
            s._create_gpe()

        if s.params.net.snr.is_active:
            s._create_snr()

        #creating parrot neurons, number equaling the length of the number of neurons in that channel
        if s.params.net.msn.is_active:
            s.parrots_str = [nest.Create('parrot_neuron', len(channel)) for channel in s.str_channels]

    def _create_msn(s):
        # Queremos crear un sustituto de este núcleo que dispare a una frecuencia determinada
        if type(s.params.net.msn.is_active) is not bool:
            #creating a parrot neuron, half the amount of msn neurons, times the reduction factor for using the poisson population
            s.msn_d1 = nest.Create('parrot_neuron', int(s.params.net.msn.d1_n * s.params.net.msn.reduction_factor))

            #creating a parrot neuron, half the amount of msn neurons, times the reduction factor for using the poisson population
            s.msn_d2 = nest.Create('parrot_neuron', int(s.params.net.msn.d2_n * s.params.net.msn.reduction_factor))

            #adding the neuron populations
            s.msn = s.msn_d1 + s.msn_d2

            #poisson generator will generate a unique spike train for all of its targets, you have to use a parrot neuron in between the poisson generator and the target
            s.msn_poisson = nest.Create('poisson_generator', 1, {'rate': s.params.net.msn.is_active})
            nest.Connect(s.msn_poisson, s.msn)

            #returning a list of lists with neurons for each channel
            s.msn_d1_channels = split_in_groups(s.msn_d1, s.params.net.p_channels)
            s.msn_d2_channels = split_in_groups(s.msn_d2, s.params.net.p_channels)

            #aggregating the number of d1 and d2 neurons for each channel
            s.msn_channels = [d1+d2 for (d1,d2) in zip(s.msn_d1_channels, s.msn_d2_channels)]
            #creating the striatal channels to look like msn
            s.str_channels = deepcopy(s.msn_channels)

        # Queremos crear el núcleo normal
        else:
            #creating msn neurons
            s.msn_d1 = nest.Create('msn', s.params.net.msn.d1_n)
            s.msn_d2 = nest.Create('msn', s.params.net.msn.d2_n)
            nest.SetStatus(s.msn_d1, s.params.neuron.msn_d1)
            nest.SetStatus(s.msn_d2, s.params.neuron.msn_d2)

            #for if curve
            # nest.SetStatus(s.msn_d1, {'I_e': s.curr})

            #combining the 2 populations of neurons
            s.msn = s.msn_d1 + s.msn_d2

            #returning a list of lists with neurons for each channel 
            s.msn_d1_channels = [s.msn_d1[:int((len(s.msn_d1)*2/5))], s.msn_d1[int((len(s.msn_d1)*2/5)):int((len(s.msn_d1)*4/5))], s.msn_d1[int((len(s.msn_d1)*4/5)):]]
            s.msn_d2_channels = [s.msn_d2[:int((len(s.msn_d2)*2/5))], s.msn_d2[int((len(s.msn_d2)*2/5)):int((len(s.msn_d2)*4/5))], s.msn_d2[int((len(s.msn_d2)*4/5)):]]

            #aggregating the number of d1 and d2 neurons for each channel 
            s.msn_channels = [d1+d2 for (d1,d2) in zip(s.msn_d1_channels, s.msn_d2_channels)]

            #assigning for the striatal channels to be equal to msn channels
            for i in range(len(s.str_channels)):
                s.str_channels[i] += s.msn_channels[i]

            # Establezco valores iniciales de las neuronas. Pongo el voltaje al mínimo
            #putting the initial voltage to minimum
            dUVms = [{"V_m": s.params.neuron.msn_d1.vp} for _ in s.msn]
            nest.SetStatus(s.msn, dUVms)

    def _create_gpe(s):

        # Queremos crear un sustituto de este núcleo que dispare a una frecuencia determinada
        if type(s.params.net.gpe.is_active) is not bool:
            s.gpe_poisson = nest.Create('poisson_generator', 1, {'rate': s.params.net.gpe.is_active})
            s.gpe = nest.Create('parrot_neuron', s.params.net.gpe_a.n + s.params.net.gpe_b.n + s.params.net.gpe_c.n)
            nest.Connect(s.gpe_poisson, s.gpe)
            s.gpe_channels = split_in_groups(s.gpe, s.params.net.p_channels)

        # Queremos crear el núcleo normal
        else:
            s.gpe = nest.NodeCollection()
            s.gpe_channels = []

            if s.params.net.gpe_a.is_active:
                s.gpe_a = nest.Create("gpe", s.params.net.gpe_a.n)
                nest.SetStatus(s.gpe_a, s.params.neuron.gpe_a)
                s.gpe_a_channels = [s.gpe_a[:int((len(s.gpe_a)*2/5))], s.gpe_a[int((len(s.gpe_a)*2/5)):int((len(s.gpe_a)*4/5))], s.gpe_a[int((len(s.gpe_a)*4/5)):]]

                s.gpe += s.gpe_a
                if len(s.gpe_channels) != 3:
                    for from_channel in s.gpe_a_channels:
                        s.gpe_channels.append(from_channel)
                else:
                    s.gpe_channels = [to_channel+from_channel for (to_channel, from_channel) in zip(s.gpe_channels, s.gpe_a_channels)]

            if s.params.net.gpe_b.is_active:
                s.gpe_b = nest.Create("gpe", s.params.net.gpe_b.n)
                nest.SetStatus(s.gpe_b, s.params.neuron.gpe_b)
                s.gpe_b_channels = [s.gpe_b[:int((len(s.gpe_b)*2/5))], s.gpe_b[int((len(s.gpe_b)*2/5)):int((len(s.gpe_b)*4/5))], s.gpe_b[int((len(s.gpe_b)*4/5)):]]

                s.gpe += s.gpe_b
                if len(s.gpe_channels) != 3:
                    for from_channel in s.gpe_b_channels:
                        s.gpe_channels.append(from_channel)
                else:
                    s.gpe_channels = [to_channel+from_channel for (to_channel, from_channel) in zip(s.gpe_channels, s.gpe_b_channels)]


            if s.params.net.gpe_c.is_active:
                s.gpe_c = nest.Create("gpe", s.params.net.gpe_c.n)
                nest.SetStatus(s.gpe_c, s.params.neuron.gpe_c)
                s.gpe_c_channels = [s.gpe_c[:int((len(s.gpe_c)*2/5))], s.gpe_c[int((len(s.gpe_c)*2/5)):int((len(s.gpe_c)*4/5))], s.gpe_c[int((len(s.gpe_c)*4/5)):]]

                s.gpe += s.gpe_c
                if len(s.gpe_channels) != 3:
                    for from_channel in s.gpe_c_channels:
                        s.gpe_channels.append(from_channel)
                else:
                    s.gpe_channels = [to_channel+from_channel for (to_channel, from_channel) in zip(s.gpe_channels, s.gpe_c_channels)]


    def _create_snr(s):
        s.snr = nest.Create("snr", s.params.net.snr.n)
        nest.SetStatus(s.snr, s.params.neuron.snr)
        s.snr_channels = [s.snr[:int((len(s.snr)*2/5))], s.snr[int((len(s.snr)*2/5)):int((len(s.snr)*4/5))], s.snr[int((len(s.snr)*4/5)):]]
        
        s.parrots_gpe = [nest.Create('parrot_neuron', len(channel)) for channel in s.snr_channels]
        s.parrots_stn = nest.Create('parrot_neuron', 1)

    def _create_connections(s):
        channels = list(range(len(s.params.net.p_channels)))

        vt_sg = nest.Create("spike_generator",
                            params={"spike_times":[2300]})

        # vt_pg = nest.Create("poisson_generator", {'rate': 30})
        # create  volume transmitter
        vt = nest.Create("volume_transmitter")
        vt_parrot = nest.Create("parrot_neuron")
        nest.Connect(vt_sg, vt_parrot)
        nest.Connect(vt_parrot, vt, syn_spec={"synapse_model": "static_synapse",
                                            "weight": 1.,
                                            "delay": 1.})   # delay is ignored!
        # nest.Connect(vt_pg, vt_parrot)
        # nest.Connect(vt_parrot, vt, syn_spec={"synapse_model": "static_synapse",
        #                                     "weight": 1.,
        #                                     "delay": 1.})   # delay is ignored!
        vt_gid = vt.get("global_id")

        for i in channels:
            nest.Connect(s.poisson[i], s.parrots_str[i], conn_spec={'rule':'pairwise_bernoulli', 'p': 0.25}, syn_spec={"delay": 1.})
        #     # nest.Connect(post_sg, post_neuron, "one_to_one", syn_spec={"delay": 1., "weight": 9999.})
        #     # nest.Connect(s.parrots_str[i], s.msn_channels[i], "all_to_all", syn_spec={'synapse_model': 'stdp_nestml_rec'})
        #     # nest.Connect(s.parrots_str[i], s.msn_channels[i], "all_to_all", syn_spec={'synapse_model': 'stdp_nestml_rec'})



        # Cortex "afferents" (poisson)
        for i in channels:
            if s.params.net.msn.is_active: nest.Connect(s.poisson[i], s.parrots_str[i])

        # MSN efferents
        if s.params.net.msn.is_active:
            nest.CopyModel('stdp_dopamine_synapse', "stdp_nestml_rec",
            {"weight": 0.5,
             "delay": 1.0,
             "receptor_type": 0,
             "vt": vt_gid,
             "Wmax": 0.5,
             "Wmin":0.5,
             "tau_n": 20,
             "tau_c": 50,
            })
            
            if type(s.params.net.msn.is_active) is bool:
                # Connect parrots (ctx) to msn
                for i in channels:
                    nest.Connect(pre=s.parrots_str[i][:len(s.msn_channels[i])], post=s.msn_channels[i],
                                    conn_spec={'rule':'one_to_one'}, syn_spec={'synapse_model': 'static_synapse', 'delay': 10.0, 'weight': 1.0166666666666666, 'receptor_type': 1})
                    nest.Connect(s.parrots_str[i][:len(s.msn_channels[i])], s.msn_channels[i],
                                    conn_spec={'rule':'one_to_one'}, syn_spec={'synapse_model': 'static_synapse', 'delay': 10.0, 'weight': 0.0190625, 'receptor_type': 2})

                    # uncomment for STDP synapses
                    # nest.Connect(pre=s.parrots_str[i][:len(s.msn_channels[i])], post=s.msn_channels[i],
                    #                 conn_spec={'rule':'one_to_one'}, syn_spec={'synapse_model': 'stdp_nestml_rec', 'receptor_type': 1, 'delay':10, 'weight':1})
                    # nest.Connect(s.parrots_str[i][:len(s.msn_channels[i])], s.msn_channels[i],
                    #                 conn_spec={'rule':'one_to_one'}, syn_spec={'synapse_model': 'stdp_nestml_rec', 'receptor_type': 2, 'delay':10, 'weight':0.01})

                # Connect msn to itself, depending on the type of connectivity defined (random or physical)
                s._connect_msn_to_msn()

                # Connect MSN to SNr, depending on the type of connectivity (plastic or static)
                if s.params.net.snr.is_active:
                    syn_spec = s.params.syn.msn_snr.syn_spec_static
                    for i in channels:
                        syn_specc={'synapse_model': 'static_synapse', 'delay': 1.0, 'weight': 57.654700925145136, 'receptor_type': 1}
                        nest.Connect(s.msn_d1_channels[i], s.snr_channels[i], conn_spec={'rule':'pairwise_bernoulli', 'p': 0.033}, syn_spec=syn_specc)
       
        # SNr efferents
        if s.params.net.snr.is_active:
            nest.Connect(s.snr, s.snr, conn_spec=params.syn.snr_snr.conn_spec, syn_spec={'synapse_model': 'static_synapse', 'delay': 1.0, 'weight': 0.32539679447089986, 'receptor_type': 2})
            for i in channels:
                nest.Connect(s.poisson_gpe[i], s.parrots_gpe[i])
                nest.Connect(pre=s.parrots_gpe[i][:len(s.snr_channels[i])], post=s.snr_channels[i],
                                    conn_spec={'rule':'pairwise_bernoulli', 'p': 0.1066}, syn_spec={'synapse_model': 'static_synapse', 'delay': 3.0, 'weight': 59.67154181984671, 'receptor_type': 3})
 
            # adding
            nest.Connect(s.poisson_stn, s.parrots_stn)
            nest.Connect(pre=s.parrots_stn, post=s.snr,
                            conn_spec={'rule':'pairwise_bernoulli', 'p': 0.3}, syn_spec={'synapse_model': 'static_synapse', 'delay': 1.5, 'weight': 3.39200118729286, 'receptor_type': 4})
            nest.Connect(pre=s.parrots_stn, post=s.snr,
                                    conn_spec={'rule':'pairwise_bernoulli', 'p': 0.3}, syn_spec={'synapse_model': 'static_synapse', 'delay': 1.5, 'weight': 0.2 * 3.39200118729286, 'receptor_type': 5})


    def _connect_msn_to_msn(s):
        nest.Connect(s.msn, s.msn, conn_spec={'rule': 'pairwise_bernoulli', 'p': 0.31762652705061084, 'allow_autapses': False}, syn_spec={'synapse_model': 'static_synapse', 'delay': nest.random.uniform(min=1.0, max=2.0), 'weight': 0.25, 'receptor_type': 3})

    def _create_detectors(s):
        s.spikedetector_global = nest.Create('spike_recorder')
        # s.spikedetector_msn_d1 = nest.Create('spike_recorder')
        # s.spikedetector_msn_d2 = nest.Create('spike_recorder')

        if s.params.net.msn.is_active:
        #create a spike detector for each of the channels
            s.spikedetector_msn_channels = [nest.Create("spike_recorder") for _ in range(len(s.params.net.p_channels))]

            #for each channel and corresponding detector, connect the channel to the detector 
            for from_channel, to_detector in zip(s.msn_channels, s.spikedetector_msn_channels):
                nest.Connect(from_channel, to_detector)

            #for each channel in the d1 and corresponding detector, connect the channel to the detector
            s.spikedetector_msn_d1_channels = [nest.Create("spike_recorder") for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.msn_d1_channels, s.spikedetector_msn_d1_channels):
                nest.Connect(from_channel, to_detector)

            #for each channel in the d2 and corresponding detector, connect the channel to the detector
            s.spikedetector_msn_d2_channels = [nest.Create("spike_recorder") for _ in range(len(s.params.net.p_channels))]
            for from_channel, to_detector in zip(s.msn_d2_channels, s.spikedetector_msn_d2_channels):
                nest.Connect(from_channel, to_detector)

            #connect the population of msn to the global spike detector
            nest.Connect(s.msn, s.spikedetector_global)
            s.spikedetector_msn = nest.Create("spike_recorder")

            #connect the population of msn to the local msn spike detector
            nest.Connect(s.msn, s.spikedetector_msn)



            if type(s.params.net.msn.is_active) is bool:
                #if msn is active, create and connect a multimeter on msn
                s.multimeter_msn = nest.Create("multimeter")
                nest.SetStatus(s.multimeter_msn, {"record_from":["V_m", "U_m"]})

                #connect the multimeter to the second neuron in the msn population
                nest.Connect(s.multimeter_msn, s.msn)

            if s.params.net.snr.is_active:
                s.spikedetector_snr_channels = [nest.Create("spike_recorder") for _ in range(len(s.params.net.p_channels))]
                for from_channel, to_detector in zip(s.snr_channels, s.spikedetector_snr_channels):
                    nest.Connect(from_channel, to_detector)
                nest.Connect(s.snr, s.spikedetector_global)
                s.spikedetector_snr = nest.Create("spike_recorder")
                nest.Connect(s.snr, s.spikedetector_snr)

            if s.params.net.gpe.is_active:
                s.spikedetector_gpe_channels = [nest.Create("spike_recorder") for _ in range(len(s.params.net.p_channels))]
                for from_channel, to_detector in zip(s.gpe_channels, s.spikedetector_gpe_channels):
                    nest.Connect(from_channel, to_detector)
                nest.Connect(s.gpe, s.spikedetector_global)
                s.spikedetector_gpe = nest.Create("spike_recorder")
                nest.Connect(s.gpe, s.spikedetector_gpe)


    def plot_activity(s, nucleus_name, channels=True, frequency=True, time_range=None, time_bin=50, n_bins=None, y_range=None):
        nucleus_name = nucleus_name.lower()
        # s.__getattribute__(nucleus_name)
        plt.rcParams['font.size'] = 10

        if type(channels) is bool and channels == False:
            detector = [d for d in s.__getattribute__('spikedetector_' + nucleus_name)]
        elif (type(channels) is bool and channels == True):
            detector = [d[0] for d in s.__getattribute__('spikedetector_' + nucleus_name + '_channels')]
            channels = [i for i, _ in enumerate(detector)]
        elif type(channels) is int:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if channels == idx]
            channels = [channels]
        elif type(channels) is list:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if idx in channels]
        else:
            raise ValueError('Parameter `channel` must be bool, int or list of ints.')

        nucleus_spike_times = [d[0]['times'] for d in [nest.GetStatus(detector[i], 'events') for i, d in enumerate(detector)]]
        n_senders_channel = [len(nest.GetConnections(target=detector[i])) for i,d in enumerate(detector)]

        time_range = (0, nest.GetKernelStatus('biological_time')) if time_range == None else time_range

        if n_bins:
            time_bin = (time_range[1] - time_range[0]) / n_bins
        else:
            n_bins = int((time_range[1] - time_range[0]) / time_bin)

        if y_range: plt.ylim(y_range)

        if frequency:
            for idx, (ch, n_senders) in enumerate(zip(nucleus_spike_times, n_senders_channel)):
                interval, count, __wtf = plt.hist(
                    ch,
                    bins=n_bins,
                    range=time_range,
                    histtype='step',
                    weights=np.ones_like(ch) / (n_senders * (time_bin / 1000.0)),
                    label='Channel {}'.format(channels[idx] + 1) if type(channels) is list else None
                )
            plt.ylabel('Hz')

        else:
            for idx, ch in enumerate(nucleus_spike_times):
                interval, count, __wtf = plt.hist(
                    ch,
                    bins=n_bins,
                    range=time_range,
                    histtype='step',
                    label='ch{}'.format(channels[idx]) if type(channels) is list else None
                )
            plt.ylabel('Spikes inside bin')
        plt.xlabel('Time (ms)\nBin size is {} ms'.format(int(time_bin)))
        if type(channels) is list: plt.legend(loc='upper left')


    def get_activity(s, nucleus_name, channels=True, frequency=True, time_range=None, time_bin=50, n_bins=None):
        nucleus_name = nucleus_name.lower()
        s.__getattribute__(nucleus_name)

        if type(channels) is bool and channels == False:
            detector = s.__getattribute__('spikedetector_' + nucleus_name)
        elif (type(channels) is bool and channels == True):
            detector = [d[0] for d in s.__getattribute__('spikedetector_' + nucleus_name + '_channels')]
        elif type(channels) is int:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if channels == idx]
        elif type(channels) is list:
            detector = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if idx in channels]
        else:
            raise ValueError('Parameter `channel` must be bool, int or list of ints.')

        nucleus_spike_times = [d['times'] for d in nest.GetStatus(detector, 'events')]
        n_senders_channel = [len(nest.GetConnections(target=[d])) for d in detector]

        time_range = (0, nest.GetKernelStatus('time')) if time_range == None else time_range

        if n_bins:
            time_bin = (time_range[1] - time_range[0]) / n_bins
        else:
            n_bins = int((time_range[1] - time_range[0]) / time_bin)

        count_channel = []

        if frequency:
            for ch, n_senders in zip(nucleus_spike_times, n_senders_channel):
                count, interval = np.histogram(
                    ch,
                    bins=n_bins,
                    range=time_range,
                    weights=np.ones_like(ch) / (n_senders * (time_bin / 1000.0))
                )
                count_channel.append(count)
                intervals = interval

        else:
            for ch in nucleus_spike_times:
                count, interval = np.histogram(
                    ch,
                    bins=n_bins,
                    range=time_range,
                )
                count_channel.append(count)
                intervals = interval

        return intervals, count_channel


    def plot_raster(s, nucleus_name, channels=False, time_range=None, time_bin=50, n_bins=None, y_range=None):
        nucleus_name = nucleus_name.lower()
        s.__getattribute__(nucleus_name)

        if type(channels) is bool and channels == False:
            detectors = [s.__getattribute__('spikedetector_' + nucleus_name)]
        elif (type(channels) is bool and channels == True):
            detectors = [d[0] for d in s.__getattribute__('spikedetector_' + nucleus_name + '_channels')]
        elif type(channels) is int:
            detectors = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if channels == idx]
        elif type(channels) is list:
            detectors = [d[0] for idx, d in
                        enumerate(s.__getattribute__('spikedetector_' + nucleus_name + '_channels'))
                        if idx in channels]
        else:
            raise ValueError('Parameter `channel` must be bool, int or list of ints.')

        time_range = (0, nest.GetKernelStatus('time')) if time_range == None else time_range

        if n_bins:
            time_bin = (time_range[1] - time_range[0]) / n_bins
        else:
            n_bins = int((time_range[1] - time_range[0]) / time_bin)

        for d in detectors:
            raster_plot.from_device(d, hist_binwidth=time_bin, ylim=y_range)


    def _compute_signals(s, nucleus='msn', time_range=None, mode='gaussian_convolution'):
        time_range = (0, 2500) if time_range is None else time_range

        if mode == 'gaussian_convolution':
            times, spikes_channels = s.get_activity(nucleus, time_bin=1)
            win = signal.gaussian(M=50, std=7.5)
            sum_win = sum(win)
            s._signals[nucleus] = [signal.convolve(s, win, mode='same')[time_range[0]:time_range[1]]/sum_win for s in spikes_channels]


    def signal(s, nucleus, channels=True, *args, **kwargs):

        try:
            s.__getattribute__('_signals')
        except AttributeError:
            s._signals = {}

        if not nucleus in s._signals:
            s._compute_signals(nucleus, *args, **kwargs)

        if type(channels) is bool:
            if channels==True:
                channels = list(range(len(s.params.net.p_channels)))
            else:
                return None 
        elif type(channels) is int:
            channels = [channels]

        return itemgetter(*channels)(s._signals[nucleus])
