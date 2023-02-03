import nest
from matplotlib import pyplot as plt
import numpy as np
import itertools

# , log_progress
from utils import raster_plot, verbose, transient_selectivity, steadystate_selectivity
# from tqdm import tqdm
from basalganglia import BasalGanglia
from bg_params import params
# from . import restart_kernel


class Experiment:

    time_bin = 20.0
    interpolation_kind = ['linear', 'cubic'][1]

    def __init__(self):
        raise NotImplementedError(
            "'Experiment' is a virtual class, you need to inherit it.")

    def run(self):
        raise NotImplementedError(
            "'Experiment' is a virtual class, you need to inherit it.")

    def plot(self):
        raise NotImplementedError(
            "'Experiment' is a virtual class, you need to inherit it.")

    def transient_selectivity(s, nucleus='msn'):
        return transient_selectivity(s.bg.signal(nucleus, 0), s.bg.signal(nucleus, 1), t_stable=(2000, 2500), t_transit=(1500, 1700))

    def steadystate_selectivity(s, nucleus='msn'):
        return steadystate_selectivity(s.bg.signal(nucleus, 1), t_pre=(1250, 1500), t_stable=(2000, 2500))


class ExperimentDA(Experiment):

    # Frecuencias de disparo de la neurona cortical
    cortex_low_freq = 2.2
    cortex_high_freq = 3.6

    def __init__(s, dopamine_level=0.3, bg=None):

        # Initial values
        # restart_kernel()
        s.bg = BasalGanglia() if bg is None else bg

        # Setting dopamine
        s.dopamine_level = dopamine_level
        s.bg.params.neuron.msn_d1.y1 = dopamine_level
        s.bg.params.neuron.msn_d1.y2 = dopamine_level
        s.bg.params.neuron.msn_d2.y1 = dopamine_level
        s.bg.params.neuron.msn_d2.y2 = dopamine_level

        # Configure the model
        s.bg.configure()

    # @verbose(_params.sim.verbose)

    def run(s, warming_up_time=1500):
        m = s.cortex_low_freq * s.bg.params.net.ctx_afferents_per_neuron
        M = s.cortex_high_freq * s.bg.params.net.ctx_afferents_per_neuron

        stn_low = 10
        stn_high = 35

        nest.SetStatus(s.bg.poisson[0], {'rate': m})
        nest.SetStatus(s.bg.poisson[1], {'rate': m})
        nest.SetStatus(s.bg.poisson_stn, {'rate': stn_low})
        # nest.SetStatus(s.bg.poisson_gpe[0], {'rate': 0.7})
        # nest.SetStatus(s.bg.poisson_gpe[1], {'rate': 0.3})
        # nest.Simulate(warming_up_time)

        nest.SetStatus(s.bg.poisson_gpe[0], {'rate': 1.2})
        nest.SetStatus(s.bg.poisson_gpe[1], {'rate': 1.35})
        nest.Simulate(500)

        nest.SetStatus(s.bg.poisson_gpe[0], {'rate': 0.8})
        nest.SetStatus(s.bg.poisson_gpe[1], {'rate': 0.8})
        nest.Simulate(1000)

        # Evidence accumulation period
        for i in range(25):
            fi = i/50.0
            nest.SetStatus(s.bg.poisson[0], {'rate': m + (M-m)*fi})
            nest.SetStatus(s.bg.poisson[1], {'rate': m + (M-m)*fi})
            nest.SetStatus(s.bg.poisson_stn, {
                           'rate': stn_low + (stn_high-stn_low)*fi})
            nest.SetStatus(s.bg.poisson_gpe[0], {'rate':  0.8})
            nest.SetStatus(s.bg.poisson_gpe[1], {'rate':  0.6})

            nest.Simulate(1.0)

        # Action selection
        for i in range(25):
            fi = i/50.0
            nest.SetStatus(s.bg.poisson[0], {'rate': (m+M)*0.5 + (M-m)*fi})
            nest.SetStatus(s.bg.poisson[1], {'rate': (m+M)*0.5 - (M-m)*fi})
            nest.SetStatus(s.bg.poisson_stn, {'rate': (
                stn_high+stn_low)*0.5 - (stn_high-stn_low)*fi})
            nest.SetStatus(s.bg.poisson_gpe[0], {'rate':  0.9})
            nest.SetStatus(s.bg.poisson_gpe[1], {'rate':  0.6})
            nest.Simulate(1.0)

        # Final 1000ms
        nest.SetStatus(s.bg.poisson_gpe[0], {'rate':  0.1})
        nest.Simulate(1000.0)


    def plot(s):
        s._plot_network()
        plt.show()

    # @verbose(_params.sim.verbose)
    def _plot_network(s):
        plt.rcParams['figure.figsize'] = [3.0, 8.0]
        for i in range(len(_params.net.p_channels)):
            raster_plot.from_device(
                s.bg.spikedetector_msn_channels[i], hist_binwidth=s.time_bin, title='Canal {}'.format(i + 1), ylim=(0, 40))
