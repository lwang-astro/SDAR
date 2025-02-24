# read AR snapshot
from .base import *
from .particle import *
from .ar import *

class HermiteParticle(SDARParticle):
    """ Hermite particle type
    keys: (class members)
        Members inherited from SDARParticle: see manual of SDARParticle
        dt    (1D): time step
        time  (1D): current time
        acc   (2D,3): acceleration x, y, z
        jerk  (2D,3): acceleration derivative x, y, z
        pot   (1D): potential
    """

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        keyword arguments:
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['dt', float_type], ['time', float_type], ['acc', (float_type, 3)], ['jerk', (float_type, 3)], ['pot', float_type]]

        SDARParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class HermiteEnergy(DictNpArrayMix):
    """ Hermite integrator energy data
    keys: (class members)
        de (1D): energy error 
        etot_ref (1D): initial total energy
        ekin (1D): kinetic energy
        epot (1D): potential energy
        epert (1D): perturbation energy
        de_change (1D): cumulative energy change
        de_interrupt (1D): interrupt energy change (binary stellar evolution)
        de_modify (1D): modify energy change (stellar evolution)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [["de",np.float64],["etot_ref",np.float64],["ekin",np.float64],["epot",np.float64],["epert",np.float64],["de_change",np.float64],["de_interrupt",np.float64],["de_modify",np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)


class HermiteProfile(DictNpArrayMix):
    """ Hermite profile
    Keys: (class members)
        h4_step_single (1D): single particle total steps
        h4_step_group (1D): AR group total steps
        ar_step (1D): ar total steps
        ar_step_tsyn (1D): ar time synchronize steps
        break_group (1D): number of break groups
        new_group (1D): number of new groups
        if (time_measure):
            prof_total (1D): total time for integration
            prot_h4_single (1D): total time for hermite single particle integration
            prof_h4_group (1D): total time for hermite group integration
            prof_adjust (1D): total time for adjust groups
            prof_init (1D): total time for initialization
            prof_modify_single (1D): total time for modify single
            prof_select_act (1D): total time for select active particles
            prof_ar (1D): total time for AR integration
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys=[["h4_step_single", np.int64], ["h4_step_group", np.int64], ["ar_step", np.int64], ["ar_step_tsyn", np.int64], ["break_group", np.int64], ["new_group", np.int64]]
        if ('time_measure' in kwargs.keys()):
            if (kwargs['time_measure']):
                keys = keys + [["prof_total", np.float64], ["prof_h4_single", np.float64], ["prof_h4_group", np.float64], ["prof_adjust", np.float64], ["prof_init", np.float64], ["prof_modify_single", np.float64], ["prof_select_act", np.float64], ["prof_ar", np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class HermiteData(DictNpArrayMix):
    """ Hermite+SDAR integrator print column data, used in petar.hard.debug
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        time_offset (1D): time offset to calculate the global time (time+time_offset)
        energy_phy (HermiteEnergy): physical energy data
        energy_sd (HermiteEnergy): slowdown energy data
        sd (SlowDownGroup): slowdown data
        profile (HermiteProfile): hermite profile
        particles (ParticleGroup): particle group, depending on the keyword argument 'member_type' and 'cm_type'
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            member_type: member particle type (HermiteParticle)
            cm_type: c.m. particle type (SDARParticle)
            N_particle: int (0)
                Number of members of one group
            N_sd: int (0)
                Number of slowdown pairs
            time_measure: bool (False)
                if True, add time measure keys in profile
        """

        if ('member_type' in kwargs.keys()):
            kwargs['member_type'] = kwargs['member_type']
        else:
            kwargs['member_type'] = HermiteParticle

        if ('cm_type' in kwargs.keys()):
            kwargs['cm_type'] = kwargs['cm_type']
        else:
            kwargs['cm_type'] = SDARParticle

        keys=[['time', np.float64], ['time_offset', np.float64], ['energy_phy', HermiteEnergy], ['energy_sd', HermiteEnergy], ['sd', SlowDownGroup]]
        keys = keys + [['profile', HermiteProfile], ['particles', ParticleGroup]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
