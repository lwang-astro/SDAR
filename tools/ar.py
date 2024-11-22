# read AR snapshot
from .base import *
from .particle import *

class SDARParticle(SimpleParticle):
    """ AR particle class
    keys: (class members)
        Members inherited from SimpleParticle: see manual of SimpleParticle
        radius (1D): radius
        id (1D): particle id
    """    

    def __init__ (self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        keyword arguments:
            float_type: type (np.float64)
                floating point data type
        """
        if ('float_type' in kwargs.keys()): float_type = kwargs['float_type']
        else: float_type = np.float64

        keys = [['radius', float_type], ['id', np.int64]]

        SimpleParticle.__init__(self, _dat, _offset, _append, **kwargs)
        DictNpArrayMix.__init__(self, keys, _dat, _offset+self.ncols, True, **kwargs)

class SDARProfile(DictNpArrayMix):
    """ SDAR profile
    Keys: (class members)
        n_step_sum (1D): total number of integration step
        n_step_tsyn_sum (1D): total number of steps for time synchronization 
        n_step (1D): number of integration step every output interval
        n_step_tsyn (1D): number of steps for time synchronization every output interval
        if (time_measure):
            prof_total (1D): total time for integration
            prot_int (1D): total time for integration step
            prof_int_tsyn (1D): total time for time synchronization
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys=[["n_step_sum", np.int64], ["n_step_tsyn_sum", np.int64], ["n_step", np.int64], ["n_step_tsyn", np.int64]]
        if ('time_measure' in kwargs.keys()):
            if (kwargs['time_measure']):
                keys = keys + [["prof_total", np.float64], ["prof_int", np.float64], ["prof_int_tsyn", np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SDARInfo(DictNpArrayMix):
    """ SDAR integrator information
    Keys: (class members)
        ds (1D): integration step
        time_offset (1D): time offset to obtain the actual time (time_offset + time)
        r_break_crit (1D): distance criterion to break group (used in Hermite)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys=[["ds", np.float64], ["time_offset", np.float64], ["r_break_crit", np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)   

class SDARData(DictNpArrayMix):
    """ SDAR integrator print column data, used in original SDAR sample code
    Keys: (class members)
        time (1D): current evolved time (counting from zero)
        de (1D): physical energy error
        etot_ref (1D): initial total energy
        ekin (1D): kinetic energy
        epot (1D): potential energy
        gt_drift (1D): time tranformation for drift step
        H (1D): extened phase space Hamiltonian
        H_approx (1D): approximated phase space Hamiltonian
        de_interrupt (1D): energy change due to interruption
        dH_interrupt (1D): H change due to interruption
        perturber (perturber_type): perturber data, depending on the keyword argument 'perturber_type'
                                    if not given, this key is not included
        info (SDARInfo): SDAR information shown as follows:
            ds (1D): integration step
            time_offset (1D): time offset to obtain the actual time (time_offset + time)
            r_break_crit (1D): distance criterion to break group (used in Hermite)
        if (keyword argument 'hybrid' == True):
            hybrid_flag (1D): if 1, hybrid method is used, else, normal method
        profile (SDARProfile): SDAR profile
        if (keyword argument 'slowdown' == True):
            de_sd (1D): slowdown energy error
            etot_sd (1D): slowdown energy
            ekin_sd (1D): slowdown kinetic energy
            epot_sd (1D): slowdown potential energy
            de_sd_change (1D): slowdown energy change
            dH_sd_change (1D): slowdown H change
            de_sd_interrupt (1D): slowdown energy change due to interruption
            dH_sd_interrupt (1D): slowdown H change due to interruption
            sd (SlowDownGroup): slowdown data
        particles (ParticleGroup): particle group, depending on the keyword argument 'member_type' and 'cm_type'

    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            member_type: member particle type (SDARParticle)
            cm_type: c.m. particle type (SDARParticle)
            perturber_type: perturber particle type
            N_particle: int (0)
                Number of particles, determined from file if not provided
            slowdown: bool (False)
                if True, add slowdown keys
            N_sd: int (0)
                Number of slowdown pair, used when slowdown='on'
            time_measure: bool (False)
                if True, add time measure keys in profile
        """
        if ('member_type' in kwargs.keys()):
            kwargs['member_type'] = kwargs['member_type']
        else:
            kwargs['member_type'] = SDARParticle

        if ('cm_type' in kwargs.keys()):
            kwargs['cm_type'] = kwargs['cm_type']
        else:
            kwargs['cm_type'] = SDARParticle

        keys=[['time', np.float64], ['de', np.float64], ["etot_ref",np.float64], ["ekin",np.float64],
              ["epot",np.float64], ['gt_drift', np.float64], ['H', np.float64], ['H_approx', np.float64],
              ['de_interrupt', np.float64], ['dH_interrupt', np.float64]]

        if ('perturber_type' in kwargs.keys()):
            keys = keys + [['perturber', kwargs['perturber_type']]]

        keys = keys + [['info', SDARInfo]]

        if ('hybrid' in kwargs.keys()):
            if (kwargs['hybrid']): 
                keys = keys + [['hybrid_flag', np.int64]]
        keys = keys + [['profile', SDARProfile]]

        if ('slowdown' in kwargs.keys()):
            if (kwargs['slowdown']):
                keys = keys + [['de_sd', np.float64], ["etot_sd",np.float64],["ekin_sd",np.float64],["epot_sd",np.float64], ['de_sd_change', np.float64], ['dH_sd_change', np.float64], ['de_sd_interrupt', np.float64], ['dH_sd_interrupt', np.float64], ['sd', (SlowDownGroup, {'with_indices':True})]]

        keys = keys + [['particles', ParticleGroup]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SDARBinary(DictNpArrayMix):
    """ SDAR binary parameter
    keys: (class members)
        semi (1D): semi-major axis
        ecc (1D): eccentricity
        incline (1D): inclination angle
        rot_horizon (1D): frame rotational angle in x-y plane (longitude of ascending node)
        rot_self (1D): frame rotational angle in orbital plane (argument of periapsis)
        t_peri (1D) : time to peri-center
        period (1D): period
        ecca (1D): eccentric anomaly (-pi, pi)
        m1   (1D): component 1 mass
        m2   (1D): component 2 mass
        rrel (1D): relative distance
        am   (2D,3): specific angular momemtum x, y, z
        stab (1D): stability factor 
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        
        keys = [['semi',np.float64],['ecc', np.float64],['incline',np.float64],['rot_horizon',np.float64],['rot_self',np.float64],['t_peri',np.float64],['period',np.float64],['ecca',np.float64],['m1',np.float64],['m2',np.float64],['rrel',np.float64],['am',(np.float64,3)],['stab',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SDARInterruptBinary(DictNpArrayMix):
    """ SDAR interrupt binary parameter
    keys: (class members)
        time_now (1D): current time
        time_end (1D): targeted integration ending time
        status (1D): binary interruption types
        cm (particle_type): binary c.m. parameter
        bin (SDARBinary): binary orbital parameter
        p1 (particle_type): binary component 1
        p2 (particle_type): binary component 2
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        ----------
        keyword arguments:
            particle_type: type (SDARParticle)
                particle data type
        """

        if (not 'particle_type' in kwargs.keys()):
            kwargs['particle_type'] = SDARParticle
        particle_type = kwargs['particle_type']
        
        keys = [['time_now',np.float64],['time_end', np.float64],['status',np.int64],['cm',particle_type],['bin', SDARBinary],['p1',particle_type],['p2',particle_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
class SlowDown(DictNpArrayMix):
    """ SDAR slowdown data
    Keys: (class members)
        if keyword argument 'with_indices' == True:
            i1 (1D): index of binary component 1
            i2 (1D): index of binary component 2
        sd (1D): slowdown factor 
        sd_org (1D): slowdown original factor without limit
        sd_max (1D): maximum slowdown factor
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            with_indices: bool (False)
                if true, first two keys (columns) are index1 and index2 of slowdown binary components
        """
        keys = []
        if ('with_indices' in kwargs.keys()):
            if kwargs['with_indices']:
                keys = [['i1',np.int64],['i2',np.int64]]
        keys = keys + [["sd", np.float64], ["sd_org", np.float64], ["sd_max", np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

class SlowDownGroup(DictNpArrayMix):
    """ Slowdown data group
    Keys: (class members)
        n (1D): number of slowdown pairs
        sd[x] (SlowDown): slowdown data, [x] indicate the indice, counting from 0
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        Parameters
        ----------
        keyword arguments:
            N_sd: int (0)
                Number of SlowDown pairs
            with_indices: bool (False)
                if true, first two keys in sd[x] are index1 and index2 of slowdown binary components
        """
        keys=[['n', np.int64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n_sd = 0
        if 'N_sd' in kwargs.keys(): n_sd = kwargs['N_sd']
        elif (type(_dat)!=type(None)) & (self.size>0): n_sd = self.n[0]

        if (n_sd>0):
            keys_sd = [['sd'+str(i), SlowDown] for i in range(n_sd)]
            DictNpArrayMix.__init__(self, keys_sd, _dat, _offset+self.ncols, True, **kwargs)

