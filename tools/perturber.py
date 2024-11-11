from .base import *
from .ar import *
from .functions import *

class HermitePerturber(DictNpArrayMix):
    """ Hermite perturber data
    Keys: (class members)
        r_min_index (1D):    nearest neighbor index for each ptcl  
        mass_min_index (1D): mimimum mass in neighbors         
        r_min_sq (1D):       nearest neighbor distance square      
        r_min_mass (1D):     nearest neighbor index for each ptcl 
        mass_min (1D):       mimimum mass in neighbors             
        r_neighbor_crit_sq (1D): neighbor radius criterion
        need_resolve_flag (1D): indicate whether the members need to be resolved for outside 
        initial_step_flag (1D): indicate whether the time step need to be initialized due to the change of neighbors
        n_neighbor_group (1D): number of group neighbor
        n_neighbor_single (1D): number of single neighbor
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """

        key = [['r_min_index', np.int64], ['mass_min_index', np.int64], ['r_min_sq', np.float64], ['r_min_mass', np.float64], ['mass_min', np.float64], ['r_neighbor_crit_sq', np.float64], ['need_resolve_flag', bool], ['initial_step_flag', bool], ['n_neighbor_group', np.int64], ['n_neighbor_single', np.int64]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    