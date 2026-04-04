import collections
from .base import *
from .functions import *
from .hermite import *
from .particle import SimpleParticle

def findPair(_dat, _G, _rmax, use_kdtree=False, simple_binary=True):
    """  Find binaries in a particle data set
    The scipy.spatial.cKDTree is used to find pairs

    Parameters
    ----------
    _dat: inhermited SimpleParticle
        Particle data set
    _G: float
        Gravitational constant
    _rmax: float
        Maximum binary separation
    use_kdtree: bool (False)
        If True, use KDtree to find all binaries (slow); otherwise use information from PeTar, only hard binaries are detected (fast)
    simple_binary: bool (True)
        If True, only calculate semi and ecc (fast); otherwise calculating all binary parameters (slow)

    Return
    ----------
    kdt: KDtree structure if use_kdtree=True
    single: type of _dat
        single particle data set
    binary: Binary(simple_mode=simple_binary, member_particle_type=type(single), G=_G)
        binary data set
    """
    if (not issubclass(type(_dat), SimpleParticle)):
        raise ValueError("Data type wrong",type(_dat)," should be subclass of ", SimpleParticle)

    if (use_kdtree):
        # create KDTree
        #print('create KDTree')
        kdt=sp.cKDTree(_dat.pos)
     
        # find all close pairs
        #pairs=kdt.query_pairs(_rmax*AU2PC)
            
        # only check nearest index
        #pair_index=np.unique(np.transpose(np.array([np.array([x[0],x[1]]) for x in pairs])),axis=0)
         
        # find pair index and distance
        #print('Get index')
        r,index=kdt.query(_dat.pos,k=2)
        pair_index=np.transpose(np.unique(np.sort(index,axis=1),axis=0))
        #pair_index = np.transpose(index)

        #index = kdt.query_pairs(_rmax,output_type='ndarray')
        #pair_index = np.transpose(index)
     
        # two members
        p1 = _dat[pair_index[0]]
        p2 = _dat[pair_index[1]]
     
        # check orbits
        #print('Create binary')
        binary = Binary(p1, p2, G=_G, simple_mode=simple_binary)
        apo =binary.semi*(binary.ecc+1.0)
     
        bsel= ((binary.semi>0) & (apo<_rmax))
        binary = binary[bsel]
        
        single_mask = np.ones(_dat.size).astype(bool)
        single_mask[pair_index[0][bsel]]=False
        single_mask[pair_index[1][bsel]]=False
        single = _dat[single_mask]
        return kdt, single, binary
    else:
        idx = _dat.status.argsort()
        dat_sort = _dat[idx]
        status, index, inverse, counts = np.unique(dat_sort.status, return_index=True, return_inverse=True, return_counts=True)
        binary_i1 = index[counts==2]
        binary_i2 = binary_i1+1
        binary = Binary(dat_sort[binary_i1], dat_sort[binary_i2], _G)
        single = dat_sort[index[-1]:]

        return single, binary

def findMultiple(_single, _binary, _G, _rmax, simple_binary=True):
    """  Find triples and quadruples from single and binary data
    The scipy.spatial.cKDTree is used to find pairs

    Parameters
    ----------
    _single: inhermited SimpleParticle
        Single particle data set
    _binary: Binary
        Binary data set
    _G: float
        Gravitational constant
    _rmax: float
        Maximum binary separation
    simple_binary: bool (True)
        If True, only calculate semi and ecc (fast); otherwise calculating all binary parameters (slow)

    Return
    ----------
    kdt: KDtree structure if use_kdtree=True
    single: type of _dat
        single particle data set
    binary: Binary(simple_mode=simple_binary, member_particle_type=type(single), G=_G)
        binary data set
    triple: Binary(p1: type(single), p2: type(binary), G=_G)
        triple data set
    quadruple: Binary(p1: type(binary), p2: type(binary), G=_G)
        quadruple (binary-binary) data set
    """
    if (not issubclass(type(_single), SimpleParticle)):
        raise ValueError("Data type wrong",type(_single)," should be subclass of ", SimpleParticle)

    single_sin = SimpleParticle(_single)
    binary_sin = SimpleParticle(_binary)
    all_sin = join(single_sin, binary_sin)

    # create KDTree
    kdt=sp.cKDTree(all_sin.pos)
     
    # find pair index and distance
    r,index=kdt.query(all_sin.pos,k=2)
    pair_index=np.transpose(np.unique(np.sort(index,axis=1),axis=0))

    bout_i1 = pair_index[0]
    bout_i2 = pair_index[1]

    Ns = _single.size
    Nb = _binary.size
    quad_pre_sel= (bout_i1>=Ns) & (bout_i2>=Ns)
    tri_pre_sel = (bout_i1<Ns) & (bout_i2>=Ns)
    bin_pre_sel = (bout_i1<Ns) & (bout_i2<Ns)

    n_quad_pre = quad_pre_sel.sum()
    n_tri_pre = tri_pre_sel.sum()
    n_bin_pre = bin_pre_sel.sum()
    if (bout_i1.size != n_quad_pre + n_tri_pre + n_bin_pre):
        raise ValueError('Error: multiple index selection size miss match: dat:',bout_i1.size,'quad:',n_quad_pre,'tri:',n_tri_pre,'bin:',n_bin_pre)

    s_del_index=np.array([]).astype(int)
    b_del_index=np.array([]).astype(int)

    quadruple = Binary(member_particle_type = [type(_single), type(_single)], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
    if (quad_pre_sel.sum()):
        q1_index = bout_i1[quad_pre_sel]-Ns
        q2_index = bout_i2[quad_pre_sel]-Ns
        quad_pre = Binary(_binary[q1_index], _binary[q2_index], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
        apo = quad_pre.semi*(quad_pre.ecc+1.0)
        quad_sel = (quad_pre.semi>0) & (apo<_rmax)
        quadruple = quad_pre[quad_sel]
        b_del_index=np.append(q1_index[quad_sel],q2_index[quad_sel])

    triple = Binary(member_particle_type_one = type(_single), 
                    member_particle_type_two = [type(_single), type(_single)], 
                    **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
    if (tri_pre_sel.sum()):
        s_index = bout_i1[tri_pre_sel]
        b_index = bout_i2[tri_pre_sel]-Ns
        tri_pre = Binary(_single[s_index], _binary[b_index], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
        apo = tri_pre.semi*(tri_pre.ecc+1.0)
        tri_sel = (tri_pre.semi>0) & (apo<_rmax)
        triple = tri_pre[tri_sel]
        b_del_index=np.append(b_del_index,b_index[tri_sel])
        s_del_index=s_index[tri_sel]
        
    bmask=np.ones(Nb).astype(bool)
    if (b_del_index.size>0): bmask[b_del_index]=False;
    binary = _binary[bmask]

    if (bin_pre_sel.sum()):
        s1_index = bout_i1[bin_pre_sel]
        s2_index = bout_i2[bin_pre_sel]
        bin_pre = Binary(_single[s1_index], _single[s2_index], **{**_single.initargs, 'G':_G, 'simple_mode':simple_binary})
        apo = bin_pre.semi*(bin_pre.ecc+1.0)
        bin_sel = (bin_pre.semi>0) & (apo<_rmax)
        binary.append(bin_pre[bin_sel])
        s_del_index = np.concatenate((s_del_index, s1_index[bin_sel], s2_index[bin_sel]))

    smask=np.ones(Ns).astype(bool)
    smask[s_del_index]=False
    single = _single[smask]

    return single, binary, triple, quadruple

class BinaryTreeSDAR(DictNpArrayMix):
    """ Binary tree data output from SDAR 
    Keys: (class members)
        semi (1D): semi-major axis
        ecc  (1D): eccentricity
        incline (1D): inclination
        rot_horizon (1D): frame rotational angle in x-y plane (longitude of ascending node)
        rot_self (1D): frame rotational angle in orbital plane (argument of periapsis)
        t_peri (1D): time to peri-center
        period (1D): period
        ecca (1D): eccentric anomaly
        m1   (1D): component 1 mass
        m2   (1D): component 2 mass
        rrel (1D): relative distance
        am   (2D,3): specific angular momemtum x, y, z
        stab (1D): stability factor (>1: unstable)
        sd   (1D): slowdown factor
        sd_org(1D): original slowdown factor based on perturbation
        sd_max(1D): maximum slowdown factor based on timescale criterion
        p1 (member_particle_type) component one
        p2 (member_particle_type) component two
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        Parameters
        ----------
        keyword arguments:
            member_particle_type: type (HermiteParticle)
                Type of component particle 
        """
        member_particle_type=HermiteParticle
        if 'member_particle_type' in kwargs.keys(): member_particle_type=kwargs['member_particle_type']

        keys = [['semi',np.float64], ['ecc',np.float64], ['incline',np.float64],['rot_horizon',np.float64],['rot_self',np.float64],['t_peri',np.float64],['period',np.float64],['ecca',np.float64],['m1',np.float64],['m2',np.float64],['r',np.float64],['am',(np.float64,3)],['stab',np.float64],['sd',np.float64],['sd_org',np.float64],['sd_max',np.float64],['p1',member_particle_type],['p2',member_particle_type]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append,**kwargs)

    def generateBinaryID(self):
        """ Use CantorPairing to map two components id to one binary id
            Add new member bid 
        """
        bid = cantorPairing(self.p1.id, self.p2.id)
        self.addNewMember('bid',bid)
        
class GroupInfo(DictNpArrayMix):
    """ Group information output from PeTar
    Keys: (class members)
        type (1D): group type, 0: new group; 1: end group
        n    (1D): number of members in group (should be consistent with keyword argument N
        time (1D): current time
        pos  (2D,3): position of the group c.m. in the framework of the global system (without shift of global system c.m. if external_mode is on)
        vel  (2D,3): velocity of the group c.m. in the framework of the global system (without shift of global system c.m.)
        bin[X] (BinaryTreeSDAR): members of the group in a hierarchical binary tree
               Here X indicates the order. 0 represents the root (outer most) binary; 1,2,3 ... are inner binaries
               For a triple, bin0 is outer binary, bin1 is inner binary.
               p2 of bin0 is the c.m. of bin1, the id of p2 is the minimum id from the two components in bin1.
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)

        May receive the warning message: 
        RuntimeWarning: invalid value encountered in cast self.__dict__[key] = _dat[:,icol].astype(parameter)
        This is due to the artificial particle mode of mass_bk and status data, which are F64 instead of S64

        Parameters
        ----------
        keyword arguments:
            member_particle_type: type (HermiteParticle)
                Type of component particle, do not change this!
            float_type: type (np.float64)
                floating point data type
            N: int (2)
                Number of members of one group
        """
        keys=[['type',np.int32],['n',np.int32],['time',np.float64],['pos',(np.float64,3)],['vel',(np.float64,3)]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)

        n=2
        if 'N' in kwargs.keys(): n = kwargs['N']
        elif (_dat!=None) & (self.size>0):
            n = int(self.n[0])

        if self.size > 0:
            n_unique = np.unique(self.n.astype(int))
            if n_unique.size > 1:
                raise ValueError(
                    'GroupInfo contains mixed member counts %s. '
                    'Please read one N-member file at a time (e.g. *.group.[rank].nN).'
                    % n_unique.tolist()
                )

        keys_bin = [['bin'+str(i),BinaryTreeSDAR] for i in range(n-1)]
        DictNpArrayMix.__init__(self, keys_bin, _dat, _offset+self.ncols, True, **kwargs)
            
    def generateBinaryID(self, i):
        """ Use CantorPairing to map two components id to one binary id for one binary group
            Add new member bid into this group

        Parameters
        ----------
           i: int 
              binary group index, counting from 0
        """
        key='bin'+str(i)
        if key in self.__dict__.keys():
            self[key].generateBinaryID()
        else:
            raise ValueError('Error: failed to find ',key)
