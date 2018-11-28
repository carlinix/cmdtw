# author rcarlini@gmail.com
from libc.stdlib cimport malloc, free

cimport numpy as np
import numpy as np

np.import_array()

from .stepstr import *

# import c fuctions
cdef extern from "cmdtw.c":
    cdef int extra_size(int dp_type)

    cdef struct t_path_element:
        int i
        int j

    cdef struct t_weights:
        double a
        double b
        double c

    cdef struct t_dtw_settings:
        int compute_path
        int dist_type
        double *dist_mapped
        int *dist_mapped_len
        int dp_type
        int window_type
        double window_param
        int norm
        int offset
        t_weights weights

    cdef double cmdtw(double * ref,
                     double * query,
                     int dim,
                     int len_ref,
                     int len_query,
                     double *cost_matrix,
                     t_path_element *path,
                     int *true_path_len,
                     t_dtw_settings dtw_settings)


    cdef int m_dist(double *a, double *b, double *r, int ra, int rb, int dim, int type)

# macros
cdef enum:
    _EUCLID = 11
    _EUCLID_SQUARED = 12
    _MANHATTAN = 13
    _INDEXMATRIX = 14
    _DPW = 20
    _DP1 = 21
    _DP2 = 22
    _DP3 = 23
    _SCP0SYM = 24
    _SCP0ASYM = 25
    _SCP1DIV2SYM = 26
    _SCP1DIV2ASYM = 27
    _SCP1SYM = 28
    _SCP1ASYM = 29
    _SCP2SYM = 210
    _SCP2ASYM = 211

    _NOWINDOW = 00
    _SCBAND = 31
    _PALIVAL = 32
    _ITAKURA = 33
    _PALIVAL_MOD = 35

cdef path_wrapper(t_path_element *cpath, int cpath_len):
    """Path wrapper

    Converts path C array to python list of tuples.

    Args:
        t_path_element *cpath: pointer to t_path_elements array
        int cpath_len: length og the cpath array

    Returns:
        list: list of tuples, each tuple contains row and columns indices

    """
    path = []
    cdef int m
    for m in range(cpath_len):
        path.append((cpath[cpath_len - m - 1].i, cpath[cpath_len - m - 1].j))
    return path


class Setting:
    def __init__(self):
        self._types = {}
        self._cur_type = str()

    def set_type(self, itype):
        """set type method"""
        if itype in self._types.keys():
            self._cur_type = itype
        else:
            print(
            "Unknown type, possible types: \n" + str(self._types.keys()))

    def get_cur_type(self):
        """get type method"""
        return self._cur_type

    def get_cur_type_code(self):
        return self._types[self._cur_type]

    def get_options(self):
        """get options method"""
        return self._types.keys()

    def __str__(self):
        """__str__ method"""
        return str(self._cur_type)


class Dist(Setting):
    """Distance type class
    Contains dintance types for dynamic time wapring algorithm. There are three 
    available distance functions at the moment: 
    'cityblock',
    'euclid',
    'euclid_squared'.
    """

    def __init__(self, dist='euclid'):
        """__init__ method
        Args:
            dist(str, optional): distance type, default is 'cityblock'
        """
        self._cur_type = dist
        self._types = {'euclidean': _EUCLID,
                       'sqeuclidean': _EUCLID_SQUARED,
                       'cityblock': _MANHATTAN,
                       'index_matrix': _INDEXMATRIX
                       }



class Step(Setting):
    """Step class

    Class containts different step patterns for dynamic time warping algorithm.
    There are folowing step patterns available at the moment:
    Usual step patterns:
        'dp1' 
        'dp2' 
        'dp3'
    Sakoe-Chiba classification:
        'p0sym':
        'p0asym':
        'p05sym':
        'p05asym':
        'p1sym':
        'p1asym':
        'p2sym':
        'p2asym':

    You can see step pattern definition using print_step method
    """

    def __init__(self, step='dp2', weights = [1, 1, 1]):
        self._cur_type = step
        self._types = {'dpw': _DPW, 'dp1': _DP1, 'dp2': _DP2, 'dp3': _DP3,
                       'p0sym': _SCP0SYM, 'p0asym': _SCP0ASYM,
                       'p05sym': _SCP1DIV2SYM, 'p05asym': _SCP1DIV2ASYM,
                       'p1sym': _SCP1SYM, 'p1asym': _SCP1ASYM,
                       'p2sym': _SCP2SYM, 'p2asym': _SCP2ASYM}
        self._weights = weights

    def set_weights(self, w):
        self._weights = w

    def get_weights(self):
        return self._weights

    def step_str(self, itype):
        return stepstr[itype]

    def __str__(self):
        return str(self._cur_type)


class Window(Setting):
    """Global constraint class.
    Available constraints: scband, itakura, palival, itakura_mod
    """

    def __init__(self, window='nowindow', param=0.0):
        self._cur_type = window
        self._types = {'scband': _SCBAND, 'palival': _PALIVAL,
                       'itakura': _ITAKURA, 'palival_mod': _PALIVAL_MOD, 'nowindow': _NOWINDOW}

        self._param = param

    def set_param(self, param):
        self._param = float(param)

    def get_param(self):
        return self._param

    def __str__(self):
        return str(self._cur_type + ', parameter: ' + str(self._param))


class Settings:
    """
    class with dtw settings
    """

    def __init__(self,
                 dim=1,
                 dist='cityblock',
                 step='dp2',
                 window='nowindow',
                 param=0.0,
                 norm=False,
                 dist_mapped=[[]],
                 compute_path=False):

        self.dist = Dist(dist)
        self.step = Step(step)
        self.window = Window(window, param)
        self.compute_path = compute_path
        self.norm = norm
        self.dim = dim
        self.dist_mapped = dist_mapped
        if dist_mapped != [[]] and dist == 'index_matrix' :
            self.dist_mapped_len = len(dist_mapped)
        else:
            self.dist_mapped_len = 0
            self.dist_mapped=[[]]

    def __str__(self):
        return str('distance function: ' + str(self.dist) + '\n'
                    'Dimensions: ' + str(self.dim) + '\n'
                    'local constraint: ' + str(self.step) + '\n'
                    'window: ' + str(self.window) + '\n'
                    'normalization: ' + str(self.norm) + '\n')


class cymdtw:
    """
    Main entry, dtw algorithm
    settings is optional parameter
    default settings are dp2 without global constraint
    """

    def __init__(self, ref, query, settings=Settings()):
        self._dist = None
        self._cost = [[]]
        self._dir = [[()]]
        self._path = [()]
        self._dtw(ref, query, settings)

    def _dtw(self, ref, query, settings):
        # sequence control
        if len(ref) == 0 or len(query) == 0 or (query.shape[1] != ref.shape[1]):
            return

        cdef np.ndarray[np.float_t, ndim = 2] d_map
        d_map = np.ascontiguousarray(settings.dist_mapped, dtype=np.float)
        cdef int d_map_len = settings.dist_mapped_len

        # map python settings to C structure dtw_settings functions
        cdef t_dtw_settings c_dtw_settings
        c_dtw_settings.compute_path = <int> settings.compute_path
        c_dtw_settings.dist_type = <int> settings.dist.get_cur_type_code()
        c_dtw_settings.dp_type = <int> settings.step.get_cur_type_code()
        c_dtw_settings.dist_mapped = <double *> d_map.data
        c_dtw_settings.dist_mapped_len = <int *> &d_map_len
        c_dtw_settings.window_type = <int> settings.window.get_cur_type_code()
        c_dtw_settings.window_param = <double> settings.window.get_param()
        c_dtw_settings.norm = <int> settings.norm
        c_dtw_settings.offset = extra_size(c_dtw_settings.dp_type)
        c_dtw_settings.weights.a = <double> settings.step.get_weights()[0]
        c_dtw_settings.weights.b = <double> settings.step.get_weights()[1]
        c_dtw_settings.weights.c = <double> settings.step.get_weights()[2]

        # allocate path
        cdef t_path_element *cpath = <t_path_element *> malloc(sizeof(
                                                               t_path_element) * (<int> (len(ref) + len(query))))

        # init true path length
        cdef int cpath_len = 0

        # init numpy arrays
        cdef np.ndarray[np.float_t, ndim = 2] cref
        cdef np.ndarray[np.float_t, ndim = 2] cquery
        cdef np.ndarray[np.float_t, ndim = 2] cost

        # contiguous c array in memory
        #cref = np.ascontiguousarray(ref, dtype=np.float)
        cref = np.ascontiguousarray(ref, dtype=np.float)
        cquery = np.ascontiguousarray(query, dtype=np.float)

        # init cost matrix
        cost = np.zeros((cref.shape[0]+c_dtw_settings.offset, cquery.shape[0]+c_dtw_settings.offset), dtype=np.float)
        # call cmdtw function (in cmdtw.c)
        self._dist = cmdtw(<double*> cref.data,
                          <double*> cquery.data,
                          <int> settings.dim,
                          <int> len(ref),  #- c_dtw_settings.offset,
                          <int> len(query),  # - c_dtw_settings.offset,
                          <double*> &cost[0,0],
                          cpath,
                          &cpath_len,
                          c_dtw_settings)

        #self._cost = cost[:, :]
        #print("XONGA")
        #pprint(cost)
        self._cost = cost[c_dtw_settings.offset:, c_dtw_settings.offset:]

        # convert c path to python path
        if settings.compute_path:
            self._path = path_wrapper(cpath, cpath_len)

        # cleaning
        free(cpath)

    def get_dist(self):
        return self._dist

    def get_cost(self):
        return self._cost

    def get_path(self):
        return self._path


def dist_matrix(a, b, fdist):

    _dist = Dist(dist=fdist)
    # init numpy arrays
    cdef np.ndarray[np.float_t, ndim = 2] c_a
    cdef np.ndarray[np.float_t, ndim = 2] c_b
    cdef np.ndarray[np.float_t, ndim = 2] c_r

    c_a = np.ascontiguousarray(a, dtype=np.float)
    c_b = np.ascontiguousarray(b, dtype=np.float)
    c_r = np.zeros((c_a.shape[0], c_b.shape[0]), dtype=np.float)

    x = m_dist(<double*> c_a.data,
               <double*> c_b.data,
               <double*> &c_r[0, 0],
               <int> a.shape[0],
               <int> b.shape[0],
               <int> a.shape[1],
               <int> _dist.get_cur_type_code()
              )
    return c_r
