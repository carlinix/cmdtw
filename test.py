import numpy as np
import cmdtw.cymdtw as cymdtw
from fastmdtw.fastmdtw import dtw
import timeit
import pickle
from clip import Clip, ClipList

from cmdtw.cymdtw import dist_matrix



from scipy.spatial.distance import sqeuclidean, euclidean, cityblock
def make_distancesx(a, b):
    return np.vstack([[sqeuclidean(a[i], b[j]) for i in range(a.__len__())] for j in range(b.__len__())])


def select(clipList, clipName):
    for clip in clipList.clips:
        if clip.idx  == clipName:
            return clip
    return None



def  make_distances(a,b) :
    len_a, len_b = np.shape(a), np.shape(b)
    distances = np.ndarray(shape=(len_b[0],len_a[0]), dtype=float)
    for j in range(len_b[0]):
        for i in range(len_a[0]):
            distances[j,i] = sqeuclidean(a[i], b[j])
    return distances



class DistanceMatrix(object):
    md = []

    def __init__(self, dimx=None, dimy=None, points=None):
        # create a distance matrix.
        """

        :rtype: object
        """
        if dimy is None:
            dimy = dimx
        self.dimx = dimx
        self.dimy = dimy
        self.points = points
        DistanceMatrix.md = np.zeros(shape=(self.dimx, self.dimy))
        for x in range(0, self.dimx):
            for y in range(x + 1, self.dimy):
                DistanceMatrix.md[x][y] = DistanceMatrix.md[y][x] = cityblock(points[x], points[y])

    @staticmethod
    def mappeddist(x, y):
        return DistanceMatrix.md[x][y]

    @staticmethod
    def get_all():
        return DistanceMatrix.md



r = np.array([[0,0],[0,0],[0,0],[0,0],[0,0]])
q = np.array([[2,0],[0,0],[0,0],[0,0],[0,0]])
r = np.array([[1,1,0],[0,0,0],[0,0,0]])
q = np.array([[1,1,0],[2,2,0],[3,3,0],[1,2,0]])

print(dist_matrix(r,q,'cityblock'))
exit(0)
#
#
#
# r = np.array([[0],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]])
# q = np.array([[2],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[3],[3]])
# r = np.vstack(np.array([1, 1, 2, 3, 2, 0]))
# q = np.vstack(np.array([0, 1, 1, 2, 3, 2, 1]))
#
# m = make_distancesx(r,q)
# print(m)
#
#
# d = cymdtw.cymdtw(r,q,cymdtw.Settings(dim = r.shape[1], step = 'dp2',#Sakoe-Chiba symmetric step with slope constraint p = 0
#                                             dist='sqeuclidean',
#                                             window = 'nowindow', #type of the window
#                                             param = 2, #window parameter
#                                             norm = False, #normalization
#                                             compute_path = True))
#
#
# #OLD
# dist, a, b = dtw(x=r, y=q, dist=sqeuclidean)
#
# print(b)
# print("---")
# print(a)
# print(d.get_dist())
# print("---")
# print(dist)
# print("---")
# print("---")
# print("---")
# print(d.get_path())
# print("---")
# print(a)
# print("---")
# print("---")
# print("---")
# print(b)
# print("---")
# print(d.get_cost())

EXPERIMENT_DUMP_FILE="/Users/rcarlini/clips/clipsQ-A.pickle"
FEATUREDIM=960
FEATUREDIM=13
# Step 1: Load data

clips = ClipList()
with open(EXPERIMENT_DUMP_FILE, 'rb') as f:
    clips = pickle.load(f)
print(" done", flush=True)
md = None

q = select(clips,'actimel-1')
r = select(clips,'actimel-9')
name_a = r.idx
name_b = q.idx
r = r.descriptors
q = q.descriptors
time1 = timeit.default_timer()
k = max(r.shape[0] * 20 // 100, q.shape[0] * 20 // 100)
print(k)
d = cymdtw.cymdtw(r,q,cymdtw.Settings(dim = r.shape[1], step = 'dp2',#Sakoe-Chiba symmetric step with slope constraint p = 0
                                    dist='sqeuclidean', window = 'nowindow', #type of the window
                                    param = k, #window parameter
                                    norm = False, #normalization
                                    compute_path = True))
time2 = timeit.default_timer()
totaltime = time2 - time1
print("Status C: Q: {} R: {} Dist: {} Time:{}s".format(name_a, name_b, d.get_dist(), totaltime))
time1 = timeit.default_timer()
dist, _, _ = dtw(x=r, y=q, dist=sqeuclidean)
time2 = timeit.default_timer()
totaltime = time2 - time1
print("Status Python: Q: {} R: {} Dist: {} Time:{}s".format(name_a, name_b, dist, totaltime))


#
#
# # Load a csv file with a label and the path of each gist file (one per clip)
# clips = ClipList()
# print("Loading... ", end="", flush=True)
#
# # TODO: if exists clipsY dump__
# # This is the Learning DataSet, these clips aren't used as query
# with open('/Users/rcarlini/clips/clipsQ-A.pickle', 'rb') as f:
#     clips = pickle.load(f)
# print(" done", flush=True)
#
# #TODO: if exists clipsX dump__
# clipsx = ClipList()
# print("Loading... ", end="", flush=True)
# with open('/Users/rcarlini/clips/clipsQ-A.pickle', 'rb') as f:
#   clipsx = pickle.load(f)
# print(" done", flush=True)
#
# for i in range(clipsx.__len__()):
#     name_a = clipsx.clips[i].idx
#     for j in range(i + 1, clipsx.__len__()):
#         name_b = clipsx.clips[j].idx
#         print(name_b,flush=True)
#         r = clipsx.clips[i].descriptors
#         q = clipsx.clips[j].descriptors
#         time1 = timeit.default_timer()
#         k = max(r.shape[0] * 20 // 100, q.shape[0] * 20 // 100)
#         print(k)
#         d = cymdtw.cymdtw(r,q,cymdtw.Settings(dim = r.shape[1], step = 'dp2',#Sakoe-Chiba symmetric step with slope constraint p = 0
#                                             dist='sqeuclidean', window = 'nowindow', #type of the window
#                                             param = k, #window parameter
#                                             norm = False, #normalization
#                                             compute_path = True))
#         time2 = timeit.default_timer()
#         totaltime = time2 - time1
#         print("Status C: Q: {} R: {} Dist: {} Time:{}s".format(name_a, name_b, d.get_dist(), totaltime))
#         time1 = timeit.default_timer()
#         dist, _, _ = dtw(x=r, y=q, dist=sqeuclidean)
#         time2 = timeit.default_timer()
#         totaltime = time2 - time1
#         print("Status Python: Q: {} R: {} Dist: {} Time:{}s".format(name_a, name_b, dist, totaltime))
#
#         #print(d.get_cost())
#         #print(d.get_path())
