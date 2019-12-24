# %%
import glob
import numpy as np
from scipy.spatial.transform import Rotation as R

def cmp(a):
    return int(a.split('_')[-1].split('.')[0])

filelist = sorted(glob.glob("./assets/cglab/odometry/*.txt"), key=cmp)
print('number of files: ', len(filelist))


#%%
np.set_printoptions(suppress=True)

a = np.zeros((4, 4))
a0 = np.zeros((4, 4))
a1 = np.zeros((4, 4))
b = np.zeros((4, 4))

# Transformation of 740 frame
a0.flat = [0.491583, 0.861202, -0.129142, 2.16786,
           -0.870831, 0.486281, -0.0720079, -0.251158,
           0.000785807, 0.147858, 0.989008, 0.264198,
           0, 0, 0, 1]
a1.flat = [0.494697,    0.859999,   -0.125202,     2.15838,
           -0.869056,    0.488879,  -0.0757544,   -0.263585,
           -0.00394018,    0.146283,    0.989235,    0.275472,
           0,           0,           0,           1]
a.flat = [0.493464,   0.860172,   -0.12883,    2.16473,
          -0.869765,   0.488211, -0.0718182,  -0.257292,
          0.00112004,   0.147491,   0.989063,   0.263581,
          0,          0,          0,          1]
# Global transformation
b.flat = [0.402660, -0.912283, 0.074861, -0.763035,
          0.896355, 0.409558, 0.169736, -2.470673,
          -0.185507, -0.001244, 0.982642, 1.245954,
          0.000000, 0.000000, 0.000000, 1.000000]

xzswap = np.asarray(
    [[0, 0, 1, 0],
     [0, 1, 0, 0],
     [1, 0, 0, 0],
     [0, 0, 0, 1]]
)
zrot90 = np.asarray(
    [[0, 1, 0, 0],
     [-1, 0, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1]]
)
yrot90 = np.asarray(
    [[0, 0, 1, 0],
     [0, 1, 0, 0],
     [-1, 0, 0, 0],
     [0, 0, 0, 1]]
)
zflip = np.asarray(
    [[1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, -1, 0],
     [0, 0, 0, 1]])
xflip = np.asarray(
    [[-1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1]])

# %%
c = b@a
print('')
print('\tb @ a\n', c)
d = zflip @ c
print('\tzflip @ b @ a\n', d)
print('\tTranslation\n', d@[0, 0, 0, 1])
print('\tRotatoin in matrix\n', d[0:3, 0:3])
print('\tRotatoin in quaternion (x,y,z,w)\n',
      R.from_matrix(d[0:3, 0:3]).as_quat())

# %%
# print(xzswap@xrot90)
# print(xrot90@xzswap)

# c = b@a@zrot90@yrot90
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

# c = zrot90@yrot90@b@a
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

#**
c = b@a@yrot90@zrot90
print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

# c = yrot90@zrot90@b@a
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

def mat2TQ(a):
    c = b@a@yrot90@zrot90
    tmp = R.from_matrix(c[0:3, 0:3]).as_quat()
    Q = np.copy(tmp)
    Q[0] = tmp[3]
    Q[1] = tmp[0]
    Q[2] = tmp[1]
    Q[3] = tmp[2]
    T = c@[0, 0, 0, 1]
    T = T[0:3]
    T[2] = -T[2]
    return [T, Q]

#%%
out = open('output.txt','w')
for filepath in filelist:
    # print(filepath.split('_')[-1].split('.')[0])
    a = np.loadtxt(filepath)
    [T,Q] = mat2TQ(a)
    m = np.concatenate((T,Q))
    m = m.reshape(1, m.shape[0])
    np.savetxt(out, m, fmt='%f')
    # break
print('finished')
out.close()