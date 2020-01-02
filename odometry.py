# %%
import glob
import numpy as np
from scipy.spatial.transform import Rotation as R

def cmp(a):
    return int(a.split('_')[-1].split('.')[0])

filelist = sorted(glob.glob("./assets/cglab/dasan613-odometry_test/odometry/*.txt"), key=cmp)
print('number of files: ', len(filelist))


# #%%
np.set_printoptions(suppress=True)

a = np.zeros((4, 4))
a106 = np.zeros((4, 4))
a0 = np.zeros((4, 4))
a1 = np.zeros((4, 4))
b = np.zeros((4, 4))
b106 = np.zeros((4, 4))

prepost = np.zeros((4, 4))
prepost.flat = [
    0, 0, -1,  0, 
    -1,  0,  0,  0, 
     0,  1, 0,  0, 
     0,  0,  0,  1
]

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
a106.flat = [0.969349,    0.110927,    0.219219,    0.222696,
             -0.112766,    0.993613, -0.00414721,   0.0551511,
             -0.218279,  -0.0207003,    0.975667,    0.205632,
             0,           0,           0,           1]
# Global transformation
b.flat = [0.402660, -0.912283, 0.074861, -0.763035,
          0.896355, 0.409558, 0.169736, -2.470673,
          -0.185507, -0.001244, 0.982642, 1.245954,
          0.000000, 0.000000, 0.000000, 1.000000]
b106.flat = [-0.344988, 0.929564, -0.129980, 5.275046,
             -0.930934, -0.356542, -0.078993, 2.774127,
             -0.119772, 0.093751, 0.988365, 1.510314,
             0.000000, 0.000000, 0.000000, 1.000000]
xzswap = np.asarray(
    [[0, 0, 1, 0],
     [0, 1, 0, 0],
     [1, 0, 0, 0],
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
yflip = np.asarray(
    [[1, 0, 0, 0],
     [0, -1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1]])

def makeRotate(axis, degree):
    rot = np.eye(4)
    rot[0:3, 0:3] = R.from_euler(axis, degree, degrees=True).as_matrix()
    return rot

xrot90 = np.eye(4)
xrot90[0:3, 0:3] = R.from_euler('x', 90, degrees=True).as_matrix()
xrot30 = np.eye(4)
xrot30[0:3, 0:3] = R.from_euler('x', 30, degrees=True).as_matrix()
xrotn90 = np.eye(4)
xrotn90[0:3, 0:3] = R.from_euler('x', -90, degrees=True).as_matrix()
yrot90 = np.eye(4)
yrot90[0:3, 0:3] = R.from_euler('y', 90, degrees=True).as_matrix()
yrotn90 = np.eye(4)
yrotn90[0:3, 0:3] = R.from_euler('y', -90, degrees=True).as_matrix()
zrot90 = np.eye(4)
zrot90[0:3, 0:3] = R.from_euler('z', 90, degrees=True).as_matrix()
zrotn90 = np.eye(4)
zrotn90[0:3, 0:3] = R.from_euler('z', -90, degrees=True).as_matrix()

# # %%
# c = b@a
# print('')
# print('\tb @ a\n', c)
# d = zflip @ c
# print('\tzflip @ b @ a\n', d)
# print('\tTranslation\n', d@[0, 0, 0, 1])
# print('\tRotatoin in matrix\n', d[0:3, 0:3])
# print('\tRotatoin in quaternion (x,y,z,w)\n',
#       R.from_matrix(d[0:3, 0:3]).as_quat())

# # %%
# print(xzswap@xrot90)
# print(xrot90@xzswap)

# c = b@a@zrot90@yrot90
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

# c = zrot90@yrot90@b@a
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

# #**
# c = b@a@xzswap
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

# c = yrot90@zrot90@b@a
# print('Rotation\n\t', repr(R.from_matrix(c[0:3, 0:3]).as_quat()))
# print('Translation\n\t', repr(c@[0, 0, 0, 1]), '\n')

prerot = np.zeros((4,4))
prerot.flat = [
    0, 1, 0, 0,
    0, 0, 1, 0,
    1, 0, 0, 0,
    0, 0, 0, 1
]
prerot = zrot90[:3,:3]@yrotn90[:3,:3]
# prerot[2][0] = -1
print('prerot:\n', prerot)
# print('yrot:\n', R.from_euler('y', 30, degrees=True).as_matrix())

def mat2TQ(t):
    c = prerot@t
    # tmp = R.from_matrix(c[0:3, 0:3]).as_quat()  # xyzw
    tmp = R.from_matrix(np.eye(3)).as_quat()
    Q = np.copy(tmp)
    Q[0] = tmp[3]
    Q[1] = tmp[0]
    Q[2] = tmp[1]
    Q[3] = tmp[2]
    T = c@[0, 0, 0, 1]
    T = T[0:3]
    T[2] = T[2]
    return [T, Q]

# # %%
out = open('output.txt', 'w')
for idx, filepath in enumerate(filelist):
    # print(filepath.split('_')[-1].split('.')[0])
    a = np.loadtxt(filepath)
    # [T, Q] = mat2TQ(a)
    # m = np.concatenate((T, Q))
    # m = m.reshape(1, m.shape[0])
    # np.savetxt(out, m, fmt='%f')
    # print('before:', a)
    # a = prerot @ a
    # print('after:', a)

    position = a@[0,0,0,1]
    position = position.reshape(1, position.shape[0])

    a = np.transpose(a)
    a[:3,:3] = prerot@a[:3,:3]
    # a[:3,:3] = np.transpose(a[:3,:3])

    if 1-np.linalg.norm(a[:3,0]) > 1e-6 or 1-np.linalg.norm(a[:3,1]) > 1e-6 or 1-np.linalg.norm(a[:3,2]) > 1e-6:
        print(np.linalg.norm(a[:3,0]))
        print(np.linalg.norm(a[:3,1]))
        print(np.linalg.norm(a[:3,2]))
        print("wrong matrix")
    print(a[:3,:3])
    # if idx is 123:
    #     print(a[:3,:3])
    #     print(np.linalg.norm(a[:3,0]))
    #     print(np.linalg.norm(a[:3,1]))
    #     print(np.linalg.norm(a[:3,2]))
    # print(a[:3,:3] @ [0,-1,0])

    # print(np.degrees(R.from_matrix(a[:3,:3]).as_euler('xyz')))
    # euler = R.from_matrix(a[:3,:3]).as_euler('xyz')
    # euler[1] = -euler[1]
    # a[:3,:3] = R.from_euler('xyz', euler).as_matrix()[:3,:3]

    # tmp = a[:3,:3]@[0,-1,0]
    # tmp2 = np.copy(tmp)
    # tmp2[2] = -tmp2[2]
    # # print(tmp, tmp2)
    # val = np.dot(tmp, tmp2)
    # if 1 - abs(val) > 1e-6:
    #     theta = np.arccos(np.dot(tmp, tmp2))
    #     rot = R.from_euler('x', theta).as_matrix()
    #     a[:3,:3] = rot @ a[:3,:3]

    np.savetxt(out, a, fmt='%f')
    # np.savetxt(out, m, fmt='%f')

print('finished')
out.close()

# %%
