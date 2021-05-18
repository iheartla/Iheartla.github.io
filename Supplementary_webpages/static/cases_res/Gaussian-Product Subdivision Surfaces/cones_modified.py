#------------------------------------------------------------------------------------------
# Gaussian-Product Subdivision Surfaces - Python Demo
# Copyright (c) 2020 Reinhold Preiner
#
# This code is licensed under the MIT license. 
# https://opensource.org/licenses/mit-license.html
#------------------------------------------------------------------------------------------


import numpy as np
import igl
import os
root_folder = ".."

def loop_gps(v, f, c, n):
    """
    :param :v : ℝ^(m×3)
    :param :f : ℤ^(t×3)
    :param :c : ℝ^(3×3)
    :param :n : ℤ
    """
    v = np.asarray(v, dtype=np.float64)
    f = np.asarray(f, dtype=np.integer)
    c = np.asarray(c, dtype=np.float64)

    m = v.shape[0]
    t = f.shape[0]
    dim_0 = c.shape[0]
    assert v.shape == (m, 3)
    assert f.shape == (t, 3)
    assert c.shape == (dim_0, 3, 3)
    assert np.ndim(n) == 0
    assert dim_0 == m 

    icf = np.zeros((dim_0, 9, ))
    for i in range(1, dim_0+1):
        icf[i-1] = np.matrix.flatten(np.linalg.inv(c[i-1]), order='F')
    qq1 = np.zeros((dim_0, 3, ))
    for i in range(1, dim_0+1):
        qq1[i-1] = np.hstack((icf[i-1][1-1], icf[i-1][2-1], icf[i-1][3-1]))
    qq2 = np.zeros((dim_0, 3, ))
    for i in range(1, dim_0+1):
        qq2[i-1] = np.hstack((icf[i-1][5-1], icf[i-1][6-1], icf[i-1][9-1]))
    qlin = np.zeros((m, 3, ))
    for i in range(1, m+1):
        qlin[i-1] = np.linalg.solve(c[i-1], v[i-1, :])

    # perform Gaussian-product subdivision
    # note: igl.loop only handles 3D subdivs, so we split the 9D meshes into three 3D ones 
    for _ in range(n):
        qq1, dmy = igl.loop(qq1, f)
        qq2, dmy = igl.loop(qq2, f)
        qlin, f  = igl.loop(qlin, f) 

    # transform back to 3D
    m = qlin.shape[0]
    v_out = np.zeros((m, 3, ))
    for i in range(1, m+1):
        v_out_i_2 = np.vstack((qq1[i-1].T.reshape(1, 3), np.hstack((qq1[i-1][2-1], qq2[i-1][1-1], qq2[i-1][2-1])).T.reshape(1, 3), np.hstack((qq1[i-1][3-1], qq2[i-1][2-1], qq2[i-1][3-1])).T.reshape(1, 3)))
        v_out[i-1] = np.linalg.solve(v_out_i_2, qlin[i-1])

    return v_out, f
#------------------------------------------------------------------------------------------

# create a cone
v = np.array([
    [0, 1.5, 0],
    [1.41421, -1.5, -1.41421],[0, -1.5, -2], [-1.41421, -1.5, -1.41421],[-2, -1.5, 0],
    [-1.41421, -1.5, 1.41421], [0, -1.5, 2], [1.41421, -1.5, 1.41421],[2, -1.5, 0]
])
f = np.array([[0, 1, 2],[0, 2, 3], [0, 3, 4], [0, 4, 5],[0, 5, 6],[0, 6, 7],[0, 7, 8],[0, 8, 1]])


# initialize with isotropic covariances
c = [np.identity(3)*0.01 for i in range(len(v))]
# print("v:{}, f:{}, c:{}".format(v.shape, f.shape, c[0].shape))
vs, fs = loop_gps(v, f, c, 5)
# print("vs:{}, shape:{}".format(vs, vs.shape))
 
# make apex covariance sharper
c[0] *= 0.001
vs2, fs2 = loop_gps(v, f, c, 5)
# print("vs2:{}".format(vs2))

# make apex covariance vertically anisotropic (reduce variance in x direction)
c[0] = [[0.001,0,0],
        [0,0.01,0],
        [0,0,0.01]]
vs3, fs3 = loop_gps(v, f, c, 5)
# print("vs3:{}".format(vs3))

# make apex covariance horizontally anisotropic (reduce variance in y direction)
c[0] = [[0.01,0,0],
        [0,0.001,0],
        [0,0,0.01]]
vs4, fs4 = loop_gps(v, f, c, 5)
# print("vs4:{}".format(vs4))

# make apex covariance horizontally anisotropic and super flat
c[0] = [[1,0,0],
        [0,0.00001,0],
        [0,0,1]]
vs5, fs5 = loop_gps(v, f, c, 5)
# print("vs5:{}".format(vs5))



# save meshes to .off files
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone.off"), v, f)
print("Wrote control mesh to", os.path.join(root_folder, "data", "cone.off"))

igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-pointy-gps.off"), vs2, fs2)
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-concave-gps.off"), vs3, fs3)
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-plateau-gps.off"), vs4, fs4)
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-cylinder-gps.off"), vs5, fs5)
print("Wrote GPS meshes to", os.path.join(root_folder, "data", "cone-<variant>-gps.off"))
