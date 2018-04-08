from decimal import *
import numpy as np
from pprint import pprint

# 返回矩阵的行数和列数
def shape(M):
    try:
        if isinstance(M, list):
            if not M:
                raise ValueError
        elif not M.any():
            raise ValueError
        cols = len(M[0])
        rows = len(M)
        return (rows, cols)
    except ValueError:
        raise ValueError('The coordinates must be nonempty')
    except TypeError:
        raise TypeError('The coordinates must be an iterable')

# 每个元素四舍五入到特定小数数位
# 直接修改参数矩阵，无返回值
def matxRound(M, decPts=4):
    for row in M:
        for i in range(len(row)):
            row[i] = round(row[i], decPts)
    # print(M)

# 计算矩阵的转置
'''
def transpose(M):
    c, r = shape(M)
    MT = []
    for i in range(r):
        MT.append([])
        for j in range(c):
            MT[i].append(M[j][i])
    print('row = {}, col = {}:\n{}'.format(r, c, MT))
    # pprint((MT))
    return MT
'''
# 这个方法更简单
def transpose(M):
    return [list(col) for col in zip(*M)]

# 计算矩阵乘法 AB，如果无法相乘则raise ValueError
def matxMultiply(A, B):
    try:
        rA, cA = shape(A) # 3x2
        rB, cB = shape(B) # 2x3
        if cA != rB:
            raise ValueError
        M = []
        for i in range(rA):
            M.append([])
            for j in range(cB):
                temp = []
                for k in range(rB):
                    temp.append(A[i][k] * B[k][j])
                M[i].append(sum(temp))
        # print('row = {}, col = {}:\n{}'.format(rA, cB, M))
        return M
    except ValueError:
        raise ValueError("Matrix A's column number doesn't equal to Matrix b's row number")

# 构造增广矩阵，假设A，b行数相同
# def augmentMatrix(A, b):
#     MZ = A.copy()
#     # MZ = A
#     r, c = shape(A)
#     for i in range(r):
#         MZ[i].append(b[i][0])
#     print('row = {}, col = {}:\n{}'.format(r, c+1, np.array(MZ)))
#     return MZ

def augmentMatrix(A, b):
    return [AA + bb for AA, bb in zip(A,b)]

# r1 <---> r2
# 直接修改参数矩阵，无返回值
'''
def swapRows(M, r1, r2):
    temp = M[r1].copy()
    M[r1] = M[r2]
    M[r2] = temp
'''
# 这个方法更简单
def swapRows(M, r1, r2):
    M[r1],M[r2] = M[r2],M[r1]

# r <--- r * scale
# scale为0是非法输入，要求 raise ValueError
# 直接修改参数矩阵，无返回值
def scaleRow(M, r, scale):
    # print('scale = {}'.format(scale))
    try:
        if abs(scale) < 1e-10:
            raise ValueError
        row = M[r]
        for i in range(len(row)):
            row[i] = row[i] * scale
    except ValueError:
        raise ValueError('scale can not be zero.')

# r1 <--- r1 + r2*scale
# 直接修改参数矩阵，无返回值
def addScaledRow(M, r1, r2, scale):
    # scaleRow(M, r2, scale)
    # row = M[r1]
    M[r1] = [x+y*scale for x,y in zip(M[r1], M[r2])]

# https://elonen.iki.fi/code/misc-notes/python-gaussj/index.html
def gauss_jordan(m, eps=1.0e-15):
  """Puts given matrix (2D array) into the Reduced Row Echelon Form.
     Returns True if successful, False if 'm' is singular.
     NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
     Written by Jarno Elonen in April 2005, released into Public Domain"""
  (h, w) = (len(m), len(m[0]))
  for y in range(0,h):
    maxrow = y
    for y2 in range(y+1, h):    # Find max pivot
      if abs(m[y2][y]) > abs(m[maxrow][y]):
        maxrow = y2
    (m[y], m[maxrow]) = (m[maxrow], m[y])
    print('m[y][y] = {}'.format(m[y][y]))
    if abs(m[y][y]) <= eps:     # Singular?
      return False
    for y2 in range(y+1, h):    # Eliminate column y
      c = m[y2][y] / m[y][y]
      for x in range(y, w):
        m[y2][x] -= m[y][x] * c
  for y in range(h-1, 0-1, -1): # Back substitute
    c  = m[y][y]
    for y2 in range(0,y):
      for x in range(w-1, y-1, -1):
        m[y2][x] -=  m[y][x] * m[y2][y] / c
    m[y][y] /= c
    for x in range(h, w):       # Normalize row y
      m[y][x] /= c
  return True

# 实现 Gaussain Jordan 方法求解 Ax = b
""" Gaussian Jordan 方法求解 Ax = b.
    参数
        A: 方阵
        b: 列向量
        decPts: 四舍五入位数，默认为4
        epsilon: 判读是否为0的阈值，默认 1.0e-16

    返回列向量 x 使得 Ax = b
    返回None，如果 A，b 高度不同
    返回None，如果 A 为奇异矩阵
"""
def gj_Solve(A, b, decPts=4, epsilon=1.0e-15):
    (r, c) = shape(A)
    (rb, cb) = shape(b)
    if r != rb:
        return None
    Ab = augmentMatrix1(A, b)
    for y in range(0, r):
        maxrow = y
        for y2 in range(y+1, r):    # Find max pivot
            if abs(Ab[y2][y]) > abs(Ab[maxrow][y]):
                maxrow = y2
        swapRows(Ab, maxrow, y)
        if abs(Ab[y][y]) <= epsilon:
            return None

        for y2 in range(y+1, r):    # Eliminate column y
            co = Ab[y2][y] / Ab[y][y]
            for x in range(y, c):
                Ab[y2][x] -= Ab[y][x] * co
    print('Ab = {}'.format(np.array(Ab)))
    for y in range(r-1, 0-1, -1): # Back substitute
        co  = Ab[y][y]
        for y2 in range(0,y):
            for x in range(c-1, y-1, -1):
                Ab[y2][x] -=  Ab[y][x] * Ab[y2][y] / co
        Ab[y][y] /= co
        for x in range(r, c):       # Normalize row y
            Ab[y][x] /= co
    print('Ab = {}'.format(Ab))
    return None

""" Gaussian Jordan 方法求解 Ax = b.
    参数
        A: 方阵
        b: 列向量
        decPts: 四舍五入位数，默认为4
        epsilon: 判读是否为0的阈值，默认 1.0e-16

    返回列向量 x 使得 Ax = b
    返回None，如果 A，b 高度不同
    返回None，如果 A 为奇异矩阵
"""
def gj_Solve1(A, b, decPts=4, epsilon=1.0e-15):
    (r, c) = shape(A)
    (rb, cb) = shape(b)
    if r != rb:
        return None
    Ab = augmentMatrix(A, b)
    for j in range(0, c): #对于Ab的每一列（最后一列除外）
        maxrow = j
        for ii in range(j+1, r): #寻找列 j中，对角线以及对角线以下所有元素的绝对值的最大值
            if abs(Ab[ii][j]) > abs(Ab[maxrow][j]):
                maxrow = ii
        # print('maxrow = {}'.format(maxrow))
        if abs(Ab[maxrow][j]) <= epsilon:#如果绝对值最大值为0，则为奇异矩阵
            return None
        swapRows(Ab, maxrow, j)#使用第一个行变换，将绝对值最大值所在行交换到对角线元素所在行（行j）

        # 使用第二个行变换，将列 j的对角线元素缩放为1
        coefficient = Ab[j][j]
        scaleRow(Ab, j, 1/coefficient)

        # 多次使用第三个行变换，将列 j的其他元素消为 0
        for i in range(0, r):
            if i == j:
                continue
            coefficient = Ab[i][j]
            addScaledRow(Ab, i, j, -coefficient)

    matxRound(Ab, decPts)
    # print('Ab = \n{}'.format(np.array(Ab)))
    # print('Ab = \n{}'.format(Ab))
    Mr = []
    for i in range(0, r):
        Mr.append([])
        Mr[i].append(Ab[i][-1])
    # Mret = transpose(Ab)
    # print('Mret = \n{}'.format(Mret))
    # vec = Mret[-1]
    return Mr

from helper import *
'''
参数：X, Y 存储着一一对应的横坐标与纵坐标的两个一维数组
返回：线性回归的系数(如上面所说的 m, b)
'''
def linearRegression2D(X, Y):
    XX,YY = [],[]
    for i in range(len(X)):
        XX.append([X[i]])
        YY.append([Y[i]])
    # print('*'*200)
    # print('XX = {}'.format(XX))
    # print('YY = {}'.format(YY))
    allOne = [[1.0]]*len(X)
    # allOne.append([1]*len(X))
    # print('\nallOne = {}'.format(allOne))

    x = augmentMatrix(XX, allOne)
    # r,c = shape(x)
    # print('*'*250)
    # print('x = \nr = {}, c = {}\n{}'.format(r,c,x))

    xt = transpose(x)
    # print('#'*250)
    # r,c = shape(xt)
    # print('xt = \nr = {}, c = {}\n{}'.format(r,c,xt))

    # b等于X的转置乘以Y
    b = matxMultiply(xt, YY)
    # r,c = shape(b)
    # print('b = \nr = {}, c = {}\n{}'.format(r,c,b))

    # A等于X的转置乘以X
    A = matxMultiply(xt, x)
    # r,c = shape(A)
    # print('A = \nr = {}, c = {}\n{}'.format(r,c,A))

    h = gj_Solve1(A, b)
    # r,c = shape(h)
    # print('h = \nr = {}, c = {}\n{}'.format(r,c,h))

    return h[0][0], h[1][0]

def linearRegression(X,Y):
    XX = X
    YY = []
    for i in range(len(X)):
        YY.append([Y[i]])

    # print('XX = {}'.format(XX))
    # print('*'*200)
    # print('YY = {}'.format(YY))

    allOne = [[1.0]]*len(X)
    # allOne.append([1]*len(X))
    # print('\nallOne = {}'.format(allOne))
    # print('len(allOne) = {}'.format(len(allOne)))

    x = augmentMatrix(XX, allOne)
    # r,c = shape(x)
    # print('*'*250)
    # print('x = \nr = {}, c = {}\n{}'.format(r,c,x))

    xt = transpose(x)
    # print('#'*250)
    # r,c = shape(xt)
    # print('xt = \nr = {}, c = {}\n{}'.format(r,c,xt))

    # b等于X的转置乘以Y
    b = matxMultiply(xt, YY)
    # r,c = shape(b)
    # print('b = \nr = {}, c = {}\n{}'.format(r,c,b))

    # A等于X的转置乘以X
    A = matxMultiply(xt, x)
    # r,c = shape(A)
    # print('A = \nr = {}, c = {}\n{}'.format(r,c,A))

    h = gj_Solve1(A, b)
    r,c = shape(h)
    print('h = \nr = {}, c = {}\n{}'.format(r,c,h))
    print('^'*160)
    # print((zip(*h)))
    hh = zip(*h)
    print(hh)
    return [h[0][0], h[1][0], h[2][0]]

seed = 1978
X_3d, Y_3d = generatePoints3D(seed)
# vs_scatter_3d(X_3d, Y_3d)
r, c = shape(X_3d)
print('X_3d = \nr = {}, c = {}\n{}'.format(r, c, X_3d))
print('len(X_3d) = {}'.format(len(X_3d)))
print('_'*259)
# r, c = shape(Y_3d)
# print('Y_3d = \nr = {}, c = {}\n{}'.format(r, c, Y_3d))
print('Y_3d = {}'.format(Y_3d))
print('len(Y_3d) = {}'.format(len(Y_3d)))
print('-'*259)
coeff = linearRegression(X_3d, Y_3d)
print('coeff = {}'.format(coeff))
# vs_scatter_3d(X_3d, Y_3d, coeff)

'''
# 下面测试二维线性回归
seed = 1978
X,Y = generatePoints2D(seed)
# vs_scatter_2d(X, Y)
# print('X = {}'.format(X))
# print('Y = {}'.format(Y))
# print('-'*260)

m2,b2 = linearRegression2D(X,Y)
print('m2 = {}, b2 = {}'.format(m2, b2))
'''

'''
r,c = np.random.randint(low=3, high=9, size=2)
# r, c = (3,3)
# matrix = np.random.randint(low = -10, high = 10, size = (r, c))
matrix0 = [[5, -10, 0],
           [1, 2,   1],
           [3, 5, 3]]
matrix1 = [[5, -10, 0, 1],
           [1,   2, 1, 1],
           [3,   5, 3, 1]]
matrix2 = [[1, 1, 3, 1],
           [-10, -4, -8, 1],
           [-5, 1, 7, 1]]
matrix = matrix0
b = [[1], [1], [1]]
r, c = shape(matrix)
print('row = {}, col = {}:\n{}'.format(r, c, matrix))
print('-'*100)
h = gj_Solve1(matrix0, b)
print('*'*100)
print('h = {}'.format(h))
'''
# print(gauss_jordan(matrix))
# matrixA = np.random.randint(low=1, high=3, size=(r,c))
# matrixB = np.random.randint(low=0, high=4, size=(c,r))
# A = np.random.randint(low=-10,high=10,size=(r,c))
# b = np.random.randint(low=-10,high=10,size=(r,1))
# Amat = A.tolist()
# if A.tolist() == Amat:
#     print('haha'*10)
# if A.tolist() == Amat:
#     print('NiNi'*10)
# print(type(matrixB))
# matrix = np.random.random((r,c)) * 2
# matrix = [[1,2,3,4],
#           [5,6,7,8],
#           [9,10,11,12],
#           [13,14,15,16],
#           [17,18,19,20]]
# matrix1 = 99
# print(shape(matrix))

# swapRows(matrix, 0, 1)
# scaleRow(matrix, 1, 2)
# addScaledRow(matrix, 0, 1, 2)
# print('row = {}, col = {}:\n{}'.format(r, c, matrix))
# print('row = {}, col = {}:\n{}'.format(r, 1, b))
# print('*'*100)
# transpose(matrix)
# matxMultiply(matrixA, matrixB)
# ab = np.hstack((A,b))
# print(ab)
# Ab = augmentMatrix1(Amat, b.tolist())
# print('row = {}, col = {}:\n{}'.format(r, c+1, np.array(Ab)))
# if A.tolist() == Amat:
#     print('heihei'*10)
# print('\\'*101)
# r,c=shape(A)
# print('row = {}, col = {}:\n{}'.format(r,c, A.tolist()))
# print('/'*101)
# r,c=shape(Amat)
# print('row = {}, col = {}:\n{}'.format(r,c, Amat))
