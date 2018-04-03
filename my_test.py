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
    print(M)

# 计算矩阵的转置
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
        print('row = {}, col = {}:\n{}'.format(rA, cB, M))
        return M
    except ValueError:
        raise ValueError("Matrix A's column number doesn't equal to Matrix b's row number")

# 构造增广矩阵，假设A，b行数相同
def augmentMatrix(A, b):
    MZ = A.copy()
    # MZ = A
    r, c = shape(A)
    for i in range(r):
        MZ[i].append(b[i][0])
    print('row = {}, col = {}:\n{}'.format(r, c+1, np.array(MZ)))
    return MZ

def augmentMatrix1(A, b):
    return [AA + bb for AA, bb in zip(A,b)]

# r1 <---> r2
# 直接修改参数矩阵，无返回值
def swapRows(M, r1, r2):
    temp = M[r1].copy()
    M[r1] = M[r2]
    M[r2] = temp

# r <--- r * scale
# scale为0是非法输入，要求 raise ValueError
# 直接修改参数矩阵，无返回值
def scaleRow(M, r, scale):
    try:
        if scale < 1e-10:
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

r,c = np.random.randint(low=3, high=9, size=2)
# r, c = (3,3)
matrix = np.random.randint(low = -10, high = 10, size = (r, c))
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
print('row = {}, col = {}:\n{}'.format(r, c, matrix))
print('-'*100)
# swapRows(matrix, 0, 1)
# scaleRow(matrix, 1, 2)
addScaledRow(matrix, 0, 1, 2)
print('row = {}, col = {}:\n{}'.format(r, c, matrix))
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
