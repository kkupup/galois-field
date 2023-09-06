import numpy as np
import itertools


class GF:
    def __init__(self, w):
        self.w = w
        #  有限域元素个数
        self.num = 1 << w
        #  初始化反对数表与对数表
        self.gfilog = [-1] * self.num
        self.gflog = [-1] * self.num
        #  生成反对数表与对数表
        self.Create_GFilog()
        self.Create_GFlog()

    #  本原多项式获取
    def Primitive_Polynomial(self):
        switch = {
            3: 0b1011,
            4: 0b10011,
            8: 0b100011101
        }
        return switch[self.w]

    #  创建反对数表
    def Create_GFilog(self):
        #  本获取原多项式
        primitive_poly = self.Primitive_Polynomial()
        num = len(self.gfilog)
        #  x^0 = 1
        self.gfilog[0] = 1
        #  生成表
        for i in range(1, num - 1):
            tmp = self.gfilog[i - 1] << 1
            #  大于num时模本原多项式
            while tmp >= num:
                tmp ^= primitive_poly << (tmp.bit_length() - primitive_poly.bit_length())
            self.gfilog[i] = tmp

    #  根据反对数表生成对数表
    def Create_GFlog(self):
        for i in range(0, len(self.gfilog) - 1):
            self.gflog[self.gfilog[i]] = i

    #  有限域内加法
    @staticmethod
    def addition(para1, para2):
        return para1 ^ para2

    #  有限域内减法
    @staticmethod
    def subtraction(para1, para2):
        return para1 ^ para2

    #  有限域内乘法
    def multiplication(self, para1, para2):
        #  乘数为零直接返回0
        if para1 == 0 or para2 == 0:
            return 0
        #  NM = gfilog( ( gflogN + gflogM ) % 2^w - 1 )
        #  运用对数简化乘法运算
        para1, para2 = self.gflog[para1], self.gflog[para2]
        sum = (para1 + para2) % (self.num - 1)
        return self.gfilog[sum]

    #  有限域内除法
    def division(self, para1, para2):
        #  除数与被除数的判别
        if para1 == 0 or para2 == 0:
            return 0
        #  与乘法同理
        para1, para2 = self.gflog[para1], self.gflog[para2]
        sum = para1 - para2
        if sum < 0:
            sum += (self.num - 1)
        sum = sum % (self.num - 1)
        return self.gfilog[sum]

    #  有限域内矩阵加法
    def matrix_add(self, mat1, mat2):
        shape = mat1.shape
        #  是否同型
        if shape != mat2.shape:
            return None
        #  初始化全零矩阵
        mat = np.zeros((shape[0], shape[1]), dtype=int)
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                mat[i, j] = self.addition(mat1[i, j], mat2[i, j])
        return mat

    #  有限域内矩阵乘法
    def matrix_multiplication(self, mat1, mat2):
        shape1, shape2 = mat1.shape, mat2.shape
        #  左矩阵列数等于右矩阵行数
        if shape1[1] != shape2[0]:
            return None
        #  初始化全零矩阵
        mat = np.zeros((shape1[0], shape2[1]), dtype=int)
        for i in range(shape1[0]):
            for j in range(shape2[1]):
                #  获取mat1的i行， mat2的j列运算
                tmp1 = mat1[i]
                tmp2 = mat2[:, j]
                num = 0
                for k in range(len(tmp1)):
                    num = self.addition(self.multiplication(tmp1[k], tmp2[k]), num)
                mat[i, j] = num
        return mat

    def determinant_value(self, mat):
        shape = np.shape(mat)
        tmp_list = []
        for i in range(shape[0]):
            tmp_list.append(i)
        # 列下标的全排列
        subscript_list = list(itertools.permutations(tmp_list, shape[0]))
        list_mul = []
        # 不同行不同列元素相乘
        for i in range(len(subscript_list)):
            mul = 1
            for j in range(shape[0]):
                k = subscript_list[i][j]
                mul = self.multiplication(mul, mat[j, k])
            list_mul.append(mul)
        # 相乘后的每一项结果相加
        value = 0
        for m in range(len(list_mul)):
            value = self.addition(value, list_mul[m])
        return value

    #  交换矩阵两行
    @staticmethod
    def exchange_row(mat, row1, row2):
        if row1 == row2:
            return mat
        mat_shape = np.shape(mat)
        for i in range(mat_shape[1]):
            temp = mat[row1, i]
            mat[row1, i] = mat[row2, i]
            mat[row2, i] = temp
        return mat

    #  有限域内的初等行变换矩阵求逆
    def inverse_matrix(self, mat):
        shape = np.shape(mat)
        #  验证是否为方阵
        if shape[0] != shape[1]:
            return None
        #  验证是否可逆
        if self.determinant_value(mat) == 0:
            return None
        #  生成单位阵拼接
        unit_matrix = np.identity(shape[0], dtype=int)
        joint_mat = np.hstack((mat, unit_matrix))
        #  初等变换生成逆矩阵
        for i in range(shape[0]):
            count = i + 1
            while joint_mat[i][i] == 0:
                self.exchange_row(joint_mat, i, count)
                count += 1
            val = joint_mat[i][i]
            for j in range(shape[0] * 2):
                joint_mat[i][j] = self.division(joint_mat[i][j], val)
            for k in range(i + 1, shape[0]):
                val = joint_mat[k][i]
                for l in range(shape[0] * 2):
                    joint_mat[k][l] = self.addition(joint_mat[k][l], self.multiplication(val, joint_mat[i][l]))
        for i in range(1, shape[0]):
            for j in range(i):
                val = joint_mat[j][i]
                for k in range(2 * shape[0]):
                    joint_mat[j][k] = self.addition(joint_mat[j][k], self.multiplication(val, joint_mat[i][k]))
        return joint_mat[:, shape[0]::]


if __name__ == '__main__':
    gf = GF(3)
    Mat1 = np.array([[0, 3, 2, 5], [0, 1, 3, 4], [4, 2, 3, 1], [4, 3, 3, 4]])
    # a = np.array([[1, 2, 3, 4],[2, 3, 4, 5],[3, 4, 5, 6],[4, 5, 6, 7]])
    b = gf.inverse_matrix(Mat1)
    print(gf.matrix_multiplication(Mat1, b))


