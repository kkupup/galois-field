package org.example;

import java.util.Arrays;
import java.util.Date;
import java.util.stream.IntStream;

/**
 * @author WZH
 * @version 1.0
 * @date 2022/3/18 16:57
 */
public class GaloisField {
    private final int omega;
    private final int num;
    private int[] gfilog;
    private int[] gflog;

    public GaloisField(int omega) {
        this.omega = omega;
        this.num = 1 << this.omega;
        this.createGfilog();
        this.createGflog();
    }

    private void createGfilog() {
        this.gfilog = new int[this.num];
        int primitivePoly = this.getPrimitivePolynomial();
        this.gfilog[0] = 1;
        for (int i = 1; i < this.num - 1; i++) {
            int tmp = this.gfilog[i - 1] << 1;
            while (tmp >= this.num) {
                tmp ^= primitivePoly << (Integer.toBinaryString(tmp).length() - Integer.toBinaryString(primitivePoly).length());
            }
            this.gfilog[i] = tmp;
        }
        this.gfilog[this.num - 1] = -1;
    }

    private void createGflog() {
        this.gflog = new int[this.num];
        this.gflog[0] = -1;
        for (int i = 0; i < this.num - 1; i++) {
            this.gflog[this.gfilog[i]] = i;
        }
    }

    private int getPrimitivePolynomial() {
        switch (this.omega) {
            case 3:
                return 11;
            case 4:
                return 19;
            case 8:
                return 285;
            default:
                return 0;
        }
    }

    public int addition(int para1, int para2) {
        return para1 ^ para2;
    }

    public int subtraction(int para1, int para2) {
        return para1 ^ para2;
    }

    public int multiplication(int para1, int para2) {
        if (para1 == 0 || para2 == 0) {
            return 0;
        } else {
            para1 = this.gflog[para1];
            para2 = this.gflog[para2];
            return this.gfilog[(para1 + para2) % (this.num - 1)];
        }
    }

    public int division(int para1, int para2) {
        if (para1 == 0 || para2 == 0) {
            return 0;
        } else {
            para1 = this.gflog[para1];
            para2 = this.gflog[para2];
            int sum = para1 - para2;
            if (sum < 0) {
                sum += (this.num - 1);
            }
            sum %= this.num - 1;
            return this.gfilog[sum];
        }
    }

    public int[][] addMatrix(int[][] mat1, int[][] mat2) {
        if (mat1.length != mat2.length || mat1[0].length != mat2[0].length) {
            return null;
        } else {
            int[][] mat = new int[mat1.length][mat1[0].length];
            for (int i = 0; i < mat1.length; i++) {
                for (int j = 0; j < mat1[0].length; j++) {
                    mat[i][j] = this.addition(mat1[i][j], mat2[i][j]);
                }
            }
            System.out.println(1);
            return mat;
        }
    }

    public int[][] multiplyMatrix(int[][] mat1, int[][] mat2) {
        if (mat1[0].length != mat2.length) {
            return null;
        } else {
            int[][] mat = new int[mat1.length][mat2[0].length];
            for (int i = 0; i < mat1.length; i++) {
                for (int j = 0; j < mat2[0].length; j++) {
                    int number = 0;
                    for (int k = 0; k < mat2.length; k++) {
                        number = this.addition(this.multiplication(mat1[i][k], mat2[k][j]), number);
                    }
                    mat[i][j] = number;
                }
            }
            return mat;
        }
    }

    public int getDeterminantValue(int[][] mat) {
        int rowNum = mat[0].length;
        int colNum = mat.length;
        if (rowNum != colNum) {
            return -1;
        } else {
            int[] tmpList = IntStream.range(0, rowNum).toArray();
            int[][] resultSet = this.getFullPermutation(tmpList);
            int[] mulNums = new int[resultSet.length];
            int i = 0;
            while (i < resultSet.length) {
                int mulNum = 1;
                for (int j = 0; j < rowNum; j++) {
                    int tempColNum = resultSet[i][j];
                    mulNum = this.multiplication(mulNum, mat[j][tempColNum]);
                }
                mulNums[i] = mulNum;
                i++;
            }
            int sumNum = 0;
            for (int mulNum : mulNums) {
                sumNum = this.addition(sumNum, mulNum);
            }
            return sumNum;
        }
    }

    public int[][] getFullPermutation(int[] elements) {
        int[][] resultSet = {{elements[0]}};
        for (int i = 1; i < elements.length; i++) {
            int factorial = 1;
            for (int j = 1; j <= i + 1; j++) {
                factorial *= j;
            }
            int[][] tmp2 = new int[factorial][i + 1];
            for (int j = 0; j < resultSet.length; j++) {
                for (int k = 0; k < i + 1; k++) {
                    int rowNum = (i + 1) * j + k;
                    tmp2[rowNum][k] = elements[i];
                    System.arraycopy(resultSet[j], 0, tmp2[rowNum], 0, k);
                    System.arraycopy(resultSet[j], k, tmp2[rowNum], k + 1, i - k);
                }
            }
            resultSet = tmp2;
        }
        return resultSet;
    }

    public int[][] exchangeRow(int[][] mat, int row1, int row2) {
        if (row1 == row2) {
            return mat;
        } else {
            int[] temp = mat[row1];
            mat[row1] = mat[row2];
            mat[row2] = temp;
        }
        return mat;
    }

    public int[][] getInverseMatrix(int[][] mat) {
        int rowNum = mat[0].length;
        int colNum = mat.length;
        if (rowNum != colNum) {
            return null;
        } else if (this.getDeterminantValue(mat) == 0) {
            return null;
        } else {
            int[][] jointMat = new int[rowNum][2 * colNum];
            for (int i = 0; i < rowNum; i++) {
                for (int j = 0; j < 2 * colNum; j++) {
                    if (j < colNum) {
                        jointMat[i][j] = mat[i][j];
                    } else {
                        if (j - colNum == i) {
                            jointMat[i][j] = 1;
                        } else {
                            jointMat[i][j] = 0;
                        }
                    }
                }
            }
            for (int i = 0; i < rowNum; i++) {
                int mark = i + 1;
                while (jointMat[i][i] == 0) {
                    this.exchangeRow(jointMat, i, mark);
                    mark += 1;
                }
                int elem = jointMat[i][i];
                for (int j = 0; j < 2 * colNum; j++) {
                    jointMat[i][j] = this.division(jointMat[i][j], elem);
                }
                for (int j = i + 1; j < colNum; j++) {
                    elem = jointMat[j][i];
                    for (int k = 0; k < 2 * colNum; k++) {
                        jointMat[j][k] = this.addition(jointMat[j][k], this.multiplication(elem, jointMat[i][k]));
                    }
                }
            }
            for (int i = 1; i < rowNum; i++) {
                for (int j = 0; j < i; j++) {
                    int elem = jointMat[j][i];
                    for (int k = 0; k < 2 * colNum; k++) {
                        jointMat[j][k] = this.addition(jointMat[j][k], this.multiplication(elem, jointMat[i][k]));
                    }
                }
            }
            int[][] inverseMatrix = new int[rowNum][colNum];
            for (int i = 0; i < rowNum; i++) {
                for (int j = 0; j < colNum; j++) {
                    inverseMatrix[i][j] = jointMat[i][j + colNum];
                }
            }
            return inverseMatrix;
        }
    }

    public static void main(String[] args) {
        GaloisField gf = new GaloisField(8);
    }
}
