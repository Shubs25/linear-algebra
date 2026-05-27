from fractions import Fraction
from copy import deepcopy
from math import sqrt

class Matrix:
    def __init__(self, matrix: object, augment: Matrix = None) -> None:
        # self.matrix = matrix
        self.rows = len(matrix)
        self.cols = len(matrix[0])
        self.matrix = [[Fraction(i) for i in _] for _ in matrix]
        self.augmentVector = augment
        if augment is None:
            self.isAugmented = False
        else:
            self.isAugmented = True

    def augment(self, vector: Matrix):
        if len(vector.matrix) != self.rows:
            print(f'Vector length {len(vector.matrix)} doesn\'t match number of columns {self.rows}')
            return

        self.augmentVector = vector
        self.isAugmented = True

    def transpose(self):
        if self.isAugmented:
            print("Transpose not allowed on Augmented Matrix")
            return self.matrix

        transposedMatrix = [[0 for _ in range(self.rows)] for _ in range(self.cols)]

        for i in range(self.cols):
            for j in range(self.rows):
                transposedMatrix[i][j] = self.matrix[j][i]

        self.matrix = transposedMatrix
        self.rows, self.cols = self.cols, self.rows
        return self.matrix

    def swapRows(self, source, target):
        if source > self.rows or target > self.rows or source <= 0 or target <= 0:
            print('Invalid source/target rows')
            return self

        self.matrix[source - 1], self.matrix[target - 1] = self.matrix[target - 1], self.matrix[source - 1]
        return self

    def rowOperation(self, sourceRow, targetRow, constant = 1):
        # matrix[source] = matrix[source] + constant*matrix[target]
        for i in range(self.cols):
            self.matrix[sourceRow - 1][i] += constant*(self.matrix[targetRow - 1][i])

        return self

    def rowScale(self, targetRow, scalingFactor):
        # matrix[target] = scalingFactor*matrix[target]

        for i in range(self.cols):
            self.matrix[targetRow - 1][i] *= scalingFactor

        return self

    def getRowEchelonForm(self):
        resultantMatrix = deepcopy(self)
        
        lead = 0
        row = 0

        while lead < resultantMatrix.cols and row < resultantMatrix.rows:

            while lead < resultantMatrix.cols and resultantMatrix.matrix[row][lead] == 0:
                for k in range(row + 1, resultantMatrix.rows):
                    if resultantMatrix.matrix[k][lead] != 0:
                        resultantMatrix.swapRows(row + 1, k + 1)
                        if resultantMatrix.isAugmented:
                            resultantMatrix.augmentVector.swapRows(row + 1, k + 1)
                        break
                else:
                    # leading entries aren't diagonally aligned
                    lead += 1

            if lead == resultantMatrix.cols:
                break

            for k in range(row + 1, resultantMatrix.rows):
                if resultantMatrix.matrix[k][lead] != 0:
                    if resultantMatrix.isAugmented:
                        resultantMatrix.augmentVector.rowOperation(k + 1, row + 1, -(resultantMatrix.matrix[k][lead]/resultantMatrix.matrix[row][lead]))
                    resultantMatrix.rowOperation(k + 1, row + 1, -(resultantMatrix.matrix[k][lead]/resultantMatrix.matrix[row][lead]))


            row += 1
            lead += 1

        return resultantMatrix


    def getReducedRowEchelonForm(self):

        # First get Row Echelon Form
        resultantMatrix = self.getRowEchelonForm()

        # Work from bottom row upward
        for row in range(resultantMatrix.rows - 1, -1, -1):

            # Find pivot column
            lead = -1
            for col in range(resultantMatrix.cols):
                if resultantMatrix.matrix[row][col] != 0:
                    lead = col
                    break

            # Skip zero rows
            if lead == -1:
                continue

            # Scale row so pivot becomes 1
            pivot = resultantMatrix.matrix[row][lead]
            resultantMatrix.rowScale(row + 1, 1 / pivot)
            if resultantMatrix.isAugmented:
                resultantMatrix.augmentVector.rowScale(row + 1, 1 / pivot)

            # Eliminate entries above pivot
            for k in range(row):
                if resultantMatrix.matrix[k][lead] != 0:
                    if resultantMatrix.isAugmented:
                        resultantMatrix.augmentVector.rowOperation(
                            k + 1,
                            row + 1,
                            -resultantMatrix.matrix[k][lead]
                        )
                    resultantMatrix.rowOperation(
                        k + 1,
                        row + 1,
                        -resultantMatrix.matrix[k][lead]
                    )


        return resultantMatrix

    def getPivotColumns(self):
        pivots = []

        for row in range(self.rows):
            for col in range(self.cols):
                if self.matrix[row][col] == 1:
                    # check if this is a leading 1
                    if all(self.matrix[row][k] == 0 for k in range(col)):
                        if all(
                                self.matrix[r][col] == 0
                                for r in range(self.rows)
                                if r != row
                        ):
                            pivots.append(col + 1)
                            break

        return pivots

    def getNonPivotColumns(self):
        pivots = self.getPivotColumns()

        return [col for col in range(1, self.cols + 1) if col not in pivots]

    def getParticularSolution(self):
        b = self.augmentVector

        rref = self.getReducedRowEchelonForm()

        pivots = rref.getPivotColumns()

        n = rref.cols
        x = [0] * n

        # free variables = 0 already (implicit)

        # compute pivot variables
        for row, pivotCol in enumerate(pivots):
            val = b.matrix[row][0]

            # subtract contributions of free variables (all set to 0 → no effect)
            x[pivotCol - 1] = val

        return x

    def nullSpaceSolution(self):
        resultantMatrix = self.getReducedRowEchelonForm()

        pivots = resultantMatrix.getPivotColumns()
        freeCols = resultantMatrix.getNonPivotColumns()

        solutions = []

        for free in freeCols:

            vec = [0] * resultantMatrix.cols
            vec[free - 1] = 1

            for row, pivotCol in enumerate(pivots):
                vec[pivotCol - 1] = -resultantMatrix.matrix[row][free - 1]

            solutions.append(vec)

        return solutions

    def isConsistent(self):

        rows = self.rows
        cols = self.cols

        for row in range(rows):

            allZero = all(
                self.matrix[row][col] == 0
                for col in range(cols)
            )

            if allZero and self.augmentVector.matrix[row][0] != 0:
                return False

        return True


    def getGeneralSolution(self, coeffs = None):

        if not self.isConsistent():
            return ["inconsistent System"]

        x_p = self.getParticularSolution()
        nullspace = self.nullSpaceSolution()

        n = len(x_p)

        # return a generic symbolic representation
        if coeffs is None:
            return {
                "particular sol": x_p,
                "nullspace sol": nullspace,
                "general_form": "x = x_p + Σ c_i v_i"
            }

        # numeric evaluation
        x = x_p[:]

        for c, v in zip(coeffs, nullspace):
            for i in range(n):
                x[i] += c * v[i]

        return x

    def getLUWithElemntaryMatrices(self):
        n = self.rows

        A = [[self.matrix[i][j] for j in range(n)] for i in range(n)]
        U = [[A[i][j] for j in range(n)] for i in range(n)]
        L = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]

        elementary_matrices = []

        for k in range(n):

            if U[k][k] == 0:
                raise ZeroDivisionError("Zero pivot encountered (no pivoting allowed for SPD case assumption).")

            for i in range(k + 1, n):

                factor = U[i][k] / U[k][k]

                # construct elementary matrix E
                E = [[Fraction(1 if r == c else 0) for c in range(n)] for r in range(n)]
                E[i][k] = -factor

                elementary_matrices.append(Matrix(E))

                # apply row operation to U
                for j in range(n):
                    U[i][j] -= factor * U[k][j]

                # update L using inverse effect
                L[i][k] += factor

        return Matrix(L), Matrix(U), elementary_matrices

    def cholesky(self):
        n = self.rows

        # L initialized with zeros
        L = [[Fraction(0) for _ in range(n)] for _ in range(n)]

        A = [[self.matrix[i][j] for j in range(n)] for i in range(n)]

        for i in range(n):
            for j in range(i + 1):

                sum_ = sum(L[i][k] * L[j][k] for k in range(j))

                if i == j:
                    val = A[i][i] - sum_

                    if val <= 0:
                        raise ValueError("Matrix is not positive definite")

                    L[i][j] = Fraction(sqrt(val))

                else:
                    if L[j][j] == 0:
                        raise ZeroDivisionError("Zero diagonal encountered")

                    L[i][j] = (A[i][j] - sum_) / L[j][j]

        return Matrix(L)

    def QRDecomposition(self):
        A = self.matrix
        m = self.rows
        n = self.cols

        # Q as list of orthonormal vectors
        Q = [[0 for _ in range(n)] for _ in range(m)]
        # R matrix
        R = [[0 for _ in range(n)] for _ in range(n)]

        # store orthonormal vectors
        q_vectors = []

        for j in range(n):
            v = [A[i][j] for i in range(m)]  # column j

            # Gram-Schmidt projection
            for i in range(j):
                R[i][j] = dot(q_vectors[i], v)
                proj = scalar_mult(q_vectors[i], R[i][j])
                v = subtract(v, proj)

            R[j][j] = norm(v)

            if R[j][j] == 0:
                raise ValueError("Columns are linearly dependent")

            q = [vi / R[j][j] for vi in v]
            q_vectors.append(q)

        # build Q matrix (columns = q_vectors)
        for j in range(n):
            for i in range(m):
                Q[i][j] = q_vectors[j][i]

        return Matrix(Q), Matrix(R)


    def __str__(self):
        resultantString = ''

        for rowNumber, row in enumerate(self.matrix):
            resultantString += ' '.join(map(str, row))
            if self.isAugmented:
                resultantString += f' | {self.augmentVector.matrix[rowNumber][0]}'
            resultantString += '\n'

        return resultantString

    def __mul__(self, other: Matrix):
        m = self.rows
        n = self.cols
        p = other.cols

        # sanity check
        assert other.rows == n, "Incompatible dimensions"

        result = [[Fraction(0) for _ in range(p)] for _ in range(m)]

        for i in range(m):
            for j in range(p):
                s = Fraction(0)
                for k in range(n):
                    s += Fraction(self.matrix[i][k]) * Fraction(other.matrix[k][j])
                result[i][j] = s

        return result


# TODO: these methods can later be encapsulated in a new vector subclass
def dot(u, v):
    return sum(ui * vi for ui, vi in zip(u, v))

def norm(v):
    return sqrt(float(dot(v, v)))

def scalar_mult(v, s):
    return [s * vi for vi in v]

def subtract(u, v):
    return [ui - vi for ui, vi in zip(u, v)]


def main():

    # X =[
    #     [1, 2, 0, 1, 3, 0, 2, 1, 4],
    #     [2, 4, 1, 3, 6, 1, 5, 2, 8],
    #     [1, 2, 1, 2, 3, 1, 3, 1, 5],
    #     [3, 6, 1, 4, 9, 1, 7, 3, 12],
    #     [0, 0, 1, 1, 0, 1, 1, 0, 1],
    #     [1, 2, 2, 3, 3, 2, 4, 1, 6]
    # ]
    #
    # augmentX = [
    #     [1],
    #     [2],
    #     [1],
    #     [3],
    #     [0],
    #     [1]
    # ]

    A = [
        [1, 2, 3, 1, 0, 2, 4, 5, 6],
        [2, 4, 6, 2, 1, 4, 8, 10, 12],
        [1, 1, 1, 0, 1, 1, 2, 2, 3],
        [3, 6, 9, 3, 1, 6, 12, 15, 18],
        [0, 1, 1, 0, 1, 1, 1, 2, 2],
        [1, 3, 4, 1, 2, 3, 5, 7, 8]
    ]

    b = [
        [1],
        [2],
        [1],
        [3],
        [1],
        [2]
    ]

    listOfMatrices = [A]
    # listOfMatrices = [A1, A2, A3, B1, B2, C1, C2, D1, D2, E1, E2, F1, F2, G1, G2, H1, H2, I1, I2, J1, J2, X]


    for mat in listOfMatrices:
        matrix = Matrix(mat)
        augmentMatrix = Matrix(b)
        matrix.augment(augmentMatrix)

        print('Matrix')
        print(matrix)
        print('-'*15)
        # augmentVector = [11, 12]
        # matrix.augment(augmentVector)

        print('REF:')
        refMatrix = matrix.getRowEchelonForm()
        print(refMatrix)
        print('-'*15)

        print('RREF')
        rrefMatrix = refMatrix.getReducedRowEchelonForm()
        print(rrefMatrix)
        print('-'*15)

        print("Pivot cols:")
        print(rrefMatrix.getPivotColumns())
        print('-' * 15)

        print("Non-pivot cols:")
        print(rrefMatrix.getNonPivotColumns())
        print('-' * 15)

        print("null space solutions:")
        print(rrefMatrix.nullSpaceSolution())
        print('-' * 15)

        print("Particular solution:")
        print(rrefMatrix.getParticularSolution())
        print('-' * 15)

        print("General solution:")
        print(rrefMatrix.getGeneralSolution())
        print('-' * 15)


        print('-' * 15 + 'END')

def main2():
    matrixA = Matrix([
        [4, 1, 2],
        [1, 3, 0],
        [2, 0, 5],
    ])
    # matrixA = Matrix([
    #     [1, 2, 0, 1, 3],
    #     [0, 1, 4, 2, 1],
    #     [2, 0, 1, 3, 5],
    #     [1, 1, 1, 0, 2],
    #     [3, 2, 5, 1, 0]
    # ])

    for stuff in matrixA.getLUWithElemntaryMatrices():
        if isinstance(stuff, list):
            print('Elementary matrices: ' + '-'*10)
            for mat in stuff:
                print(mat)
            print('-'*15)
        else:
            print(stuff)

    print('*--------------------cholesky----------------*')
    print(matrixA.cholesky())
    print('*--------------------cholesky----------------*')
    print('*--------------------QR----------------*')
    for mat in matrixA.QRDecomposition():
        print(mat)
    print('*--------------------QR----------------*')



if __name__ == '__main__':
    main2()

