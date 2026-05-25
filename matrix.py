from fractions import Fraction
from copy import deepcopy

class Matrix:
    def __init__(self, matrix: object, augment: object = None) -> None:
        # self.matrix = matrix
        self.rows = len(matrix)
        self.cols = len(matrix[0])
        self.matrix = [[Fraction(i) for i in _] for _ in matrix]
        self.augmentVector = augment
        if augment is None:
            self.isAugmented = False
        else:
            self.isAugmented = True

    def print_matrix(self):
        for rowNumber, row in enumerate(self.matrix, start=1):
            for colNumber, col in enumerate(row, start=1):
                print(f'{f"| {col}" if colNumber == self.cols + 1 else col}', end=' ')
            print()

    def augment(self, vector: list):
        if len(vector) != self.rows:
            print(f'Vector length {len(vector)} doesn\'t match number of columns {self.rows}')
            return

        self.augmentVector = vector

        for rowNumber, row in enumerate(self.matrix, start=1):
            row.append(self.augmentVector[rowNumber - 1])
        else:
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
        resultantMatrix = deepcopy(self.matrix)

        if source > self.rows or target > self.rows or source < 0 or target < 0:
            print('Invalid source/target rows')
            return resultantMatrix

        self.matrix[source - 1], self.matrix[target - 1] = self.matrix[target - 1], self.matrix[source - 1]
        return self.matrix

    def rowOperation(self, sourceRow, targetRow, constant = 1):
        # matrix[source] = matrix[source] + constant*matrix[target]
        for i in range(self.cols):
            self.matrix[sourceRow - 1][i] += constant*self.matrix[targetRow - 1][i]

        return self.matrix

    def rowScale(self, targetRow, scalingFactor):
        # matrix[target] = scalingFactor*matrix[target]

        for i in range(self.cols):
            self.matrix[targetRow - 1][i] *= scalingFactor

        return self.matrix

    def getRowEchelonForm(self):

        lead = 0
        row = 0

        while lead < self.cols and row < self.rows:

            while lead < self.cols and self.matrix[row][lead] == 0:
                for k in range(row + 1, self.rows):
                    if self.matrix[k][lead] != 0:
                        self.swapRows(row + 1, k + 1)
                        break
                else:
                    # leading entries aren't diagonally aligned
                    lead += 1

            if lead == self.cols:
                break

            for k in range(row + 1, self.rows):
                if self.matrix[k][lead] != 0:
                    self.rowOperation(k + 1, row + 1, -(self.matrix[k][lead]/self.matrix[row][lead]))

            row += 1
            lead += 1

    def getReducedRowEchelonForm(self):

        # First get Row Echelon Form
        self.getRowEchelonForm()

        # Work from bottom row upward
        for row in range(self.rows - 1, -1, -1):

            # Find pivot column
            lead = -1
            for col in range(self.cols):
                if self.matrix[row][col] != 0:
                    lead = col
                    break

            # Skip zero rows
            if lead == -1:
                continue

            # Scale row so pivot becomes 1
            pivot = self.matrix[row][lead]
            self.rowScale(row + 1, 1 / pivot)

            # Eliminate entries above pivot
            for k in range(row):
                if self.matrix[k][lead] != 0:
                    self.rowOperation(
                        k + 1,
                        row + 1,
                        -self.matrix[k][lead]
                    )

        return self.matrix

    def getPivotColumns(self):
        pivots = []

        for row in range(self.rows):
            for col in range(self.cols):
                if self.matrix[row][col] == 1:
                    # check if this is a leading 1
                    if all(self.matrix[row][k] == 0 for k in range(col)):
                        pivots.append(col)
                        break

        return pivots

    def getNonPivotColumns(self):
        pivots = self.getPivotColumns()

        return [col for col in range(self.cols) if col not in pivots]

    def nullSpaceSolutions(self):
        self.getReducedRowEchelonForm()

        pivots = self.getPivotColumns()
        freeCols = self.getNonPivotColumns()

        solutions = []

        for free in freeCols:

            vec = [0] * self.cols
            vec[free] = 1

            for row, pivotCol in enumerate(pivots):
                vec[pivotCol] = -self.matrix[row][free]

            solutions.append(vec)

        return solutions





def main():
    A1 = [
        [1, 2],
        [3, 4]
    ]

    A2 = [
        [2, 1, 3],
        [4, 1, 6],
        [2, 0, 2]
    ]

    A3 = [
        [1, 0, 2],
        [0, 1, 3],
        [0, 0, 1]
    ]

    B1 = [
        [1, 2, 3],
        [0, 1, 4],
        [0, 0, 1]
    ]

    B2 = [
        [1, 0, 0, 5],
        [0, 1, 0, 6],
        [0, 0, 1, 7]
    ]

    C1 = [
        [1, 2, 3],
        [0, 0, 0],
        [0, 0, 0]
    ]

    C2 = [
        [0, 0, 0],
        [1, 2, 3],
        [0, 0, 0]
    ]

    D1 = [
        [1, 2, 3],
        [2, 4, 6],
        [3, 6, 9]
    ]

    D2 = [
        [1, 2, 3, 4],
        [2, 4, 6, 8],
        [1, 1, 1, 1]
    ]

    E1 = [
        [0, 2, 1],
        [1, 1, 0],
        [2, 3, 4]
    ]

    E2 = [
        [0, 0, 1],
        [0, 2, 3],
        [1, 0, 0]
    ]

    F1 = [
        [-1, 2, -3],
        [2, -4, 6],
        [-3, 6, -9]
    ]

    F2 = [
        [0, -2, 1],
        [-1, 3, -4],
        [2, -1, 5]
    ]

    G1 = [
        [1, 2, 0, 1, 3],
        [2, 4, 1, 3, 7],
        [1, 2, 1, 2, 4]
    ]

    G2 = [
        [1, 2, 3],
        [2, 4, 7],
        [1, 1, 1],
        [3, 5, 9],
        [2, 3, 4]
    ]

    H1 = [
        [1, 3, 2],
        [2, 6, 5],
        [1, 3, 4]
    ]

    H2 = [
        [2, 4, 8],
        [3, 6, 12],
        [1, 2, 3]
    ]

    I1 = [
        [1, 2, 3],
        [1, 2, 3],
        [1, 2, 3]
    ]

    I2 = [
        [2, 1, 0],
        [2, 1, 0],
        [4, 2, 0]
    ]

    J1 = [
        [0, 0, 1, 2],
        [0, 1, 0, 3],
        [1, 0, 0, 4],
        [2, 0, 0, 5]
    ]

    J2 = [
        [0, 0, 0, 1],
        [0, 2, 3, 4],
        [1, 0, 5, 6],
        [0, 1, 0, 0]
    ]

    X =[
        [1, 0, 2, -1],
        [0, 1, 3, 4],
        [0, 0, 0, 0]
    ]

    listOfMatrices = [X]
    # listOfMatrices = [A1, A2, A3, B1, B2, C1, C2, D1, D2, E1, E2, F1, F2, G1, G2, H1, H2, I1, I2, J1, J2, X]

    for mat in listOfMatrices:
        matrix = Matrix(mat)

        print('Matrix')
        matrix.print_matrix()
        print('-'*15)
        # augmentVector = [11, 12]
        # matrix.augment(augmentVector)

        print('REF:')
        matrix.getRowEchelonForm()

        matrix.print_matrix()
        print('-'*15)

        print('RREF')
        matrix.getReducedRowEchelonForm()

        matrix.print_matrix()

        print('-'*15)

        print("Pivot cols:")
        print(matrix.getPivotColumns())

        print('-' * 15)

        print("Non-pivot cols:")
        print(matrix.getNonPivotColumns())

        print('-' * 15)

        print("null space solutions:")
        print(matrix.nullSpaceSolutions())

        print('-' * 15 + 'END')



if __name__ == '__main__':
    main()

