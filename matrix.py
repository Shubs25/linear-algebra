from fractions import Fraction
from copy import deepcopy

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

    def print_matrix(self):
        for rowNumber, row in enumerate(self.matrix, start=1):
            for colNumber, col in enumerate(row, start=1):
                print(f'{f"| {col}" if colNumber == self.cols + 1 else col}', end=' ')
            print()

    def augment(self, vector: Matrix):
        if len(vector.matrix) != self.rows:
            print(f'Vector length {len(vector.matrix)} doesn\'t match number of columns {self.rows}')
            return

        self.augmentVector = vector
        self.isAugmented = True

        # for rowNumber, row in enumerate(self.matrix, start=1):
        #     row.append(self.augmentVector[rowNumber - 1])
        # else:
        #     self.isAugmented = True


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
                        pivots.append(col + 1)
                        break

        return pivots

    def getNonPivotColumns(self):
        pivots = self.getPivotColumns()

        return [col + 1 for col in range(self.cols) if col not in pivots]

    def getParticularSolution(self, nullSpaceSolution: bool = False):
        # b is a vector

        b = Matrix([[0] for _ in range(self.rows)])

        if not nullSpaceSolution:
            b = self.augmentVector

        resultantMatrix = self.getReducedRowEchelonForm()

        pivots = resultantMatrix.getPivotColumns()
        freeCols = resultantMatrix.getNonPivotColumns()

        solutions = []

        for free in freeCols:

            vec = [0] * resultantMatrix.cols
            vec[free - 1] = 1

            for row, pivotCol in enumerate(pivots):
                vec[pivotCol - 1] = b.matrix[row][0] - resultantMatrix.matrix[row][free - 1]

            solutions.append(vec)

        return solutions

    def nullSpaceSolution(self):
        nullMatrix = Matrix([[0] for _ in range(self.rows)])

        return self.getParticularSolution(True)

    def __str__(self):
        resultantString = ''
        for rowNumber, row in enumerate(self.matrix, start=1):
            for colNumber, col in enumerate(row, start=1):
                resultantString += f'{f"| {col}" if colNumber == self.cols + 1 else col} '
            if self.isAugmented:
                resultantString += f'| {self.augmentVector.matrix[rowNumber - 1][0]}'
            resultantString += '\n'

        return resultantString



def main():
    '''
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
    '''

    X =[
        [1, 2, 0, 1, 3, 0, 2, 1, 4],
        [2, 4, 1, 3, 6, 1, 5, 2, 8],
        [1, 2, 1, 2, 3, 1, 3, 1, 5],
        [3, 6, 1, 4, 9, 1, 7, 3, 12],
        [0, 0, 1, 1, 0, 1, 1, 0, 1],
        [1, 2, 2, 3, 3, 2, 4, 1, 6]
    ]


    listOfMatrices = [X]
    # listOfMatrices = [A1, A2, A3, B1, B2, C1, C2, D1, D2, E1, E2, F1, F2, G1, G2, H1, H2, I1, I2, J1, J2, X]
    augmentX = [
        [1],
        [2],
        [1],
        [3],
        [0],
        [1]
    ]

    for mat in listOfMatrices:
        matrix = Matrix(mat)
        augmentMatrix = Matrix(augmentX)
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


        print('-' * 15 + 'END')



if __name__ == '__main__':
    main()

