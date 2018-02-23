"""Assignment 1 {Due Monday, February 5, 2018}: 
    Write a function that given matrix A, 
    will return row echelon form of A

example
-------
    A = [[x, y, z],
         [0, y, z],
         [0, 0, z]]

    

"""

import unittest 
from numpy import array
from pprint import pprint
from time import time
import functools

def row_echelonize(A):
    """puts the matrix A into row echelon form.

    example (matrix A)
    ------------------
         2  1  1      2  1  1
         4 -6  0  ->  0 -8 -2
        -2  7  2      0  0  1

        steps
        -----
        1. Subtract  2 times the first equation from the second -> A[1, :] = A[1, :] -  2 * A[0, :]
        2. Subtract -1 times the first equation from the third  -> A[2, :] = A[2, :] - -1 * A[0, :]
        3. Subtract -1 times the second equation from the third -> A[3, :] = A[3, :] - -1 * A[0, :]

        note
        ----
        looks like the pattern is subtract the (row value/the pivot)*first_now from the row
    """

    try:
        rows, cols = len(A), len(A[0])
    except TypeError:
        return A

    A = array(A, dtype=float)

    row, col = 0, 0
    while row < rows - 1 and col < cols - 1:
        try:
            pivot_row = next(k for k, val in enumerate(A[row:, col], row) if val != 0)
        except StopIteration:
            col += 1
            continue

        if pivot_row > row:
            # swap the rows
            A[pivot_row, :], A[row, :] = A[row, :], A[pivot_row, :].copy()

        pivot = A[row, col]
        for k in range(row + 1, rows):
            front = A[k, col]
            # do the matrix math
            A[k, :] = A[k, :] - front/pivot*A[row, :]

        row, col = row + 1, col + 1

    return A.tolist()

def time_it(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        time_start = time()
        value = func(self, *args, **kwargs)
        time_total = time() - time_start
        print("function: % 30s ran in %0.10f." % (func.__name__, time_total))
        return value
    return wrapper
       
try:
    from c_gaussian import c_row_echelonize
except ImportError:
    c_row_echelonize = row_echelonize

try:
    import matlab.engine
except ImportError:
    pass
    

class RowEchelonTestCase(unittest.TestCase):

    matlab = None

    @classmethod
    def setUp(cls):
        try:
            cls.matlab = matlab.engine.start_matlab()
        except NameError:
            pass

    @time_it
    def row_echelonize(self, A):
        return row_echelonize(A)
    
    @time_it
    def c_row_echelonize(self, A):
        return c_row_echelonize(A)
    
    def matlab_row_echelonize(self, A):
        if self.matlab != None:
            U = self.matlab.row_echelonize(matlab.double(A))
            print("function: % 30s ran in %0.10f." % ("matlab_row_echelonize", U['time']))
            A = U['A']
            try:
                return [i for i in A[0]] if len(A) == 1 else [[i for i in j] for j in A]  
            except TypeError:
                return [A]
        else:
            return row_echelonize(A)
    
    def test_row_echelon_1(self):
        A = [
            [2.0,  1.0,  1.0],
            [4.0, -6.0,  0.0],
            [-2.0, 7.0,  2.0]
        ]

        B = [
            [2,  1,  1],
            [0, -8, -2],
            [0,  0,  1]
        ]
        self.assertEqual(self.row_echelonize(A), B)
        self.assertEqual(self.c_row_echelonize(A), B)
        self.assertEqual(self.matlab_row_echelonize(A), B)
    def test_row_echelon_2(self):
        A = [0]

        self.assertEqual(self.row_echelonize(A), A)
        self.assertEqual(self.c_row_echelonize(A), A)
        self.assertEqual(self.matlab_row_echelonize(A), A)

    def test_row_echelon_3(self):
        A = [0.0, 1.0, 2.0, 3.0, 1.0, 2.0]
        self.assertEqual(self.row_echelonize(A), A)
        self.assertEqual(self.c_row_echelonize(A), A) 
        self.assertEqual(self.matlab_row_echelonize(A), A)

    def test_row_echelon_5(self):
        A = [
            [0.0, 2.0,  1.0,  1.0],
            [0.0, 4.0, -6.0,  0.0],
            [0.0, -2.0, 7.0,  2.0]
        ]

        B = [
            [0, 2,  1,  1],
            [0, 0, -8, -2],
            [0, 0,  0,  1]
        ]

        self.assertEqual(self.row_echelonize(A), B)
        self.assertEqual(self.c_row_echelonize(A), B)
        self.assertEqual(self.matlab_row_echelonize(A), B)


    def test_row_echelon_4(self):
        A = [
            [1, 1, 1],
            [1, 1, 3],
            [2, 5, 8]
        ]
        B = [
            [1, 1, 1],
            [0, 3, 6],
            [0, 0, 2]
        ]
        self.assertEqual(self.row_echelonize(A), B)
        self.assertEqual(self.c_row_echelonize(A), B)
        self.assertEqual(self.matlab_row_echelonize(A), B)


if __name__ == '__main__':
    unittest.main()