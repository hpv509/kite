import numpy as np
import os

def mycentral(Lx, Ly):
    div = [2, 2]
    Ndiv = [Lx / div[0], Ly / div[1]]
    N_unit_cells = Lx * Ly
    imp_posA = np.zeros([1, 2], dtype = 'int')
    imp_posA[0, 0] = Lx / 2
    imp_posA[0, 1] = Ly / 2
    # imp_pos = np.zeros([1, 7], dtype="int")
    # x = Lx / 2
    # y = Ly / 2
    # n = x + y * Lx
    # o = 0

    # imp_pos[0, 0] = 0
    # imp_pos[0, 1] = n
    # imp_pos[0, 2] = x
    # imp_pos[0, 3] = y
    # imp_pos[0, 4] = o

    # imp_pos[0, 5] = x // div[0]
    # imp_pos[0, 6] = y // div[1]

    # new_imp_indexA = []
    # new_imp_indexB = []
    # row = imp_pos[0, :]
    # if row[4] == 0:
    #     new_imp_indexA.append(row[1] % N_unit_cells)
    # if row[4] == 1:
    #     new_imp_indexB.append(row[1] % N_unit_cells)
    return imp_posA
