import networkx as nx
import random
import numpy as np
import pickle
import matplotlib.pyplot as plt
from numba import jit
import os
import multiprocessing
from multiprocessing import Process, Manager
import functools
import time
import math
import scipy.io as scio

# Numerical simulations of the fixation probability given network adjacent matrix and individual update rates

@jit(nopython=True)
def rand_pick_list(pick_list, prob_list):
    x = random.uniform(0, 1)
    cumulative_probability = 0.0
    for item, item_probability in zip(pick_list, prob_list):
        cumulative_probability += item_probability
        if x <= cumulative_probability:
            break
    return item


def edge_list_array(edge_list):
    edge_mat = np.zeros([len(edge_list), 2], np.int)
    for i in range(len(edge_list)):
        edge_mat[i, :] = np.array(edge_list[i])

    return edge_mat


def nbr_dict_mat(nbr_dict):
    """
    turn nbr_dict to mat
    :param nbr_dict:
    :return:
    """
    nodesnum = len(nbr_dict)
    nbr_mat = np.zeros([nodesnum, nodesnum], np.int)
    deg_array = np.zeros(nodesnum, np.int)
    for i, nbrs in nbr_dict.items():
        deg_array[i] = len(nbrs)
        if len(nbrs) > 0:
            nbr_mat[i][:len(nbrs)] = np.array(nbrs)
    return nbr_mat, deg_array


@jit(nopython=True)
def single_round(state_dict, payoff_array, game_matrix, edge_mat, deg_array):
    """
    play game on each group
    """

    for i in range(edge_mat.shape[0]):
        nodex = edge_mat[i, 0]
        nodey = edge_mat[i, 1]
        payoff_array[nodex] += game_matrix[state_dict[nodex]][state_dict[nodey]]
        payoff_array[nodey] += game_matrix[state_dict[nodey]][state_dict[nodex]]

    return payoff_array / deg_array


@jit(nopython=True)
def replicate_dynamic(state_array, payoff_array, nbr_mat, deg_array, nodesnum, w, rate_normalize):
    """
    replicator dynamic after single round game: DB

    """

    update_node = rand_pick_list(np.arange(nodesnum), rate_normalize)
    nbrs_num = deg_array[update_node]
    nbrs_array = nbr_mat[update_node][:nbrs_num]

    fitness_group = 1 + w * payoff_array[nbrs_array]
    prob_array = fitness_group / np.sum(fitness_group)
    state_array[update_node] = state_array[rand_pick_list(nbrs_array, prob_array)]

    return state_array


@jit(nopython=True)
def evolution(game_matrix, edge_mat, nbrs_mat, deg_array, nodesnum, w, rate_normalize):
    """
    whole process of evolution for 10000 times of generation
    """

    total_generation = int(1e10)
    payoff_array = np.zeros(nodesnum, dtype=np.float_)
    state_array = np.zeros(nodesnum, dtype=np.int_)
    coop_ini = np.random.choice(nodesnum)
    state_array[coop_ini] = 1

    for time in range(total_generation):
        payoff_array = single_round(state_array, payoff_array, game_matrix, edge_mat, deg_array)
        state_array = replicate_dynamic(state_array, payoff_array, nbrs_mat, deg_array, nodesnum, w, rate_normalize)
        payoff_array[:] = 0
        coord = np.sum(state_array)
        if coord > nodesnum - 1:
            return 1
        if coord == 0:
            return 0

    return coord / nodesnum


@jit(nopython=True)
def process(core, b, edge_mat, nbrs_mat, deg_array, nodesnum, rate_normalize):
    w = 0.01

    game_matrix = np.zeros((2, 2))
    game_matrix[0][0] = 0  # P defect--defect
    game_matrix[0][1] = b  # T d-c
    game_matrix[1][0] = -1  # S
    game_matrix[1][1] = b - 1  # R

    repeat_time = int(1e6)
    repeat_array = np.zeros(repeat_time)

    for rep in range(repeat_time):
        coord_freq = evolution(game_matrix, edge_mat, nbrs_mat, deg_array, nodesnum, w, rate_normalize)
        repeat_array[rep] = coord_freq

    return np.sum(repeat_array == 1) / (np.sum(repeat_array == 1) + np.sum(repeat_array == 0))


if __name__ == "__main__":
    # with open("./result/sf_100_k6_0.pk", 'rb') as f:
    #     static_matrix = pickle.load(f)
    static_matrix = scio.loadmat("sf_100_k6.mat")["A_sf"]

    # ---------------------------  graph construction  --------------------------
    graph = nx.from_numpy_matrix(static_matrix)
    edge_list = list(graph.edges())
    edge_mat = edge_list_array(edge_list)
    nbrs_dict = nx.to_dict_of_lists(graph)
    nbrs_mat, deg_array = nbr_dict_mat(nbrs_dict)
    nodesnum = 100

    rate_array = 1 / deg_array
    # rate_array = deg_array
    # rate_array = np.ones(static_matrix.shape[1], dtype=np.float_)
    rate_normalize = np.squeeze(rate_array) / np.sum(np.squeeze(rate_array))

    b_para_list = [5.2, 5.4, 5.6, 5.8, 6]
    for b_para in b_para_list:
        core_list = np.arange(20)  # 64-cpu core

        pool = multiprocessing.Pool()
        t1 = time.time()

        pt = functools.partial(process, b=b_para, edge_mat=edge_mat, nbrs_mat=nbrs_mat, deg_array=deg_array,
                               nodesnum=nodesnum, rate_normalize=rate_normalize)
        coor_freq_list = pool.map(pt, core_list)

        coor_freq_core = sum(coor_freq_list) / len(coor_freq_list)

        pool.close()
        pool.join()
        t2 = time.time()
        print("Total time:" + (t2 - t1).__str__())
        print((b_para, coor_freq_core))
        # file = "./result/sf_100_k6_0_static_db_heterRate_uniMutant_1_w_b" + str(b_para) + "_1.pk"
        #
        # with open(file, 'wb') as f:
        #     pickle.dump(coor_freq_core, f)
