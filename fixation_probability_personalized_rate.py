import networkx as nx
import random
import numpy as np
import pickle
# import matplotlib.pyplot as plt
from numba import jit
import os
import multiprocessing
# from multiprocessing import Process, Manager
import functools
import time
import math
import scipy.io as scio



# Numerical simulations of the fixation probability of cooperation (rho_c) given
# a network adjacent matrix and individual update rates

@jit(nopython=True)
def single_round(state_array, payoff_array, game_matrix, edge_mat, deg_array):
    """
    each individual plays a single round of game with all its neighbors, and returen te average payoff in this round
    Note: all the inputs variables are numpy array to enable the jit acceleration
    :param state_array: the state of each individual
    :param payoff_array: the payoff of each individual
    :param game_matrix: game matrix for the donation game
    :param edge_mat: the numpy array for edges
    :param deg_array: numpy array of degree for each node
    :return: average payoff for each individual
    """

    for i in range(edge_mat.shape[0]):
        nodex = edge_mat[i, 0]  # iterate each edge on the network
        nodey = edge_mat[i, 1]
        # individuals on both ends of the edge play games and accumulate payoffs
        payoff_array[nodex] += game_matrix[state_array[nodex]][state_array[nodey]]
        payoff_array[nodey] += game_matrix[state_array[nodey]][state_array[nodex]]

    return payoff_array / deg_array  # average payoff for each individual


@jit(nopython=True)
def strategy_update(state_array, payoff_array, nbr_mat, deg_array, nodesnum, w, rate_normalize):
    """
    strategy update with personalized individual update rate
    :param state_array: the state of each individual
    :param payoff_array: the payoff of each individual
    :param nbr_mat: numpy array of neighbors for each individual
    :param deg_array: the numpy array for nodes' degree
    :return: the state array after strategy update

    """

    update_node = rand_pick_list(np.arange(nodesnum),
                                 rate_normalize)  # choose the update node with probability proportional to its update rate
    nbrs_num = deg_array[update_node]  # the number of neighbors of this update node
    nbrs_array = nbr_mat[update_node][:nbrs_num]  # numpy array of the neighbors of the update node

    fitness_group = 1 + w * payoff_array[nbrs_array]  # fitness of the update node's neighbors
    prob_array = fitness_group / np.sum(
        fitness_group)  # probability for the update node imitating strategy from one of its neighbors
    state_array[update_node] = state_array[rand_pick_list(nbrs_array, prob_array)]  # strategy update

    return state_array


@jit(nopython=True)
def evolution(game_matrix, edge_mat, nbrs_mat, deg_array, nodesnum, w, rate_normalize):
    """
    the evolutionary process of cooperation starting from a single cooperator and end when the population reaches
     full cooperation or full defection
    :param game_matrix: game matrix for the donation game
    :param edge_mat: numpy array of network edges
    :param nbrs_mat: numpy array of neighbors for each node
    :param deg_array: numpy array of degree for each node
    :param nodesnum: network size
    :param w: selection intensity
    :param rate_normalize: normalized update rates for individuals with total rate of 1
    :return: frequency of cooperators in a realization of evolutionary process
    """

    total_generation = int(1e10)  # the maximum generation (round) of evolution
    payoff_array = np.zeros(nodesnum, dtype=np.float_)  # array of individual payoff
    state_array = np.zeros(nodesnum, dtype=np.int_)  # array of individual state: 1-->cooperation; 0-->defection
    coop_ini = np.random.choice(nodesnum)
    state_array[coop_ini] = 1  # randomly choose a single node as cooperator

    for time in range(total_generation):
        payoff_array = single_round(state_array, payoff_array, game_matrix, edge_mat,
                                    deg_array)  # all individuals plays around of game
        state_array = strategy_update(state_array, payoff_array,
                                      nbrs_mat, deg_array,
                                      nodesnum, w,
                                      rate_normalize)  # an individual updates its strategy according to its update rate
        payoff_array[:] = 0  # clear the payoff for each individual for the new round of game
        coord = np.sum(state_array)  # count the number of cooperators
        if coord > nodesnum - 1:
            return 1  # end the evolution when system reaches full cooperation
        if coord == 0:
            return 0  # end the evolution when system reaches full defection

    return coord / nodesnum  # return the frequency of cooperators if the system does not reach absorption state within the given maximum generations


@jit(nopython=True)
def process(core, b, edge_mat, nbrs_mat, deg_array, nodesnum, rate_normalize):
    """
    Repeat the simulations to compute the fixation probability of cooperation by calculate the fraction of simulations
    in which the population reaches full cooperation out of realizations reaching absorption.
    :param b: the benefit in the donation game
    :param edge_mat: numpy array of network edges
    :param nbrs_mat: numpy array of neighbors for each node
    :param deg_array: numpy array of degree for each node
    :param nodesnum: network size
    :param rate_normalize: normalized update rates for individuals with total rate of 1
    :return: fixation probability of cooperation
    """

    w = 0.01

    game_matrix = np.zeros((2, 2))  # game matrix of the donation game in the main text
    game_matrix[0][0] = 0  # P defect--defect
    game_matrix[0][1] = b  # T defect-cooperate
    game_matrix[1][0] = -1  # S cooperate-defect
    game_matrix[1][1] = b - 1  # R cooperate-cooperate

    repeat_time = int(1e6)  # total number of realizations
    repeat_array = np.zeros(repeat_time)

    for rep in range(repeat_time):
        freq_c = evolution(game_matrix, edge_mat, nbrs_mat, deg_array, nodesnum, w,
                           rate_normalize)  # frequency of cooperators in a realization of evolutionary process
        repeat_array[rep] = freq_c

    return np.sum(repeat_array == 1) / (np.sum(repeat_array == 1) + np.sum(
        repeat_array == 0))  # only count processes which reach absorption state


@jit(nopython=True)
def rand_pick_list(pick_list, prob_list):
    """
    choose an element according to a given probability distribution, which performs the same as "numpy.random.choice",
    since the function "numpy.random.choice" can not run under numba jit acceleration
    :param pick_list: the element array
    :param prob_list: the array of probability distribution
    :return: the chosen element
    """
    x = random.uniform(0, 1)
    cumulative_probability = 0.0
    for item, item_probability in zip(pick_list, prob_list):
        cumulative_probability += item_probability
        if x <= cumulative_probability:
            break
    return item


def edge_list_array(edge_list):
    """
    convert the list of edges to numpy array to enable jit acceleration
    :param edge_list: the list of edges
    :return: the numpy array (number of edges *2) of edges
    """
    edge_mat = np.zeros([len(edge_list), 2], np.int)
    for i in range(len(edge_list)):
        edge_mat[i, :] = np.array(edge_list[i])

    return edge_mat


def nbr_dict_mat(nbr_dict):
    """
    convert the dict of neighbors to numpy array to enable jit acceleration
    :param nbr_dict: a dict of neighbors (value) of each node (key)
    :return: a n*n numpy array of neighbors
    """
    nodesnum = len(nbr_dict)
    nbr_mat = np.zeros([nodesnum, nodesnum], np.int)
    deg_array = np.zeros(nodesnum, np.int)
    for i, nbrs in nbr_dict.items():
        deg_array[i] = len(nbrs)
        if len(nbrs) > 0:
            nbr_mat[i][:len(nbrs)] = np.array(nbrs)
    return nbr_mat, deg_array


if __name__ == "__main__":
    static_matrix = scio.loadmat("sf_100_k6.mat")["A_sf"]  # import the adjacent matrix from the datatype of mat
    # with open("./result/sf_100_k6_0.pk", 'rb') as f:      # or import the adjacent matrix from a pickle file
    #     static_matrix = pickle.load(f)

    # obtain the numpy array of edges, neighbors, and nodes' degree
    graph = nx.from_numpy_matrix(static_matrix)
    edge_list = list(graph.edges())
    edge_mat = edge_list_array(edge_list)  # numpy array of edges
    nbrs_dict = nx.to_dict_of_lists(graph)
    nbrs_mat, deg_array = nbr_dict_mat(nbrs_dict)  # numpy array of neighbors and degree for each node
    nodesnum = static_matrix.shape[0]  # network size

    # define the individual update rates and normalize to total rate of 1
    rate_array = 1 / deg_array  # inversely proportional to nodes' degree
    # rate_array = deg_array     # proportional to nodes' degree
    # rate_array = np.ones(nodesnum, dtype=np.float_)   # identical update rates
    rate_normalize = np.squeeze(rate_array) / np.sum(np.squeeze(rate_array))

    b_array = [5.2, 5.4, 5.6, 5.8, 6]  # array of benefit b in donation game, example in Fig. 2b
    cpu_cores_num = 10  # 10-cpu core, this number should match to cpu-core running this programme
    rhoc_array = []
    for b_para in b_array:
        core_list = np.arange(cpu_cores_num)

        # start parallel computing the fixation probability of cooperation rho_c
        pool = multiprocessing.Pool()
        t1 = time.time()

        pt = functools.partial(process, b=b_para, edge_mat=edge_mat, nbrs_mat=nbrs_mat, deg_array=deg_array,
                               nodesnum=nodesnum,
                               rate_normalize=rate_normalize)  # perform the process function on each core
        rho_c_list = pool.map(pt, core_list)  # return the fixation probability rho_c computed on each cpu core

        rho_c = sum(rho_c_list) / len(rho_c_list)  # calculate the average rho_c
        rhoc_array.append(rho_c)
        pool.close()
        pool.join()
        t2 = time.time()
        print("b="+str(b_para)+",", "rho_c="+str(rho_c))
        print("Total time:" + (t2 - t1).__str__())

    file = "sf_n100_k6_PersonalizedRate_lambda_1_k.mat"

    # output: array of the benefit b, and array of the corresponding fixation probability of cooperation (rho_c)
    # note: here the cost c=1 (see the game matrix), therefore b/c is equal to b
    scio.savemat(file, {'b_array': b_array,
                        'rhoc_array': rhoc_array})


