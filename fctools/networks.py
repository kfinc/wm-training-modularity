# -*- coding: utf-8 -*-


"""
Created on Tue Jun 12 2018
Last edit: Tue Jun 12 2018
@author: kfinc

"""

import numpy as np
import pandas as pd

def calculate_lsn_edges (A, labels):
    """Function calculates number of edges between and within predefined large-scale networks (LSNs).
    The function takes binary symetrical adjacency matrix, module assignment of each ROI and calculate number of edges between and within each
    large-scale network.

    Parameters
    ------------
    array: N x N binary ajdacency matrix
    array: N-length vector with module assignment for each node

    Returns
    ------------
    array: M x M matrix with number of edges between each module

    """
    columns = np.unique(labels)
    lsn_matrix = np.zeros((len(labels), len(columns)))
    lsn_edges = np.zeros((len(columns), len(columns)))

    for col in range(len(columns)):
        module = columns[col, ]
        for row in range(len(labels)):
            if (labels[row, ] == module):
                lsn_matrix[row, col] = 1
            else:
                lsn_matrix[row, col] = 0

    lsn_edges = lsn_matrix.T @ A @ lsn_matrix
    return lsn_edges


def allegiance_matrix(M):
    """Calculates n x n allegiance matrix P from n x t matrix M, where n represents nodes, and t represents time windows. 
    Each value on allegiance matrix P represents the probability that node i and node j have been assigned to the same 
    community (Bassett et al., 2014; Mattar et al., 2015).
    
    Parameters
    ------------
    array: n x t module assignment matrix t

    Returns
    ------------
    array: n x n allegiance matrix P """
    
    n_nodes = M.shape[0]
    n_slices = M.shape[1]
    T = np.zeros((n_nodes, n_nodes))
    
    for i in range(n_nodes):
        for j in range(i):
            if i == j:
                continue
            else:
                t = np.sum(M[i, :] == M[j, :])
                T[i, j] = t
    
    P = (T + T.T)/n_slices
    np.fill_diagonal(P, 1)
    
    return P


def allegiance_matrices_4d(M):
    """Calculates 4D array composed of allegiance matrices for each subject and each condition/session.
    
    Parameters
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x M(window) 

    Returns
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x N(node), where N x N is allegiance matrix 
        
    """
    n_roi = len(M[0,0,:,0])  
    n_sub = len(M[:,0,0,0])
    n_ses = len(M[0,:,0,0])
    AM = np.zeros((n_sub, n_ses, n_roi, n_roi))
    
    for sub in range(n_sub):
        for ses in range(n_ses):
            AM[sub, ses, :, :] = allegiance_matrix(M[sub, ses, :, :])
    return AM    
    
def allegiance_matrix_opti(M):
    """Calculates n x n allegiance matrix P from o x n x t matrix M, where o represents module community detection 
    optimizations, n represents nodes, and t represents time windows. Each value on allegiance matrix PO represents the 
    probability that node i and node j have been assigned to the same 
    community (Bassett et al., 2014; Mattar et al., 2015).
    
    Parameters
    ------------
    array: o x n x t module assignment matrix t

    Returns
    ------------
    array: n x n allegiance matrix P (mean across all optimizations)"""
    
    n_optimizations = M.shape[0]
    n_nodes = M.shape[1]
    P = np.zeros((n_optimizations, n_nodes, n_nodes))
    PO = np.zeros((n_nodes, n_nodes))
    
    for i in range(n_optimizations):
        P[i, :, :] = allegiance_matrix(M[i, :, :])
        PO = P.mean(axis = 0)
    return PO

def allegiance_matrix_opti_5d(M):
    
    n_subs = M.shape[0]
    n_sess = M.shape[1]
    n_optimizations = M.shape[2]
    n_nodes = M.shape[3]
    
    P = np.zeros((n_subs, n_sess, n_optimizations, n_nodes, n_nodes))
    
    for i in range(n_subs):
        for j in range(n_sess):
            P[i, j, :, :, :] = allegiance_matrix_opti_3d(M[i, j, :, :, :])
    return P

def sort_matrices_4d(M, idx):
    """Sorts matrices according to predefinded index.
    
    Parameters
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x N(node) (unsorted)
    array: N-length vector with index to sort matrix

    Returns
    ------------
    array: 4D array S(subjects) x C(condition/session) x N(node) x N(node) (sorted)"""
    M1 = M[:,:,:,idx]
    M2 = M1[:,:,idx,:]
    return M2


def dumming(labels):
    """Gennerate vectors with dummy variables from single vector with categorical variable
    
    Parameters
    ------------
    array: np.array, N-length vector with categorical variable

    Returns
    ------------
    array: 3D array with M rows and N columns with binary variables"""

    columns = np.unique(labels)
    dummies = np.zeros((len(labels), len(columns)))

    for col in range(len(columns)):
        module = columns[col]
        for row in range(len(labels)):
            if (labels[row] == module):
                dummies[row, col] = 1
            else:
                dummies[row, col] = 0
    return dummies


def dumming_pd(labels):
    """Gennerate vectors with dummy variables from single vector with categorical variable
    
    Parameters
    ------------
    array: np.array, N-length vector with categorical variable

    Returns
    ------------
    array: pd.DataFrame with M rows and N columns with binary variables
    
    """
    columns = np.unique(labels)
    return pd.DataFrame(dumming(labels), columns=columns)

def upper_tri_masking(A):
    """Getting values of upper triangle of matrix without diagonal"""
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:,None] < r
    return A[mask]


def fc_cartography(M, modules):
    """Function which calculates mean integration and recruitment values from sorted allegiance matrices"""
    dummy_networks = dumming(sorted(modules))
    roi_n = len(dummy_networks)
    net = np.size(dummy_networks,1)
    diagnostics = np.zeros((net, net))

    for i in range(net):
        for j in range(net):
            vec1 = dummy_networks[:,i].astype('bool')
            vec2 = dummy_networks[:,j].astype('bool')

            L = M[vec1,:]
            P = L[:,vec2]

            #if i == j:
             #   m = upper_tri_masking(P).mean()
            #else:
            m = P.mean()

            diagnostics[i, j] = m
    return diagnostics


'''def fc_cartography(M, modules):
    """Function which calculates mean integration and recruitment values from sorted allegiance matrices"""
    dummy_networks = dumming(sorted(modules))
    n_nod = len(dummy_networks)
    n_net = np.size(dummy_networks,1)
    diagnostics = np.zeros((n_net, n_net))

    for i in range(n_net):
        for j in range(n_net):
            vec1 = dummy_networks[:,i].astype('bool')
            vec2 = dummy_networks[:,j].astype('bool')

            L = M[vec1,:]
            P = L[:,vec2]

            #if i == j:
            #    m = upper_tri_masking(P).mean()
            #else:
            m = P.mean()

            diagnostics[i, j] = m
    return diagnostics
    '''



def fc_cartography_4d(M, modules):
    """Function which calculates mean integration and recruitment values from sorted allegiance matrices. 4D version"""
    sub_n = len(M[:,0,0,0])
    ses_n = len(M[0,:,0,0])
    mod_n = len(np.unique(modules))

    fc_cart = np.zeros((sub_n, ses_n, mod_n, mod_n))

    for i in range(sub_n):
        for j in range(ses_n):
            X = M[i,j,:,:]
            fc_cart[i, j, :, :] = fc_cartography(X, modules)
    return fc_cart


def net_filter(M, modules, net_list):
    modules = np.sort(np.asarray(modules))
    net_filt = [True if elem in net_list else False for elem in modules]
    M1 = M[net_filt, :]
    M2 = M1[:, net_filt]
    return M2

def net_filter_4d(M, modules, net_list):
    modules = np.sort(np.asarray(modules))
    net_filt = [True if elem in net_list else False for elem in modules]
    M1 = M[:, :, net_filt, :]
    M2 = M1[:, :, :, net_filt]
    return M2

def single_window_allegiance(t):
    n_nod = len(t)
    T = np.zeros((n_nod, n_nod))

    for i, row in enumerate(t):
        for j, col in enumerate(t):
            if row == col:
                T[i,j] = 1
            else:
                continue
    return T

def single_window_allegiance_mean(M):
    n_opt = M.shape[0]
    n_nod = M.shape[1]
    T = np.zeros((n_opt, n_nod, n_nod))

    for i in range(n_opt):
        T[i, :, :] = single_window_allegiance(M[i, :])
    
    return T.mean(axis=0)

def all_window_allegiance_mean(M):

    n_win = M.shape[2]
    n_nod = M.shape[1]

    T = np.zeros((n_win, n_nod, n_nod))

    for i in range(n_win):
        T[i, :, :] =  single_window_allegiance_mean(M[:,:,i])
    
    return T

def node_flexibility(node_vector):
    n_win = len(node_vector)
    count = 0
    for i in range(n_win-1):
        if node_vector[i] != node_vector[i+1]:
            count += 1
    return count/n_win


def flexibility(matrix):
    n_nod = matrix.shape[0]
    n_win = matrix.shape[1]

    node_flex = [node_flexibility(elem) for elem in matrix]
    return node_flex



def randmio_und(R, itr):
    '''
    Edited from aestrivex/bctpy with corrected matrix symetry check
    
    This function randomizes an undirected network, while preserving the
    degree distribution. The function does not preserve the strength
    distribution in weighted networks.
    Parameters
    ----------
    W : NxN np.ndarray
        undirected binary/weighted connection matrix
    itr : int
        rewiring parameter. Each edge is rewired approximately itr times.
    Returns
    -------
    R : NxN np.ndarray
        randomized network
    eff : int
        number of actual rewirings carried out
    '''
    if not np.allclose(R, R.T):
        raise BCTParamError("Input must be undirected")
    R = R.copy()
    n = len(R)
    i, j = np.where(np.tril(R))
    k = len(i)
    itr *= k

    # maximum number of rewiring attempts per iteration
    max_attempts = np.round(n * k / (n * (n - 1)))
    # actual number of successful rewirings
    eff = 0

    for it in range(int(itr)):
        att = 0
        while att <= max_attempts:  # while not rewired
            while True:
                e1, e2 = np.random.randint(k, size=(2,))
                while e1 == e2:
                    e2 = np.random.randint(k)
                a = i[e1]
                b = j[e1]
                c = i[e2]
                d = j[e2]

                if a != c and a != d and b != c and b != d:
                    break  # all 4 vertices must be different

            if np.random.random() > .5:
                i.setflags(write=True)
                j.setflags(write=True)
                i[e2] = d
                j[e2] = c  # flip edge c-d with 50% probability
                c = i[e2]
                d = j[e2]  # to explore all potential rewirings

            # rewiring condition
            if not (R[a, d] or R[c, b]):
                R[a, d] = R[a, b]
                R[a, b] = 0
                R[d, a] = R[b, a]
                R[b, a] = 0
                R[c, b] = R[c, d]
                R[c, d] = 0
                R[b, c] = R[d, c]
                R[d, c] = 0

                j.setflags(write=True)
                j[e1] = d
                j[e2] = b  # reassign edge indices
                eff += 1
                break
            att += 1

    return R, eff