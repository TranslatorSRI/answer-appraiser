
# LICENSE
# Copyright 2024 Flannick Lab

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

# imports
import scipy
import scipy.sparse as sparse
import scipy.stats
from scipy.stats import truncnorm, norm
from scipy.special import psi as digamma
from scipy.special import erfc
import numpy as np
from numpy.random import gamma
from numpy.random import normal
from numpy.random import exponential
from sklearn.decomposition import NMF
import json

import dcc.dcc_utils as dutils 
import dcc.matrix_utils as mutils 

# TODO - process notes
# 1. done - X matrix is sitting in memory (something like 20K genes x 40K gene sets)
# 2. done - Request comes in with Y vector, need to reorder the Y vector to match the order of your X matrix
# 3. done - p_values = compute_beta_tildes (X, Y)
# 4. done - Filter V = X[:,p_values < 0.05] (edited) 
# 5. done - Then V = V[X.sum(axis=1) > 0,:] (edited) 
# 6. Then factor with V passed in as V0 (edited) (_bayes_nmf_l2())
# Y = np.zeros(X.shape[0])

# W is going to be gene x factor loadings
# K is gene set x factor loadings
# The factors are your groups
# The top scoring gene sets are your labels for the factors
# Then top scoring genes are the ones in the factor
# (Each factor is a cluster)

# constants
logger = dutils.get_logger(__name__)


BATCH_SIZE = 4500

class RunFactorException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

# methods
def calculate_factors(matrix_gene_sets_gene_original, list_gene, list_system_genes, map_gene_index, map_gene_set_index, mean_shifts, scale_factors, p_value=0.05, log=False):
    '''
    will produce the gene set factors and gene factors
    '''
    # initialize
    list_factor = [] 
    list_factor_genes = []
    list_factor_gene_sets = []
    gene_factor = None
    gene_set_factor = None
    map_lowest_factor_per_gene = {}

    # step 1/2: get the gene vector from the gene list
    vector_gene, list_input_gene_indices = mutils.generate_gene_vector_from_list(list_gene=list_gene, map_gene_index=map_gene_index)

    # log
    if log:
        print("step 1: got gene set matrix of shape: {}".format(matrix_gene_sets_gene_original.shape))
        print("step 1: got mean_shifts of shape: {}".format(mean_shifts.shape))
        print("step 1: got scale_factors of shape: {}".format(scale_factors.shape))
        print("step 2: got gene vector of shape: {}".format(vector_gene.shape))

    # step 3: get the p_values by gene set
    vector_gene_set_pvalues = compute_beta_tildes(X=matrix_gene_sets_gene_original, Y=vector_gene, scale_factors=scale_factors, mean_shifts=mean_shifts)

    if log:
        print("step 3: got p values vector of shape: {}".format(vector_gene_set_pvalues.shape))
        print("step 3: filtering gene sets using p_value: {}".format(p_value))

    # step 4: filter the gene set columns based on computed pvalue for each gene set
    matrix_gene_set_filtered_by_pvalues, selected_gene_set_indices = filter_matrix_columns(matrix_input=matrix_gene_sets_gene_original, vector_input=vector_gene_set_pvalues, 
                                                                                           cutoff_input=p_value, log=log)
    # matrix_gene_set_filtered_by_pvalues, selected_gene_set_indices = filter_matrix_columns(matrix_input=matrix_gene_sets_gene_original, vector_input=vector_gene_set_pvalues, 
    #                                                                                        cutoff_input=0.5, log=log)

    if log:
        print("step 4: got gene set filtered (col) matrix of shape: {}".format(matrix_gene_set_filtered_by_pvalues.shape))
        print("step 4: got gene set filtered indices of length: {}".format(len(selected_gene_set_indices)))
        print("step 4: got gene set filtered indices: {}".format(selected_gene_set_indices))

    # step 5: filter gene rows by only the genes that are part of the remaining gene sets from the filtered gene set matrix
    matrix_gene_filtered_by_remaining_gene_sets, selected_gene_indices = filter_matrix_rows_by_sum_cutoff(matrix_to_filter=matrix_gene_set_filtered_by_pvalues, 
                                                                                                          matrix_to_sum=matrix_gene_set_filtered_by_pvalues, log=log)
    # check how many genes are left out
    list_input_genes_filtered_out_indices = [item for item in list_input_gene_indices if item not in selected_gene_indices.tolist()]    

    if log:
        print("step 5: ===> got input gene filtered out of length: {}".format(len(list_input_genes_filtered_out_indices)))
        print("step 5: got gene filtered indices of length: {}".format(len(selected_gene_indices)))
        print("step 5: ===> got gene filtered (rows) matrix of shape: {} to start bayes NMF".format(matrix_gene_filtered_by_remaining_gene_sets.shape))
        # print("step 5: got gene filtered indices of length: {}".format(selected_gene_indices.shape))

    if not all(dim > 0 for dim in matrix_gene_filtered_by_remaining_gene_sets.shape):
        print("step 6: ===> skipping due to pre bayes NMF matrix of shape".format(matrix_gene_filtered_by_remaining_gene_sets.shape))

    else:
        # step 6: from this double filtered matrix, compute the factors
        gene_factor, gene_set_factor, _, _, exp_lambda, _ = _bayes_nmf_l2(V0=matrix_gene_filtered_by_remaining_gene_sets)
        # gene_factor, gene_set_factor = run_nmf(matrix_input=matrix_gene_filtered_by_remaining_gene_sets, log=log)

        if log:
            print("step 6: got gene factor matrix of shape: {}".format(gene_factor.shape))
            print("step 6: got gene set factor matrix of shape: {}".format(gene_set_factor.shape))
            print("step 6: got lambda matrix of shape: {} with data: {}".format(exp_lambda.shape, exp_lambda))

        # step 7: find and rank the gene and gene set groups
        list_factor, list_factor_genes, list_factor_gene_sets, updated_gene_factors = rank_gene_and_gene_sets(X=None, Y=None, exp_lambdak=exp_lambda, exp_gene_factors=gene_factor, exp_gene_set_factors=gene_set_factor.T,
                                                                        list_system_genes=list_system_genes, map_gene_set_index=map_gene_set_index, 
                                                                        list_gene_mask=selected_gene_indices, list_gene_set_mask=selected_gene_set_indices, log=log)

        # step 7a - get the lowest factor per gene
        map_lowest_factor_per_gene = get_lowest_gene_factor_by_gene(exp_gene_factors=updated_gene_factors, list_system_genes=list_system_genes, list_gene_mask=selected_gene_indices, log=False)
        # print(json.dumps(map_lowest_factor_per_gene, indent=2))

        if log:
            print("step 7: got factor list: {}".format(list_factor))
            print("step 7: got gene list:")
            for row in list_factor_genes: 
                print (row)
            print("step 7: got gene set list:")
            for row in list_factor_gene_sets: 
                print (row)


    # only return the gene factors and gene set factors
    return list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, gene_set_factor, map_lowest_factor_per_gene


def group_factor_results(list_factor, list_factor_genes, list_factor_gene_sets, log=False):
    '''
    will group the results that come out of the computation steps
    '''
    list_result = []

    # loop through the factors
    for index, row in enumerate(list_factor):
        list_result.append({'top_set': row, 'gene_sets': list_factor_gene_sets[index], 'genes': list_factor_genes[index]})

    # return
    return list_result


def compute_beta_tildes(X, Y, scale_factors, mean_shifts, y_var=1, resid_correlation_matrix=None, log=False):
    '''
    get the scale factors and mean shifts from _calc_X_shift_scale()
    TODO - check if y_var default 1 is ok (Jason)
    '''

    logger.info("Calculating beta tildes")

    # Y can be a matrix with dimensions:
    # number of parallel runs x number of gene sets
    if len(Y.shape) == 2:
        len_Y = Y.shape[1]
        Y = (Y.T - np.mean(Y, axis=1)).T
    else:
        len_Y = Y.shape[0]
        Y = Y - np.mean(Y)

    dot_product = np.array(X.T.dot(Y.T) / len_Y).T

    variances = np.power(scale_factors, 2)

    # avoid divide by 0 only
    variances[variances == 0] = 1

    #multiply by scale factors because we store beta_tilde in units of scaled X
    beta_tildes = scale_factors * dot_product / variances

    if len(Y.shape) == 2:
        ses = np.outer(np.sqrt(y_var), scale_factors)
    else:
        ses = np.sqrt(y_var) * scale_factors

    ses /= (np.sqrt(variances * (len_Y - 1)))

    # FIXME: implement exact SEs
    # rather than just using y_var as a constant, calculate X.multiply(beta_tildes)
    # then, subtract out Y for non-zero entries, sum square, sum total
    # then, add in square of Y for zero entries, add in total
    # use these to calculate the variance term

    se_inflation_factors = None
    if resid_correlation_matrix is not None:
        logger.info("Adjusting standard errors for correlations")
        # need to multiply by inflation factors: (X * sigma * X) / variances

        # SEs and betas are stored in units of centered and scaled X
        # we do not need to scale X here, however, because cor_variances will then be in units of unscaled X
        # since variances are also in units of unscaled X, these will cancel out

        r_X = resid_correlation_matrix.dot(X)
        r_X_col_means = r_X.multiply(X).sum(axis=0).A1 / X.shape[0]
        cor_variances = r_X_col_means - np.square(r_X_col_means)
        
        # never increase significance
        cor_variances[cor_variances < variances] = variances[cor_variances < variances]

        # both cor_variances and variances are in units of unscaled X
        se_inflation_factors = np.sqrt(cor_variances / variances)

    # MPD - change so only return pvalues
    # return finalize_regression(beta_tildes=beta_tildes, ses=ses, se_inflation_factors=se_inflation_factors)
    (cal_beta_tildes, cal_ses, cal_z_scores, cal_p_values, cal_se_inflation_factors) = finalize_regression(beta_tildes=beta_tildes, ses=ses, se_inflation_factors=se_inflation_factors)
    return cal_p_values

def finalize_regression(beta_tildes, ses, se_inflation_factors):

    if se_inflation_factors is not None:
        ses *= se_inflation_factors

    if np.prod(ses.shape) > 0:
        # empty mask
        empty_mask = np.logical_and(beta_tildes == 0, ses <= 0)
        max_se = np.max(ses)

        ses[empty_mask] = max_se * 100 if max_se > 0 else 100

        # if no y var, set beta tilde to 0

        beta_tildes[ses <= 0] = 0

    z_scores = np.zeros(beta_tildes.shape)
    ses_positive_mask = ses > 0
    z_scores[ses_positive_mask] = beta_tildes[ses_positive_mask] / ses[ses_positive_mask]
    if np.any(~ses_positive_mask):
        logger.info("There were %d gene sets with negative ses; setting z-scores to 0" % (np.sum(~ses_positive_mask)))

    p_values = 2*scipy.stats.norm.cdf(-np.abs(z_scores))

    # pvalues is what I want
    return (beta_tildes, ses, z_scores, p_values, se_inflation_factors)



# NOTE - this code is adapted from https://github.com/gwas-partitioning/bnmf-clustering
# NOTE - paper is: Smith, K., Deutsch, A.J., McGrail, C. et al. Multi-ancestry polygenic mechanisms of type 2 diabetes. Nat Med 30, 1065–1074 (2024). https://doi.org/10.1038/s41591-024-02865-3
# 
# def _bayes_nmf_l2(V0, n_iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0):
def _bayes_nmf_l2(V0, n_iter=10000, a0=10, tol=1e-3, K=15, K0=15, phi=1.0):
    '''
    example?
        result = _bayes_nmf_l2(matrix, a0=alpha0, K=max_num_factors, K0=max_num_factors)
        exp_lambdak = result[4]
        exp_gene_factors = result[1].T
        exp_gene_set_factors = result[0]
    '''
    # seed the random call
    np.random.seed(42)

    # convert the sparse matrix tto a dense array
    V0 = V0.toarray()

    # Bayesian NMF with half-normal priors for W and H
    # V0: input z-score matrix (variants x traits)
    # n_iter: Number of iterations for parameter optimization
    # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
    # tol: Tolerance for convergence of fitting procedure
    # K: Number of clusters to be initialized (algorithm may drive some to zero)
    # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
    # phi: Scaling parameter

    eps = 1.e-50
    delambda = 1.0
    #active_nodes = np.sum(V0, axis=0) != 0
    #V0 = V0[:,active_nodes]
    V = V0 - np.min(V0)
    Vmin = np.min(V)
    Vmax = np.max(V)
    N = V.shape[0]
    M = V.shape[1]

    W = np.random.random((N, K)) * Vmax #NxK
    H = np.random.random((K, M)) * Vmax #KxM

    I = np.ones((N, M)) #NxM
    V_ap = W.dot(H) + eps #NxM

    # print("1 V is: {} of dimension: {}".format(type(V), V.shape))
    temp = np.std(V)
    # print("2 V is: {} of dimension: {}".format(type(V), V.shape))
    temp = np.power(temp, 2) 
    # print("3 V is: {} of dimension: {}".format(type(V), V.shape))
    temp = temp * phi
    # print("4 V is: {} of dimension: {}".format(type(V), V.shape))

    phi = np.power(np.std(V), 2) * phi
    C = (N + M) / 2 + a0 + 1
    b0 = 3.14 * (a0 - 1) * np.mean(V) / (2 * K0)
    lambda_bound = b0 / C
    lambdak = (0.5 * np.sum(np.power(W, 2), axis=0) + 0.5 * np.sum(np.power(H, 2), axis=1) + b0) / C
    lambda_cut = lambda_bound * 1.5

    n_like = [None]
    n_evid = [None]
    n_error = [None]
    n_lambda = [lambdak]
    it = 1
    count = 1
    while delambda >= tol and it < n_iter:
        H = H * (W.T.dot(V)) / (W.T.dot(V_ap) + phi * H * np.repeat(1/lambdak, M).reshape(len(lambdak), M) + eps)
        V_ap = W.dot(H) + eps
        W = W * (V.dot(H.T)) / (V_ap.dot(H.T) + phi * W * np.tile(1/lambdak, N).reshape(N, len(lambdak)) + eps)
        V_ap = W.dot(H) + eps
        lambdak = (0.5 * np.sum(np.power(W, 2), axis=0) + 0.5 * np.sum(np.power(H, 2), axis=1) + b0) / C
        delambda = np.max(np.abs(lambdak - n_lambda[it - 1]) / n_lambda[it - 1])
        like = np.sum(np.power(V - V_ap, 2)) / 2
        n_like.append(like)
        n_evid.append(like + phi * np.sum((0.5 * np.sum(np.power(W, 2), axis=0) + 0.5 * np.sum(np.power(H, 2), axis=1) + b0) / lambdak + C * np.log(lambdak)))
        n_lambda.append(lambdak)
        n_error.append(np.sum(np.power(V - V_ap, 2)))
        if it % 100 == 0:
            logger.info("Iteration=%d; evid=%.3g; lik=%.3g; err=%.3g; delambda=%.3g; factors=%d; factors_non_zero=%d" % (it, n_evid[it], n_like[it], n_error[it], delambda, np.sum(np.sum(W, axis=0) != 0), np.sum(lambdak >= lambda_cut)))
        it += 1

    return W, H, n_like[-1], n_evid[-1], n_lambda[-1], n_error[-1]
        #W # Variant weight matrix (N x K)
        #H # Trait weight matrix (K x M)
        #n_like # List of reconstruction errors (sum of squared errors / 2) pe_bayes_nmf_l2r iteration
        #n_evid # List of negative log-likelihoods per iteration
        #n_lambda # List of lambda vectors (shared weights for each of K clusters, some ~0) per iteration
        #n_error # List of reconstruction errors (sum of squared errors) per iteration

        # The key is that the lambdk (n_lambda) values tell you which factors remain vs. are zeroed out
        # And then genes (H) / gene sets (W) are assigned to the factors for which they have the highest weights
        # There is then logic to “score” each cluster and also annotate them with an LMM, but that is probably not necessary (yet)

def rank_gene_and_gene_sets(X, Y, exp_lambdak, exp_gene_factors, exp_gene_set_factors, list_gene_mask, list_gene_set_mask, list_system_genes, map_gene_set_index, cutoff=1e-5, log=False):
    '''
    will rank the gene sets and gene factors
    '''
    # self.exp_lambdak = result[4]
    # self.exp_gene_factors = result[1].T
    # self.exp_gene_set_factors = result[0]

    # log
    if log:
        print("got lambda of shape: {}".format(exp_lambdak.shape))
        print("got gene factor of shape: {}".format(exp_gene_factors.shape))
        print("got gene set factor of shape: {}".format(exp_gene_set_factors.shape))

    # subset_down
    # GUESS: filter and keep if exp_lambdak > 0 and at least one non zero factor for a gene and gene set; then filter by cutoff
    factor_mask = exp_lambdak != 0 & (np.sum(exp_gene_factors, axis=0) > 0) & (np.sum(exp_gene_set_factors, axis=0) > 0)
    factor_mask = factor_mask & (np.max(exp_gene_set_factors, axis=0) > cutoff * np.max(exp_gene_set_factors))

    if log:
        print("end up with factor mask of shape: {} and true count: {}".format(factor_mask.shape, np.sum(factor_mask)))

    # TODO - QUESTION
    # filter by factors; why invert factor_mask?
    if np.sum(~factor_mask) > 0:
        exp_lambdak = exp_lambdak[factor_mask]
        exp_gene_factors = exp_gene_factors[:,factor_mask]
        exp_gene_set_factors = exp_gene_set_factors[:,factor_mask]
    # elif self.betas_uncorrected is not None:
    #     gene_set_values = self.betas_uncorrected

    if log:
        print("got NEW shrunk lambda of shape: {}".format(exp_lambdak.shape))
        print("got NEW shrunk gene factor of shape: {}".format(exp_gene_factors.shape))
        print("got NEW shrunk gene set factor of shape: {}".format(exp_gene_set_factors.shape))

    # gene_values = None
    # if self.combined_prior_Ys is not None:
    #     gene_values = self.combined_prior_Ys
    # elif self.priors is not None:
    #     gene_values = self.priors
    # elif self.Y is not None:
    #     gene_values = self.Y

    # if gene_set_values is not None:
    #     self.factor_gene_set_scores = self.exp_gene_set_factors.T.dot(gene_set_values[self.gene_set_factor_gene_set_mask])
    # else:
    #     self.factor_gene_set_scores = self.exp_lambdak

    factor_gene_set_scores = exp_lambdak

    # if gene_values is not None:
    #     self.factor_gene_scores = self.exp_gene_factors.T.dot(gene_values[self.gene_factor_gene_mask]) / self.exp_gene_factors.shape[0]
    # else:
    #     self.factor_gene_scores = self.exp_lambdak

    factor_gene_scores = exp_lambdak

    # get the indices with factors sorted in descending order
    reorder_inds = np.argsort(-factor_gene_set_scores)
    exp_lambdak = exp_lambdak[reorder_inds]

    factor_gene_set_scores = factor_gene_set_scores[reorder_inds]
    factor_gene_scores = factor_gene_scores[reorder_inds]
    exp_gene_factors = exp_gene_factors[:,reorder_inds]
    exp_gene_set_factors = exp_gene_set_factors[:,reorder_inds]

    #now label them
    # TODO - QUESTION - do I feed in the filtered ids here?
    # gene_set_factor_gene_set_inds = np.where(self.gene_set_factor_gene_set_mask)[0]
    # gene_factor_gene_inds = np.where(self.gene_factor_gene_mask)[0]

    # TODO - could make top count factors returned a variable; currently constant in code
    num_top = 5

    # get the top count for gene set and genes
    top_gene_inds = np.argsort(-exp_gene_factors, axis=0)[:num_top,:]               # not used
    top_gene_set_inds = np.argsort(-exp_gene_set_factors, axis=0)[:num_top,:]

    factor_labels = []
    top_gene_sets = []
    top_genes = []
    factor_prompts = []

    # log
    if log:
        print("looping through factor gene set scores of size: {} and data: \n{}".format(len(factor_gene_set_scores), factor_gene_set_scores))
        print("got top pathway ids type: {} and data: {}".format(type(top_gene_set_inds), top_gene_set_inds))
        print("got top gene ids: {}".format(top_gene_inds))

    for i in range(len(factor_gene_set_scores)):
        # orginal for reference
        # top_gene_sets.append([self.gene_sets[i] for i in np.where(self.gene_set_factor_gene_set_mask)[0][top_gene_set_inds[:,i]]])
        # top_genes.append([self.genes[i] for i in np.where(self.gene_factor_gene_mask)[0][top_gene_inds[:,i]]])
        
        # list_temp = top_gene_inds[:,i].tolist()
        # print("top gene ids: {}".format(list_temp))

        # build the list of genes and gene sets that were filtered first by p_value, then the factor process in this method
        list_gene_index_factor = get_referenced_list_elements(list_referenced=list_gene_mask, list_index=top_gene_inds[:,i].tolist(), log=False)
        list_gene_set_index_factor = get_referenced_list_elements(list_referenced=list_gene_set_mask, list_index=top_gene_set_inds[:,i].tolist(), log=False)

        # print("map type: {}".format(type(map_gene_set_index)))
        # print("list type: {}".format(type(list_gene_set_index_factor)))
        # print("list: {}".format(list_gene_set_index_factor))

        # if log:
        #     print("got pathway indexes: {}".format(list_gene_set_index_factor))
        #     print("got gene indexes: {}".format(list_gene_index_factor))

        # build the list of groupings of gene sets
        top_gene_sets.append([map_gene_set_index.get(gs_index) for gs_index in list_gene_set_index_factor])

        # build the list of groupings of genes (cutoff of 0.01)
        # top_genes.append([list_system_genes[g_index] for g_index in list_gene_index_factor])
        # top_genes.append([list_system_genes[g_index] for index_local, g_index in enumerate(list_gene_index_factor) if exp_gene_factors[top_gene_inds[index_local], i] > 0.01])
 
        list_temp = []
        for index_local, g_index in enumerate(list_gene_index_factor):
            score_gene = exp_gene_factors[top_gene_inds[index_local, i], i]
            if score_gene > 0.01:
                list_temp.append({'gene': list_system_genes[g_index], 'score': score_gene})
        top_genes.append(list_temp)

        factor_labels.append(top_gene_sets[i][0] if len(top_gene_sets[i]) > 0 else "")
        factor_prompts.append(",".join(top_gene_sets[i]))

    # return the 3 grouping data structures, then the updated (filtered) gene factors
    return factor_labels, top_genes, top_gene_sets, exp_gene_factors


def get_lowest_gene_factor_by_gene(exp_gene_factors, list_system_genes, list_gene_mask, log=False):
    '''
    will return the lowest factor per gene - to be used as novelty calculation for ARS
    '''
    # initialize
    map_result = {}

    if all(dim > 0 for dim in exp_gene_factors.shape):    
        # log
        if log:
            print("lowest factor - got gene factor of shape: {}".format(exp_gene_factors.shape))
            # print("lowest factor - got filtered gene mask of size: {} and data: \n{}".format(len(list_gene_mask), list_gene_mask))

        # get the lowest value per row
        min_per_row = np.min(exp_gene_factors, axis=1)

        if log:
            print("lowest factor - got gene factor MINIMUM of shape: {} and type: {}".format(min_per_row.shape, type(min_per_row)))
            for index in range(len(list_gene_mask)):
                print("lowest factor - for gene: {} get factor : {}".format(list_system_genes[list_gene_mask[index]], exp_gene_factors[index]))

        # build the map
        if min_per_row is not None:
            for index, row_factor in enumerate(min_per_row.tolist()):
                index_system = list_gene_mask[index]
                gene_name = list_system_genes[index_system]
                map_result[gene_name] = row_factor

        # if log:
        #     print("lowest factor - got gene factor MINIMUM map of size: {} and data: {}".format(len(map_result), map_result))

        logger.info("returning lowest factor (novelty) gene map is size: {}".format(len(map_result)))

    else:
        logger.info("returning empty map for exp_gene_factor of shape: {}".format(exp_gene_factors.shape))


    # return
    return map_result


def get_referenced_list_elements(list_referenced, list_index, log=False):
    '''
    will return a list of the referenced elements from the list
    '''
    list_result = []

    # log
    if log:
        print("ref list: {}".format(list_referenced))
        print("index list: {}".format(list_index))

    # get the elements
    list_result = [list_referenced[i] for i in list_index]

    # return
    return list_result



def _calc_X_shift_scale(X, y_corr_cholesky=None):
    '''
    returns the mean shifts and scale factors for the initial pathwayxgene matrix X
    '''
    if y_corr_cholesky is None:
        mean_shifts = X.sum(axis=0).A1 / X.shape[0]
        scale_factors = np.sqrt(X.power(2).sum(axis=0).A1 / X.shape[0] - np.square(mean_shifts))

        # NOTE - delete? mods not necessary if don't convert cscsparse.toarray()
        # mean_shifts = np.ravel(X.sum(axis=0)) / X.shape[0]
        # scale_factors = np.sqrt(X.power(2).sum(axis=0).A1 / X.shape[0] - np.square(mean_shifts))
    else:
        scale_factors = np.array([])
        mean_shifts = np.array([])
        for X_b, begin, end, batch in _get_X_blocks_internal(X, y_corr_cholesky):
            (cur_mean_shifts, cur_scale_factors) = _calc_shift_scale(X_b)
            mean_shifts = np.append(mean_shifts, cur_mean_shifts)
            scale_factors = np.append(scale_factors, cur_scale_factors)

    # log
    logger.info("returning mean shifts matrix of shape: {}".format(mean_shifts.shape))
    logger.info("returning scale factors matrix of shape: {}".format(scale_factors.shape))

    return (mean_shifts, scale_factors)


def _calc_shift_scale(X_b):
    '''
    used in _calc_X_shift_scale()
    '''
    mean_shifts = []
    scale_factors = []
    for i in range(X_b.shape[1]):
        X_i = X_b[:,i]
        mean_shifts.append(np.mean(X_i))
        scale_factor = np.std(X_i)
        if scale_factor == 0:
            scale_factor = 1
        scale_factors.append(scale_factor)
    return (np.array(mean_shifts), np.array(scale_factors))


def _get_X_blocks_internal(X_orig, y_corr_cholesky, whiten=True, full_whiten=False, start_batch=0, mean_shifts=None, scale_factors=None):
    ''' 
    used in _calc_X_shift_scale()
    '''

    if y_corr_cholesky is None:
        #explicitly turn these off to help with caching
        whiten = False
        full_whiten = False

    num_batches = _get_num_X_blocks(X_orig)

    consider_cache = X_orig is X_orig and num_batches == 1 and mean_shifts is None and scale_factors is None

    for batch in range(start_batch, num_batches):
        # logger.info("Getting X%s block batch %s (%s)" % ("_missing" if X_orig is X_orig_missing_gene_sets else "", batch, "fully whitened" if full_whiten else ("whitened" if whiten else "original")), TRACE)
        begin = batch * BATCH_SIZE
        end = (batch + 1) * BATCH_SIZE
        if end > X_orig.shape[1]:
            end = X_orig.shape[1]

        if last_X_block is not None and consider_cache and last_X_block[1:] == (whiten, full_whiten, begin, end, batch):
            logger.info("Using cache!")
            yield (last_X_block[0], begin, end, batch)
        else:
            X_b = X_orig[:,begin:end].toarray()
            if mean_shifts is not None:
                X_b = X_b - mean_shifts[begin:end]
            if scale_factors is not None:
                X_b = X_b / scale_factors[begin:end]

            if whiten or full_whiten:
                X_b = _whiten(X_b, y_corr_cholesky, whiten=whiten, full_whiten=full_whiten)

            #only cache if we are accessing the original X
            if consider_cache:
                last_X_block = (X_b, whiten, full_whiten, begin, end, batch)
            else:
                last_X_block = None

            yield (X_b, begin, end, batch)


def _whiten(matrix, corr_cholesky, whiten=True, full_whiten=False):
    '''
    used in _get_X_blocks_internal()
    '''
    if full_whiten:
        #fully whiten, by sigma^{-1}; useful for optimization
        matrix = scipy.linalg.cho_solve_banded((corr_cholesky, True), matrix, overwrite_b=True)
    elif whiten:
        #whiten X_b by sigma^{-1/2}
        matrix = scipy.linalg.solve_banded((corr_cholesky.shape[0]-1, 0), corr_cholesky, matrix, overwrite_ab=True)
    return matrix


def _get_num_X_blocks(X_orig, batch_size=None):
    '''
    used in _get_X_blocks_internal()
    '''
    if batch_size is None:
        batch_size = BATCH_SIZE
    return int(np.ceil(X_orig.shape[1] / batch_size))


def filter_matrix_columns(matrix_input, vector_input, cutoff_input=0.05, log=False):
    '''
    will filter the matrix based on the vector and cutoff
    '''

    # REFERENCE
    # keep_metrics = np.where(np.any(metric_p_values_m < 0.05, axis=0))[0]

    # this works
    # selected_column_indices = np.where(np.any(vector_input < cutoff_input, axis=0))[0]
    # TODO this does not
    # selected_column_indices = np.where(mask)[0]

    # log
    if log:
        print("got matrix to filter of shape: {} and type: {}".format(matrix_input.shape, type(matrix_input)))
        print("got filter vector of shape: {} and type: {}".format(vector_input.shape, type(vector_input)))

    selected_column_indices = np.where(np.any(vector_input < cutoff_input, axis=0))[0]
    matrix_result = matrix_input[:, selected_column_indices]

    # log
    if log:
        print("got filtered column list of length: {}".format(len(selected_column_indices)))
        print("got resulting shape from column filters from: {} to {}".format(matrix_input.shape, matrix_result.shape))
        # print("example filtered: {}".format(matrix_result[11205]))

    # return
    return matrix_result, selected_column_indices


def filter_matrix_rows_by_sum_cutoff(matrix_to_filter, matrix_to_sum, cutoff_input=0, log=False):
    '''
    will filter the matrix based on sum of each row and cutoff
    '''
    # mask = matrix_to_sum.sum(axis=1) > cutoff_input
    # selected_indices = np.where(mask)[0]
    # keep_metrics = np.where(np.any(mask, axis=0))[0]
    # matrix_result = matrix_to_filter[keep_metrics, :]
    # # matrix_result = matrix_to_filter[mask, :]

    if log:
        print("got matrix to filter of shape: {} and type: {}".format(matrix_to_filter.shape, type(matrix_to_filter)))
        print("got matrix to sum of shape: {} and type: {}".format(matrix_to_sum.shape, type(matrix_to_sum)))

    mask = matrix_to_sum.sum(axis=1) > cutoff_input
    # selected_indices = np.where(mask)[0]
    selected_indices = np.where(np.any(mask, axis=1))[0]
    matrix_result = matrix_to_filter[selected_indices, :]
    # matrix_result = matrix_to_filter[mask, :]

    # log
    if log:
        print("got resulting shape from row sum filters from: {} to {}".format(matrix_to_filter.shape, matrix_result.shape))
        # print("got filter rows indices: {}".format(selected_indices))
        # print("example matrix to sum: {}".format(matrix_to_sum.toarray()[2]))

    # return
    return matrix_result, selected_indices


def run_nmf(matrix_input, num_components=15, log=False):
    '''
    run the sklearn NMF
    '''
    # Initialize the NMF model
    model = NMF(n_components=num_components, random_state=42)

    # Fit the model and transform the matrix X
    W = model.fit_transform(matrix_input)
    H = model.components_    

    # log
    if log:
        print("for gene factor of shape: {}".format(W.shape))
        print("for gene set factor of shape: {}".format(H.shape))

    # return
    return W, H


# main
if __name__ == "__main__":
    pass








# old

# # The sparse matrix X needs to be genes (rows) x gene sets (columns)
# # There is code read_X that does this, but it is probably overkill
# # 
# # You’ll want to read into the sparse matrix at boot up of the server, and then store some auxiliary structures that map gene names to rows and gene set names to columns
# # (b) The user request will then be a list of genes and values
# # You’ll load that into a y vector of length equal to X.shape[0] with genes in the same order
# # (c) You can then run compute_beta_tildes which takes the X and y vectors as inputs (plus some quantities that can be easily computed from X)
# # This will return some values that allow you to filter the columns of the X matrix (I can show you how tomorrow)

# def run_factor(X_orig, max_num_factors=15, alpha0=10, beta0=1, gene_set_filter_type=None, gene_set_filter_value=None, gene_filter_type=None, gene_filter_value=None, gene_set_multiply_type=None, gene_multiply_type=None, run_transpose=True, max_num_iterations=100, rel_tol=0.01, lmm_auth_key=None):

#     BETA_UNCORRECTED_OPTIONS = ["beta_uncorrected", "betas_uncorrected"]
#     BETA_OPTIONS = ["beta", "betas"]
#     P_OPTIONS = ["p_value", "p", "p_value"]
#     NONE_OPTIONS = ["none"]

#     # self._record_params({"max_num_factors": max_num_factors, "alpha0": alpha0, "gene_set_filter_type": gene_set_filter_type, "gene_set_filter_value": gene_set_filter_value, "gene_filter_type": gene_filter_type, "gene_filter_value": gene_filter_value, "gene_set_multiply_type": gene_set_multiply_type, "gene_multiply_type": gene_multiply_type, "run_transpose": run_transpose})

#     gene_set_mask = np.full(self.X_orig.shape[1], True)
#     if gene_set_filter_type is not None and gene_set_filter_value is not None:
#         if gene_set_filter_type.lower() in BETA_UNCORRECTED_OPTIONS:
#             if self.betas_uncorrected is None:
#                 raise RunFactorException("Can't run filtering in factor; betas_uncorrected was not loaded")
#             gene_set_mask = self.betas_uncorrected > gene_set_filter_value
#         elif gene_set_filter_type.lower() in BETA_OPTIONS:
#             if self.betas is None:
#                 raise RunFactorException("Can't run filtering in factor; betas was not loaded")
#             gene_set_mask = self.betas > gene_set_filter_value
#         elif gene_set_filter_type.lower() in P_OPTIONS:
#             if self.p_values is None:
#                 raise RunFactorException("Can't run filtering in factor; p was not loaded")
#             gene_set_mask = self.p_values < gene_set_filter_value
#         elif gene_set_filter_type.lower() in NONE_OPTIONS:
#             pass
#         else:
#             raise RunFactorException("Valid values for --gene-set-filter-type are beta_uncorrected|beta|p_value")

#     gene_set_vector = np.ones(self.X_orig.shape[1])
#     if gene_set_multiply_type is not None:
#         if gene_set_multiply_type.lower() in BETA_UNCORRECTED_OPTIONS:
#             if self.betas_uncorrected is None:
#                 raise RunFactorException("Can't run multiply in factor; betas_uncorrected was not loaded")
#             gene_set_vector = self.betas_uncorrected
#         elif gene_set_multiply_type.lower() in BETA_OPTIONS:
#             if self.betas is None:
#                 raise RunFactorException("Can't run multiply in factor; betas was not loaded")
#             gene_set_vector = self.betas
#         elif gene_set_multiply_type.lower() in P_OPTIONS:
#             if self.p_values is None:
#                 raise RunFactorException("Can't run multiply in factor; p was not loaded")
#             gene_set_vector = -np.log(self.p_values)
#         elif gene_set_multiply_type.lower() in NONE_OPTIONS:
#             pass
#         else:
#             bail("Valid values for --gene-set-multiply-type are beta_uncorrected|beta|p_value")


#     PRIOR_OPTIONS = ["prior", "priors"]
#     COMBINED_OPTIONS = ["combined", "combined_prior_Y", "combined_prior_Ys"]
#     Y_OPTIONS = ["y", "log_bf"]

#     gene_mask = np.full(self.X_orig.shape[0], True)
#     if gene_filter_type is not None and gene_filter_value is not None:
#         if gene_filter_type.lower() in PRIOR_OPTIONS:
#             if self.priors is None:
#                 raise RunFactorException("Can't run filtering in factor; priors were not loaded")
#             gene_mask = self.priors > gene_filter_value
#         elif gene_filter_type.lower() in COMBINED_OPTIONS:
#             if self.combined_prior_Ys is None:
#                 raise RunFactorException("Can't run filtering in factor; combined was not loaded")
#             gene_mask = self.combined_prior_Ys > gene_filter_value
#         elif gene_filter_type.lower() in Y_OPTIONS:
#             if self.Y is None:
#                 raise RunFactorException("Can't run filtering in factor; log_bf was not loaded")
#             gene_mask = self.Y > gene_filter_value
#         elif gene_filter_type.lower() in NONE_OPTIONS:
#             pass
#         else:
#             raise RunFactorException("Valid values for --gene-filter-type are prior|combined|log_bf")

#     gene_vector = np.ones(self.X_orig.shape[0])
#     if gene_multiply_type is not None:
#         if gene_multiply_type.lower() in PRIOR_OPTIONS:
#             if self.priors is None:
#                 raise RunFactorException("Can't run multiply in factor; priors were not loaded")
#             gene_vector = self.priors
#         elif gene_multiply_type.lower() in COMBINED_OPTIONS:
#             if self.combined_prior_Ys is None:
#                 raise RunFactorException("Can't run multiply in factor; combined was not loaded")
#             gene_vector = self.combined_prior_Ys
#         elif gene_multiply_type.lower() in Y_OPTIONS:
#             if self.Y is None:
#                 raise RunFactorException("Can't run multiply in factor; log_bf was not loaded")
#             gene_vector = self.Y
#         elif gene_multiply_type.lower() in NONE_OPTIONS:
#             pass
#         else:
#             raise RunFactorException("Valid values for --gene-multiply-type are prior|combined|log_bf")

#     #make sure everything is positive
#     gene_set_vector[gene_set_vector < 0] = 0
#     gene_vector[gene_vector < 0] = 0

#     logger.info("Running matrix factorization")

#     matrix = self.X_orig[:,gene_set_mask][gene_mask,:].multiply(gene_set_vector[gene_set_mask]).toarray().T * gene_vector[gene_mask]
#     logger.info("Matrix shape: (%s, %s)" % (matrix.shape))

#     if not run_transpose:
#         matrix = matrix.T


#     #this code is adapted from https://github.com/gwas-partitioning/bnmf-clustering
#     def _bayes_nmf_l2(V0, n_iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0):

#         # Bayesian NMF with half-normal priors for W and H
#         # V0: input z-score matrix (variants x traits)
#         # n_iter: Number of iterations for parameter optimization
#         # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
#         # tol: Tolerance for convergence of fitting procedure
#         # K: Number of clusters to be initialized (algorithm may drive some to zero)
#         # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
#         # phi: Scaling parameter

#         eps = 1.e-50
#         delambda = 1.0
#         #active_nodes = np.sum(V0, axis=0) != 0
#         #V0 = V0[:,active_nodes]
#         V = V0 - np.min(V0)
#         Vmin = np.min(V)
#         Vmax = np.max(V)
#         N = V.shape[0]
#         M = V.shape[1]

#         W = np.random.random((N, K)) * Vmax #NxK
#         H = np.random.random((K, M)) * Vmax #KxM

#         I = np.ones((N, M)) #NxM
#         V_ap = W.dot(H) + eps #NxM

#         phi = np.power(np.std(V), 2) * phi
#         C = (N + M) / 2 + a0 + 1
#         b0 = 3.14 * (a0 - 1) * np.mean(V) / (2 * K0)
#         lambda_bound = b0 / C
#         lambdak = (0.5 * np.sum(np.power(W, 2), axis=0) + 0.5 * np.sum(np.power(H, 2), axis=1) + b0) / C
#         lambda_cut = lambda_bound * 1.5

#         n_like = [None]
#         n_evid = [None]
#         n_error = [None]
#         n_lambda = [lambdak]
#         it = 1
#         count = 1
#         while delambda >= tol and it < n_iter:
#             H = H * (W.T.dot(V)) / (W.T.dot(V_ap) + phi * H * np.repeat(1/lambdak, M).reshape(len(lambdak), M) + eps)
#             V_ap = W.dot(H) + eps
#             W = W * (V.dot(H.T)) / (V_ap.dot(H.T) + phi * W * np.tile(1/lambdak, N).reshape(N, len(lambdak)) + eps)
#             V_ap = W.dot(H) + eps
#             lambdak = (0.5 * np.sum(np.power(W, 2), axis=0) + 0.5 * np.sum(np.power(H, 2), axis=1) + b0) / C
#             delambda = np.max(np.abs(lambdak - n_lambda[it - 1]) / n_lambda[it - 1])
#             like = np.sum(np.power(V - V_ap, 2)) / 2
#             n_like.append(like)
#             n_evid.append(like + phi * np.sum((0.5 * np.sum(np.power(W, 2), axis=0) + 0.5 * np.sum(np.power(H, 2), axis=1) + b0) / lambdak + C * np.log(lambdak)))
#             n_lambda.append(lambdak)
#             n_error.append(np.sum(np.power(V - V_ap, 2)))
#             if it % 100 == 0:
#                 log("Iteration=%d; evid=%.3g; lik=%.3g; err=%.3g; delambda=%.3g; factors=%d; factors_non_zero=%d" % (it, n_evid[it], n_like[it], n_error[it], delambda, np.sum(np.sum(W, axis=0) != 0), np.sum(lambdak >= lambda_cut)), TRACE)
#             it += 1

#         return W, H, n_like[-1], n_evid[-1], n_lambda[-1], n_error[-1]
#             #W # Variant weight matrix (N x K)
#             #H # Trait weight matrix (K x M)
#             #n_like # List of reconstruction errors (sum of squared errors / 2) per iteration
#             #n_evid # List of negative log-likelihoods per iteration
#             #n_lambda # List of lambda vectors (shared weights for each of K clusters, some ~0) per iteration
#             #n_error # List of reconstruction errors (sum of squared errors) per iteration


#     result = _bayes_nmf_l2(matrix, a0=alpha0, K=max_num_factors, K0=max_num_factors)
    
#     self.exp_lambdak = result[4]
#     self.exp_gene_factors = result[1].T
#     self.exp_gene_set_factors = result[0]


#     #subset_down
#     factor_mask = self.exp_lambdak != 0 & (np.sum(self.exp_gene_factors, axis=0) > 0) & (np.sum(self.exp_gene_set_factors, axis=0) > 0)
#     factor_mask = factor_mask & (np.max(self.exp_gene_set_factors, axis=0) > 1e-5 * np.max(self.exp_gene_set_factors))

#     if np.sum(~factor_mask) > 0:
#         self.exp_lambdak = self.exp_lambdak[factor_mask]
#         self.exp_gene_factors = self.exp_gene_factors[:,factor_mask]
#         self.exp_gene_set_factors = self.exp_gene_set_factors[:,factor_mask]

#     self.gene_factor_gene_mask = gene_mask
#     self.gene_set_factor_gene_set_mask = gene_set_mask

#     gene_set_values = None
#     if self.betas is not None:
#         gene_set_values = self.betas
#     elif self.betas_uncorrected is not None:
#         gene_set_values = self.betas_uncorrected

#     gene_values = None
#     if self.combined_prior_Ys is not None:
#         gene_values = self.combined_prior_Ys
#     elif self.priors is not None:
#         gene_values = self.priors
#     elif self.Y is not None:
#         gene_values = self.Y

#     if gene_set_values is not None:
#         self.factor_gene_set_scores = self.exp_gene_set_factors.T.dot(gene_set_values[self.gene_set_factor_gene_set_mask])
#     else:
#         self.factor_gene_set_scores = self.exp_lambdak

#     if gene_values is not None:
#         self.factor_gene_scores = self.exp_gene_factors.T.dot(gene_values[self.gene_factor_gene_mask]) / self.exp_gene_factors.shape[0]
#     else:
#         self.factor_gene_scores = self.exp_lambdak


#     reorder_inds = np.argsort(-self.factor_gene_set_scores)
#     self.exp_lambdak = self.exp_lambdak[reorder_inds]
#     self.factor_gene_set_scores = self.factor_gene_set_scores[reorder_inds]
#     self.factor_gene_scores = self.factor_gene_scores[reorder_inds]
#     self.exp_gene_factors = self.exp_gene_factors[:,reorder_inds]
#     self.exp_gene_set_factors = self.exp_gene_set_factors[:,reorder_inds]

#     #now label them
#     gene_set_factor_gene_set_inds = np.where(self.gene_set_factor_gene_set_mask)[0]
#     gene_factor_gene_inds = np.where(self.gene_factor_gene_mask)[0]

#     num_top = 5
#     top_gene_inds = np.argsort(-self.exp_gene_factors, axis=0)[:num_top,:]
#     top_gene_set_inds = np.argsort(-self.exp_gene_set_factors, axis=0)[:num_top,:]

#     self.factor_labels = []
#     self.top_gene_sets = []
#     self.top_genes = []
#     factor_prompts = []
#     for i in range(len(self.factor_gene_set_scores)):
#         self.top_gene_sets.append([self.gene_sets[i] for i in np.where(self.gene_set_factor_gene_set_mask)[0][top_gene_set_inds[:,i]]])
#         self.top_genes.append([self.genes[i] for i in np.where(self.gene_factor_gene_mask)[0][top_gene_inds[:,i]]])
#         self.factor_labels.append(self.top_gene_sets[i][0] if len(self.top_gene_sets[i]) > 0 else "")
#         factor_prompts.append(",".join(self.top_gene_sets[i]))

#     if lmm_auth_key is not None:
#         prompt = "Print a label to assign to each group: %s" % (" ".join(["%d. %s" % (j+1, ",".join(self.top_gene_sets[j])) for j in range(len(self.top_gene_sets))]))
#         log("Querying LMM with prompt: %s" % prompt, DEBUG)
#         response = query_lmm(prompt, lmm_auth_key)
#         if response is not None:
#             try:
#                 responses = response.strip().split("\n")
#                 responses = [x for x in responses if len(x) > 0]
#                 if len(responses) == len(self.factor_labels):
#                     for i in range(len(self.factor_labels)):
#                         self.factor_labels[i] = responses[i]
#                 else:
#                     raise Exception
#             except Exception:
#                 log("Couldn't decode LMM response %s; using simple label" % response)
#                 pass


#     logger.info("Found %d factors" % len(self.exp_lambdak))
