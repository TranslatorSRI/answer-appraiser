
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
import numpy as np
from scipy.sparse import csc_matrix

import dcc.file_utils as futils 
import dcc.dcc_utils as dutils 

# constants
logger = dutils.get_logger(__name__)

# session vars
# matrix_gene_set = None
# map_gene_set_indexes_result = None
# if not matrix_gene_set:
#     map_gene_set, map_gene_set_indexes = load_geneset_set_matrix()



# methods
def load_geneset_matrix(map_gene_index, list_gene_set_files, path_gene_set_files, log=False):
    '''  
    will create the sparse matrix used for the pigean calculation
    '''
    # initialize
    list_row = []
    list_columns = []
    list_data = []
    count_column = 0
    data_value = 1
    matrix_result = None 
    map_gene_set_indexes_result = {}

    # read the file
    for file in list_gene_set_files:
        map_gene_set = futils.read_tab_delimited_file(path_gene_set_files + file)

        # log
        logger.info("read gene set file: {} with num entries: {}".format(file, len(map_gene_set)))

        # for each line, add 
        for gene_set_from_file, list_genes_from_file in map_gene_set.items():
            # log
            if log:
                print("\nfor set: {}, got geneset list: {}".format(gene_set_from_file, list_genes_from_file))

            # for each gene, add data
            for gene in list_genes_from_file:
                if map_gene_index.get(gene) is not None:
                    list_row.append(map_gene_index.get(gene))
                    list_columns.append(count_column)
                    list_data.append(data_value)

                    # log
                    if log:
                        print("adding gene: {} for row: {} and column: {}".format(gene, map_gene_index.get(gene), count_column))
                
                else:
                    if log:
                        print("NOT adding gene: {} for row: {} and column: {}".format(gene, map_gene_index.get(gene), count_column))

        
            # add the gene set to the list
            map_gene_set_indexes_result[count_column] = gene_set_from_file

            # increment the column count for each gene set
            count_column = count_column + 1

    # log
    if log:
        logger.info("rows: {}".format(list_row))
        logger.info("cols: {}".format(list_columns))
        logger.info("data: {}".format(list_data))

        print("\nrows: {}".format(list_row))
        print("cols: {}".format(list_columns))
        print("data: {}".format(list_data))

    # create the matrix
    # NOTE - remove toarray due to errors in compute_utils._calc_X_shift_scale() for .A1 and .power()
    # matrix_result = csc_matrix((list_data, (list_row, list_columns)), shape=(len(map_gene_index), len(map_gene_set_indexes_result))).toarray()
    matrix_result = csc_matrix((list_data, (list_row, list_columns)), shape=(len(map_gene_index), len(map_gene_set_indexes_result)))

    # log
    logger.info("returning gene set matrix of shape: {}".format(matrix_result.shape))
    
    # return
    return matrix_result, map_gene_set_indexes_result
    

def generate_gene_vector_from_list(list_gene, map_gene_index, log=False):
    '''
    will generate the gene vector from the list of genes
    '''
    # initialize
    vector_result = None 
    list_row = []
    list_columns = []
    list_data = []
    data_value = 1
    
    # from the genes, create a sparse vector
    for gene in list_gene:
        if map_gene_index.get(gene) is not None:
            # was getting (m, 1) vector, when wanted (1,m) vector
            # list_row.append(map_gene_index.get(gene))
            # list_columns.append(0)
            list_columns.append(map_gene_index.get(gene))
            list_row.append(0)
            list_data.append(data_value)

    # create the vector
    vector_result = csc_matrix((list_data, (list_row, list_columns)), shape=(1, len(map_gene_index))).toarray()

    # log
    if log:
        print("\nrows: {}".format(list_row))
        print("cols: {}".format(list_columns))
        print("data: {}".format(list_data))

    # return 
    return vector_result, list_columns

def sum_of_gene_row(sparse_matrix, gene_index, log=False):
    '''
    Function to compute the sum of all entries for a specific gene row
    '''
    return sparse_matrix[gene_index].sum()

# main
if __name__ == "__main__":
    pass