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


# NOTE - this file is to be used to interface with the gene nmf code outside of the Flask app

# imports
import os

from .dcc import startup_utils as sutils
from .dcc import matrix_utils as mutils
from .dcc import compute_utils as cutils
from .dcc import dcc_utils as dutils
from .dcc import sql_utils as sql_utils
from .dcc import gui_utils as gutils


# constants
P_VALUE_CUTOFF = 0.3
MAX_NUMBER_GENE_SETS_FOR_COMPUTATION = 100
current_dir = os.path.dirname(os.path.abspath(__file__))
# DIR_CONF = "conf/"
DIR_CONF = os.path.join(current_dir, "conf/")

# logger
logger = dutils.get_logger(__name__)
logger.info("using configuration directory: {}".format(DIR_CONF))
logger.info("using p_value curtoff of: {}".format(P_VALUE_CUTOFF))

# in memory compute variables
map_conf = sutils.load_conf(dir=DIR_CONF)

# load the database
db_file = os.path.join(current_dir, map_conf.get("root_dir") + map_conf.get("db_file"))
logger.info("loading database file: {}".format(db_file))
sql_connection = sql_utils.db_sqlite_get_connection(db_path=db_file)

# map_gene_index, list_system_genes = futils.load_gene_file_into_map(file_path=map_conf.get('root_dir') + map_conf.get('gene_file'))
map_gene_index, list_system_genes, map_gene_ontology = (
    sql_utils.db_load_gene_table_into_map(conn=sql_connection)
)
matrix_gene_sets, map_gene_set_index = mutils.load_geneset_matrix(
    map_gene_index=map_gene_index,
    list_gene_set_files=map_conf.get("gene_set_files"),
    path_gene_set_files=os.path.join(current_dir, map_conf.get("root_dir")),
    log=False,
)

# get the other 2 cached matrices
(mean_shifts, scale_factors) = cutils._calc_X_shift_scale(X=matrix_gene_sets)

print("================ Bayes NMF data structures LOADED! ===========================")


# methods
def get_gene_nmf_novelty_for_gene_list(
    list_input_genes,
    p_value_cutoff=P_VALUE_CUTOFF,
    max_num_gene_sets=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION,
    log=False,
):
    """
    'will process the gene nmf call for the gene list given and return the gene novelty
    """
    map_result = {}

    # get the calculated data
    map_gene_novelty, list_input_translated = process_genes_novelty(
        list_input_genes=list_input_genes,
        p_value_cutoff=p_value_cutoff,
        max_num_gene_sets=max_num_gene_sets,
    )

    # log result
    logger.info("got novelty result map of size: {}".format(len(map_gene_novelty)))

    # format the data
    map_result = gutils.gui_build_novelty_results_map(
        map_gene_ontology=map_gene_ontology,
        list_input_gene_names=list_input_translated,
        map_gene_index=map_gene_index,
        matrix_gene_sets=matrix_gene_sets,
        map_gene_novelty=map_gene_novelty,
    )

    # return
    return map_result


def get_gene_full_nmf_for_gene_list(
    list_input_genes,
    p_value_cutoff=P_VALUE_CUTOFF,
    max_num_gene_sets=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION,
    log=False,
):
    """
    'will process the gene nmf call for the gene list given and return the full results
    """
    map_result = {}

    # get the calculated data
    (
        list_factor,
        list_factor_genes,
        list_factor_gene_sets,
        map_gene_novelty,
        list_input_translated,
    ) = process_genes_full(
        list_input_genes=list_input_genes,
        p_value_cutoff=p_value_cutoff,
        max_num_gene_sets=max_num_gene_sets,
    )

    # log result
    logger.info("got novelty result map of size: {}".format(len(map_gene_novelty)))

    # format the data
    # map_result = gutils.gui_build_novelty_results_map(map_gene_ontology=map_gene_ontology, list_input_gene_names=list_input_translated, map_gene_index=map_gene_index,
    #                                           matrix_gene_sets=matrix_gene_sets, map_gene_novelty=map_gene_novelty)
    map_result = gutils.gui_build_results_map(
        list_factor=list_factor,
        list_factor_gene_sets=list_factor_gene_sets,
        list_factor_genes=list_factor_genes,
        map_gene_ontology=map_gene_ontology,
        list_input_gene_names=list_input_translated,
        map_gene_index=map_gene_index,
        matrix_gene_sets=matrix_gene_sets,
        map_gene_novelty=map_gene_novelty,
    )

    # return
    return map_result


def process_genes_novelty(
    list_input_genes, p_value_cutoff, max_num_gene_sets, log=False
):
    """
    processes the input genes
    """
    # initialize
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)

    # preprocess
    # translate the genes into what the system can handle
    logger.info("got raw gene inputs of size: {}".format(len(list_input_genes)))
    list_input_translated = sql_utils.db_get_gene_names_from_list(
        conn=sql_conn_query, list_input=list_input_genes
    )
    logger.info(
        "got translated gene inputs of size: {}".format(len(list_input_translated))
    )

    # do the calculations
    (
        list_factor,
        list_factor_genes,
        list_factor_gene_sets,
        gene_factor,
        gene_set_factor,
        map_gene_novelty,
        logs_process,
    ) = cutils.calculate_factors(
        matrix_gene_sets_gene_original=matrix_gene_sets,
        p_value=p_value_cutoff,
        max_num_gene_sets=max_num_gene_sets,
        list_gene=list_input_translated,
        list_system_genes=list_system_genes,
        map_gene_index=map_gene_index,
        map_gene_set_index=map_gene_set_index,
        mean_shifts=mean_shifts,
        scale_factors=scale_factors,
        log=True,
    )

    # log
    for row in logs_process:
        logger.info(row)

    # return
    return map_gene_novelty, list_input_translated


def process_genes_full(list_input_genes, p_value_cutoff, max_num_gene_sets, log=False):
    """
    processes the input genes
    """
    # initialize
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)

    # preprocess
    # translate the genes into what the system can handle
    logger.info("got raw gene inputs of size: {}".format(len(list_input_genes)))
    list_input_translated = sql_utils.db_get_gene_names_from_list(
        conn=sql_conn_query, list_input=list_input_genes
    )
    logger.info(
        "got translated gene inputs of size: {}".format(len(list_input_translated))
    )

    # do the calculations
    (
        list_factor,
        list_factor_genes,
        list_factor_gene_sets,
        gene_factor,
        gene_set_factor,
        map_gene_novelty,
        logs_process,
    ) = cutils.calculate_factors(
        matrix_gene_sets_gene_original=matrix_gene_sets,
        p_value=p_value_cutoff,
        max_num_gene_sets=max_num_gene_sets,
        list_gene=list_input_translated,
        list_system_genes=list_system_genes,
        map_gene_index=map_gene_index,
        map_gene_set_index=map_gene_set_index,
        mean_shifts=mean_shifts,
        scale_factors=scale_factors,
        log=True,
    )

    # log
    for row in logs_process:
        logger.info(row)

    # return
    return (
        list_factor,
        list_factor_genes,
        list_factor_gene_sets,
        map_gene_novelty,
        list_input_translated,
    )


# main
if __name__ == "__main__":
    pass
