"""Test Clinical Evidence function."""

# imports
import pytest 
import logging
import sys
import os

# add path of file
import sys
import os

# Determine the path to startup_utils.py
startup_utils_path = os.path.abspath("app/novelty/gene_nmf/gene_nmf_adapter.py")
startup_utils_dir = os.path.dirname(startup_utils_path)
sys.path.insert(0, startup_utils_dir)
print("added to pyton path the following directory: {}".format(startup_utils_dir))

# import file top test
import app.novelty.gene_nmf.gene_nmf_adapter as adapter
import app.novelty.gene_nmf.dcc.dcc_utils as dutils

# logger
logger = dutils.get_logger(name=__name__)
logger.info("added to pyton path the following directory: {}".format(startup_utils_dir))

# vars
list_gene_test = [
        "NCBIGene:5142",
        "NCBIGene:4853",
        "NCBIGene:4000",
        "NCBIGene:5664",
        "NCBIGene:9804",
        "NCBIGene:22891",
        "NCBIGene:5728",
        "NCBIGene:3087",
        "NCBIGene:5950",
        "NCBIGene:6319",
        "NCBIGene:114815",
        "NCBIGene:6934",
        "NCBIGene:3630",
        "NCBIGene:3784",
        "NCBIGene:3767",
        "NCBIGene:3992",
        "NCBIGene:9415",
        "NCBIGene:7351",
        "NCBIGene:7352",
        "NCBIGene:4544",
        "NCBIGene:10599",
        "NCBIGene:6660",
        "NCBIGene:7480",
        "NCBIGene:8091",
        "NCBIGene:7103",
        "NCBIGene:6927",
        "NCBIGene:2308",
        "NCBIGene:5663",
        "NCBIGene:9369",
        "NCBIGene:4088",
        "NCBIGene:10239",
        "NCBIGene:9055",
        "NCBIGene:79068",
        "NCBIGene:6137",
        "NCBIGene:7157",
        "NCBIGene:6720",
        "NCBIGene:6928",
        "NCBIGene:8877",
        "NCBIGene:9021",
        "NCBIGene:6925",
        "NCBIGene:4160",
        "NCBIGene:3643",
        "NCBIGene:7137",
        "NCBIGene:129787",
        "NCBIGene:338",
        "NCBIGene:2646",
        "NCBIGene:63892",
        "NCBIGene:53335",
        "NCBIGene:5966",
        "NCBIGene:2888",
        "NCBIGene:4760",
        "NCBIGene:3667",
        "NCBIGene:57761",
        "NCBIGene:5335",
        "NCBIGene:3172",
        "NCBIGene:5770",
        "NCBIGene:7022",
        "NCBIGene:7332",
        "NCBIGene:5468",
        "NCBIGene:7325",
        "NCBIGene:7384",
        "NCBIGene:56999",
        "NCBIGene:7018",
        "NCBIGene:285195",
        "NCBIGene:6514",    
        "NCBIGene:10644",
        "NCBIGene:9370",
        "NCBIGene:10296",
        "NCBIGene:7466",
        "NCBIGene:4214",
        "NCBIGene:85480",
        "NCBIGene:6239",
        "NCBIGene:54901",
        "NCBIGene:1432",
        "NCBIGene:7422",
        "NCBIGene:7021",
        "NCBIGene:5314",
        "NCBIGene:222546",
        "NCBIGene:5167",
        "NCBIGene:7128",
        "NCBIGene:387082",
        "NCBIGene:221895",
        "NCBIGene:2645",
        "NCBIGene:5506",
        "NCBIGene:5078",
        "NCBIGene:136259",
        "NCBIGene:79660",
        "NCBIGene:94241",
        "NCBIGene:169026",
        "NCBIGene:7436",
        "NCBIGene:169792",
        "NCBIGene:6194",
        "NCBIGene:140803",
        "NCBIGene:29979",
        "NCBIGene:6812"
    ]

def test_get_gene_nmf_novelty_for_gene_list():
    """
    Test that the gene nmf adaptare novelty function works.
    """
    # initialize
    map_result = {}

    # call method
    map_result = adapter.get_gene_nmf_novelty_for_gene_list(list_input_genes=list_gene_test, log=True)

    # logger
    logger.info("got map result of size: {}".format(len(map_result.get('gene_results'))))

    # test
    assert map_result is not None
    assert map_result.get('gene_results') is not None
    assert len(map_result.get('gene_results')) > 0
    assert len(map_result.get('gene_results')) == len(list_gene_test)


