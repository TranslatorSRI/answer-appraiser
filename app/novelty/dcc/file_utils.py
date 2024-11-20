
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
import json

import dcc.dcc_utils as dutils 


# constants
logger = dutils.get_logger(__name__)


# methods
def read_tab_delimited_file(file_path, log=False):
    '''
    will read the pathway file and build a map of gene lists
    '''
    # initialize
    result = {}

    # read in the data
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            fields = line.strip().split('\t')
            stripped_list = [item.strip() for item in fields]

            key = fields[0]
            values = stripped_list[1:]
            if key in result:
                result[key].extend(values)
            else:
                result[key] = values

    # return
    return result


def load_gene_file_into_map(file_path, log=False):
    '''
    will read the gene file and return a map of gene to position in array
    '''
    # DEPRECATED - now get it from sqlite db
    # initialize
    map_result = {}
    list_temp = []

    # read in the data
    with open(file_path, 'r') as file:
        for line in file:
            gene = line.strip()
            list_temp.append(gene)
            # result[gene] = count
            # count = count + 1

    # create the map from the list; make sure no duplicates
    list_unique_gene = list(set(list_temp))
    list_unique_gene.sort()
    map_result = {value: index for index, value in enumerate(list_unique_gene)}

    # log
    logger.info("loaded gene file: {} into num count map: {}".format(file_path, len(map_result)))

    # return
    return map_result, list_unique_gene


def get_all_files_in_dir(path_dir, log=False):
    '''
    will return the files in the directory 
    '''
    list_result = []

    # get the files

    # return
    return list_result
    





if __name__ == "__main__":

    file_path = '/home/javaprog/Code/TranslatorWorkspace/PigeanFlask/python-flask-server/conf/gene_lists/gene_set_list_msigdb_h.txt'
    data = read_tab_delimited_file(file_path)
    print(json.dumps(data, indent=2))