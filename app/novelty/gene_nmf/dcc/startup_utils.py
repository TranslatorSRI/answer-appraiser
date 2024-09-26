
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
from dcc.dcc_utils import get_logger

# constants
# DIR_CONF = "python-flask-server/conf/"
DIR_CONF = "conf/"
FILE_CONF = "startup_files.json"

# variables
matrix_gene_set = None
vector_genes = None
map_conf = None
logger = get_logger(__name__)

# methods
def load_conf(dir=DIR_CONF, log=True):
    '''
    will load the configuration
    '''
    file_conf = dir + FILE_CONF
    json_data = None

    # log
    if log:
        logger.info("loading startup conf file: {}".format(file_conf))

    # load the map
    with open(file_conf) as f:
         json_data = json.load(f)

    # return
    return json_data

def get_conf():
    '''
    will return the configuration
    '''
    if not map_conf:
        map_conf = load_conf()

    # return
    return map_conf


# main
if __name__ == "__main__":
    pass
