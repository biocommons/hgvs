# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals


from hgvs.decorators.lru_cache import lru_cache
from hgvs.exceptions import HGVSDataNotAvailableError
from hgvs.dataproviders.interface import Interface


import urllib2
import json


class Client(Interface):
    """REST client to the UTA database
    """

    def __init__(self, url='api.biocommons.org', api_version=0):
        if url.find('://') == -1:
            url = 'http://' + url
        self.prefix = url + '/hgvs/v' + str(api_version) + '/'


    def _get_response(self, url):
        try:
            res = urllib2.urlopen(url)
            return json.load(res)
        except urllib2.HTTPError:
            raise HGVSDataNotAvailableError(url)


    @lru_cache(maxsize=1)
    def data_version(self):
        url = self.prefix + 'data_version'
        res = self._get_response(url)
        return res['data_version']


    @lru_cache(maxsize=1)
    def schema_version(self):
        url = self.prefix + 'schema_version'
        res = self._get_response(url)
        return res['schema_version']


    @lru_cache(maxsize=128)
    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        url = '{0}tx_exons?tx_ac={1}&alt_ac={2}&alt_aln_method={3}'.format(self.prefix, tx_ac, alt_ac, alt_aln_method)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        url = '{0}tx_info?tx_ac={1}&alt_ac={2}&alt_aln_method={3}'.format(self.prefix, tx_ac, alt_ac, alt_aln_method)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def fetch_seq(self, ac, start_i=None, end_i=None):
        if start_i and end_i:
            url = '{0}sequence?ac={1}&start={2}&end={3}'.format(self.prefix, ac, start_i, end_i)
        else:
            url = '{0}sequence?ac={1}'.format(self.prefix, ac)
        res = self._get_response(url)
        return res['seq']
    
    @lru_cache(maxsize=128)
    def get_tx_seq(self, ac):
        url = '{0}sequence?ac={1}'.format(self.prefix, ac)
        res = self._get_response(url)
        return res

    @lru_cache(maxsize=128)
    def get_tx_for_gene(self, gene):
        url = '{0}tx_for_gene?gene={1}'.format(self.prefix, gene)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        url = '{0}tx_for_region?alt_ac={1}&alt_aln_method={2}&start={3}&end={4}'.format(self.prefix, alt_ac, alt_aln_method, start_i, end_i)
        res = self._get_response(url)
        return res

    @lru_cache(maxsize=128)
    def get_acs_for_protein_seq(self, seq):
        url = '{0}acs_for_protein_seq?seq={1}'.format(self.prefix, seq)
        res = self._get_response(url)
        return [item['ac'] for item in res]


    @lru_cache(maxsize=128)
    def get_gene_info(self, gene):
        url = '{0}gene_info?gene={1}'.format(self.prefix, gene)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def get_tx_mapping_options(self, tx_ac):
        url = '{0}tx_mapping_options?tx_ac={1}'.format(self.prefix, tx_ac)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def get_tx_identity_info(self, tx_ac):
        url = '{0}tx_identity_info?tx_ac={1}'.format(self.prefix, tx_ac)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def get_similar_transcripts(self, tx_ac):
        url = '{0}similar_transcripts?tx_ac={1}'.format(self.prefix, tx_ac)
        res = self._get_response(url)
        return res


    @lru_cache(maxsize=128)
    def get_pro_ac_for_tx_ac(self, tx_ac):
        url = '{0}pro_ac_for_tx_ac?tx_ac={1}'.format(self.prefix, tx_ac)
        res = self._get_response(url)
        return res['pro_ac']




## <LICENSE>
## Copyright 2015 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
