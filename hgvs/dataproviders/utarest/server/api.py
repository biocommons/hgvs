# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import os

from flask import g
from flask_restful import Resource
from flask_restful import reqparse
from flask_restful import fields
from flask_restful import marshal

from hgvs.dataproviders.utarest.server import Server

import hgvs.dataproviders.uta
from hgvs.exceptions import HGVSError, HGVSDataNotAvailableError


_default_db_url = os.environ.get('UTA_DB_URL', hgvs.dataproviders.uta._uta_urls['public'])
if _default_db_url.startswith('http://') or _default_db_url.startswith('https://'):
    _default_db_url = hgvs.dataproviders.uta._uta_urls['public']



@Server.before_request
def before_request():
    g.dp = hgvs.dataproviders.uta.connect(db_url=_default_db_url)





class DataVersion(Resource):
    def get(self):
        return {'data_version' : g.dp.data_version()}


class SchemaVersion(Resource):
    def get(self):
        return {'schema_version' : g.dp.schema_version()}


class TxExons(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        
        self.reqparse.add_argument('tx_ac', type = str, required = True,
                        help = 'No tx_ac provided', location = 'args')
        self.reqparse.add_argument('alt_ac', type = str, required = True,
                        help = 'No alt_ac provided', location = 'args')
        self.reqparse.add_argument('alt_aln_method', type = str, required = True,
                        help = 'No alt_aln_method provided', location = 'args')
        
        self.resource_fields = {
            'tes_exon_set_id' : fields.Integer(attribute='tx_exon_set_id'),
            'aes_exon_set_id' : fields.Integer(attribute='alt_exon_set_id'),
            'tx_ac'           : fields.String,
            'alt_ac'          : fields.String,
            'alt_strand'      : fields.Integer,
            'alt_aln_method'  : fields.String,
            'ord'             : fields.Integer,
            'tx_exon_id'      : fields.Integer,
            'alt_exon_id'     : fields.Integer,
            'tx_start_i'      : fields.Integer,
            'tx_end_i'        : fields.Integer,
            'alt_start_i'     : fields.Integer,
            'alt_end_i'       : fields.Integer,
            'cigar'           : fields.String,
        }
        
        super(TxExons, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_tx_exons(args['tx_ac'], args['alt_ac'], args['alt_aln_method'])
            res = [dict(item) for item in res]
            return marshal(res, self.resource_fields)
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class TxInfo(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        
        self.reqparse.add_argument('tx_ac', type = str, required = True,
                        help = 'No tx_ac provided', location = 'args')
        self.reqparse.add_argument('alt_ac', type = str, required = True,
                        help = 'No alt_ac provided', location = 'args')
        self.reqparse.add_argument('alt_aln_method', type = str, required = True,
                        help = 'No alt_aln_method provided', location = 'args')
        
        super(TxInfo, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_tx_info(args['tx_ac'], args['alt_ac'], args['alt_aln_method'])
            return dict(res)
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404
        except HGVSError as e:
            return {'error': str(e)}, 400


class FetchSeq(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('ac', type = str, required = True,
                        help = 'No ac provided', location = 'args')
        self.reqparse.add_argument('start', type = int, required = False,
                        help = 'No start provided', location = 'args', default=None)
        self.reqparse.add_argument('end', type = int, required = False,
                        help = 'No end provided', location = 'args', default=None)
        super(FetchSeq, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.fetch_seq(args['ac'], args['start'], args['end'])
            return {'seq' : res}
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class TxForGene(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('gene', type = str, required = True,
                        help = 'No gene provided', location = 'args')
        super(TxForGene, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_tx_for_gene(args['gene'])
            res = [dict(item) for item in res]
            return res
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class TxForRegion(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('alt_ac', type = str, required = True,
                        help = 'No alt_ac provided', location = 'args')
        self.reqparse.add_argument('alt_aln_method', type = str, required = True,
                        help = 'No alt_aln_method provided', location = 'args')
        self.reqparse.add_argument('start', type = str, required = True,
                        help = 'No start provided', location = 'args')
        self.reqparse.add_argument('end', type = str, required = True,
                        help = 'No end provided', location = 'args')
        super(TxForRegion, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_tx_for_region(args['alt_ac'], args['alt_aln_method'], args['start'] , args['end'])
            res = [dict(item) for item in res]
            return res
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class AcsForProteinSeq(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('seq', type = str, required = True,
                        help = 'No seq provided', location = 'args')
        super(AcsForProteinSeq, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_acs_for_protein_seq(args['seq'])
            res = [{'ac' : item} for item in res]
            return res
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class GeneInfo(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        
        self.reqparse.add_argument('gene', type = str, required = True,
                        help = 'No gene provided', location = 'args')
        
        self.resource_fields = {
            'hgnc'    : fields.String,
            'maploc'  : fields.String,
            'descr'   : fields.String,
            'summary' : fields.String,
            'aliases' : fields.String,
            'added'   : fields.DateTime,
        }
        
        super(GeneInfo, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_gene_info(args['gene'])
            return marshal(dict(res), self.resource_fields)
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class TxMappingOptions(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('tx_ac', type = str, required = True,
                        help = 'No tx_ac provided', location = 'args')
        super(TxMappingOptions, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_tx_mapping_options(args['tx_ac'])
            res = [dict(item) for item in res]
            return res
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class TxIdentityInfo(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('tx_ac', type = str, required = True,
                        help = 'No tx_ac provided', location = 'args')
        super(TxIdentityInfo, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_tx_identity_info(args['tx_ac'])
            return dict(res)
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class SimilarTranscripts(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('tx_ac', type = str, required = True,
                        help = 'No tx_ac provided', location = 'args')
        super(SimilarTranscripts, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_similar_transcripts(args['tx_ac'])
            res = [dict(item) for item in res]
            return res
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404


class ProAcForTxAc(Resource):
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('tx_ac', type = str, required = True,
                        help = 'No tx_ac provided', location = 'args')
        super(ProAcForTxAc, self).__init__()
    
    def get(self):
        try:
            args = self.reqparse.parse_args()
            res = g.dp.get_pro_ac_for_tx_ac(args['tx_ac'])
            if res:
                return {'pro_ac' : res}
            else:
                raise HGVSDataNotAvailableError('No pro_ac found')
        except HGVSDataNotAvailableError as e:
            return {'error': str(e)}, 404




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