# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

from hgvs.dataproviders.utarest import apis
from hgvs.dataproviders.utarest import api

api_version = '0'

apis.add_resource(api.DataVersion,        '/hgvs/v'+api_version+'/data_version')
apis.add_resource(api.SchemaVersion,      '/hgvs/v'+api_version+'/schema_version')
apis.add_resource(api.TxExons,            '/hgvs/v'+api_version+'/tx_exons')
apis.add_resource(api.TxInfo,             '/hgvs/v'+api_version+'/tx_info')
apis.add_resource(api.FetchSeq,           '/hgvs/v'+api_version+'/sequence')
apis.add_resource(api.TxForGene,          '/hgvs/v'+api_version+'/tx_for_gene')
apis.add_resource(api.TxForRegion,        '/hgvs/v'+api_version+'/tx_for_region')
apis.add_resource(api.AcsForProteinSeq,   '/hgvs/v'+api_version+'/acs_for_protein_seq')
apis.add_resource(api.GeneInfo,           '/hgvs/v'+api_version+'/gene_info')
apis.add_resource(api.TxMappingOptions,   '/hgvs/v'+api_version+'/tx_mapping_options')
apis.add_resource(api.TxIdentityInfo,     '/hgvs/v'+api_version+'/tx_identity_info')
apis.add_resource(api.SimilarTranscripts, '/hgvs/v'+api_version+'/similar_transcripts')
apis.add_resource(api.ProAcForTxAc,       '/hgvs/v'+api_version+'/pro_ac_for_tx_ac')


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