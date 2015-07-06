# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import abc


class Interface(object):
    """Variant mapping and validation requires access to external data,
    specifically exon structures, transcript alignments, and protein
    accessions.  In order to isolate the hgvs package from the myriad
    choices and tradeoffs, these data are provided through an
    implementation of the (abstract) HGVS Data Provider Interface.
    
    As of June 2014, the only available data provider implementation
    uses the Universal Transcript Archive (`UTA`_), a sister project
    that provides access to transcripts and genome-transcript
    alignments.  `Invitae`_ provides a public UTA database instance
    that is used by default; see the `UTA`_ page for instructions on
    installing your own PostgreSQL or SQLite version.  In the future,
    other implementations may be availablefor other data sources.

    Pure virtural class for the HGVS Data Provider Interface.  Every
    data provider implementation should be a subclass (possibly
    indirect) of this class.

    .. _UTA: http://bitbucket.org/biocommons/uta
    .. _Invitae: http://invitae.com/
    """

    __metaclass__ = abc.ABCMeta

    def interface_version(self):
        return 1

    def __init__(self):
        def _split_version_string(v):
            versions = map(int, v.split("."))
            if len(versions) < 2:
                versions += [0]
            assert len(versions) == 2
            return versions

        assert self.required_version is not None

        rv = _split_version_string(self.required_version)
        av = _split_version_string(self.schema_version())

        if av[0] == rv[0] and av[1] >= rv[1]:
            return True

        raise RuntimeError(
            "Incompatible versions: {k} requires schema version {rv}, but {self.url} provides version {av}".format(
                k=type(self).__name__,
                self=self,
                rv=self.required_version,
                av=self.schema_version()))

    # required_version: what version of the remote schema is required
    # by the subclass? This value is compared to the result of
    # schema_version, which must be implemented by the class.
    required_version = None

    @abc.abstractmethod
    def data_version(self):
        pass

    @abc.abstractmethod
    def schema_version(self):
        pass

    @abc.abstractmethod
    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        pass

    @abc.abstractmethod
    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        pass

    @abc.abstractmethod
    def get_tx_seq(self, ac):
        pass

    @abc.abstractmethod
    def get_tx_for_gene(self, gene):
        pass

    @abc.abstractmethod
    def get_acs_for_protein_seq(seq):
        pass

    @abc.abstractmethod
    def get_gene_info(self, gene):
        pass

    @abc.abstractmethod
    def get_tx_mapping_options(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_tx_identity_info(self, tx_ac):
        pass

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
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
