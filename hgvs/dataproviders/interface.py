# -*- coding: utf-8 -*-
"""Defines the abstract data provider interface

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import abc

import hgvs

from ..decorators.lru_cache import lru_cache, LEARN, RUN, VERIFY
from ..utils.PersistentDict import PersistentDict
from six.moves import map
import six


class Interface(six.with_metaclass(abc.ABCMeta, object)):
    """Variant mapping and validation requires access to external data,
    specifically exon structures, transcript alignments, and protein
    accessions.  In order to isolate the hgvs package from the myriad
    choices and tradeoffs, these data are provided through an
    implementation of the (abstract) HGVS Data Provider Interface.

    As of June 2014, the only available data provider implementation
    uses the `Universal Transcript Archive (UTA)`_, a sister project
    that provides access to transcripts and genome-transcript
    alignments.  `Invitae`_ provides a public UTA database instance
    that is used by default; see the `UTA`_ page for instructions on
    installing your own PostgreSQL or SQLite version.  In the future,
    other implementations may be availablefor other data sources.

    Pure virtural class for the HGVS Data Provider Interface.  Every
    data provider implementation should be a subclass (possibly
    indirect) of this class.

    .. _`Universal Transcript Archive (UTA)`: https://github.com/biocommons/uta
    .. _Invitae: http://invitae.com/
    """

    def interface_version(self):
        return 4

    def __init__(self, mode=None, cache=None):
        """
        :param mode: cache mode (None[default lru cache], 'learn', 'run', 'verify')
        :type mode: str
        :param cache: local cache file name
        :type cache: str
        """
        self.mode = None
        if mode == 'learn':
            self.mode = LEARN
        elif mode == 'run':
            self.mode = RUN
        elif mode == 'verify':
            self.mode = VERIFY

        self.cache = None
        if self.mode is not None:
            if self.mode == LEARN:
                self.cache = PersistentDict(cache, flag='c')
            else:
                self.cache = PersistentDict(cache, flag='r')

        self.data_version = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.data_version)
        self.schema_version = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.schema_version)
        self.get_acs_for_protein_seq = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_acs_for_protein_seq)
        self.get_gene_info = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_gene_info)
        self.get_pro_ac_for_tx_ac = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_pro_ac_for_tx_ac)
        self.get_seq = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_seq)
        self.get_similar_transcripts = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_similar_transcripts)
        self.get_tx_exons = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_tx_exons)
        self.get_tx_for_gene = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_tx_for_gene)
        self.get_tx_for_region = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_tx_for_region)
        self.get_tx_identity_info = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_tx_identity_info)
        self.get_tx_info = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_tx_info)
        self.get_tx_mapping_options = lru_cache(
            maxsize=hgvs.global_config.lru_cache.maxsize, mode=self.mode,
            cache=self.cache)(self.get_tx_mapping_options)

        def _split_version_string(v):
            versions = list(map(int, v.split(".")))
            if len(versions) < 2:
                versions += [0]
            assert len(versions) == 2
            return versions

        assert self.required_version is not None

        rv = _split_version_string(self.required_version)
        av = _split_version_string(self.schema_version())

        if av[0] == rv[0] and av[1] >= rv[1]:
            return

        raise RuntimeError(
            "Incompatible versions: {k} requires schema version {rv}, but {self.url} provides version {av}"
            .format(
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
    def get_acs_for_protein_seq(self, seq):
        pass

    @abc.abstractmethod
    def get_assembly_map(self, assembly_name):
        pass

    @abc.abstractmethod
    def get_gene_info(self, gene):
        pass

    @abc.abstractmethod
    def get_pro_ac_for_tx_ac(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_seq(self, ac, start_i=None, end_i=None):
        pass

    @abc.abstractmethod
    def get_similar_transcripts(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        pass

    @abc.abstractmethod
    def get_tx_for_gene(self, gene):
        pass

    @abc.abstractmethod
    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        pass

    @abc.abstractmethod
    def get_tx_identity_info(self, tx_ac):
        pass

    @abc.abstractmethod
    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        pass

    @abc.abstractmethod
    def get_tx_mapping_options(self, tx_ac):
        pass


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
