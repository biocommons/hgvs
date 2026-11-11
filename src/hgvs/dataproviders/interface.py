"""Defines the abstract data provider interface"""

import abc
import logging
import os
from typing import ClassVar

import hgvs

from ..decorators.lru_cache import LEARN, RUN, VERIFY, lru_cache
from ..utils.PersistentDict import PersistentDict
from .models import TranscriptInfo

_logger = logging.getLogger(__name__)


class Interface(metaclass=abc.ABCMeta):
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
        return 5

    @classmethod
    def from_v4(cls, provider):
        """Bridge a data provider written against interface version 4.

        Version 4 providers return plain dicts keyed by the UTA column
        names.  The bridge presents such a provider through this
        interface: methods whose contract is unchanged are delegated
        as-is, and methods that now return models convert the dict.

        >>> hdp = Interface.from_v4(legacy_provider)  # doctest: +SKIP

        Caching, schema checking and any extra attributes (``url``,
        ``seqfetcher``, ...) remain those of the wrapped provider.
        """
        return _V4Bridge(provider)

    def __init__(self, mode=None, cache=None):
        """
        :param mode: cache mode (None[default lru cache], 'learn', 'run', 'verify')
        :type mode: str
        :param cache: local cache file name
        :type cache: str
        """
        self.mode = None
        if mode == "learn":
            self.mode = LEARN
        elif mode == "run":
            self.mode = RUN
        elif mode == "verify":
            self.mode = VERIFY

        self.cache = None
        if self.mode is not None:
            if self.mode == LEARN:
                self.cache = PersistentDict(cache, flag="c")
            else:
                self.cache = PersistentDict(cache, flag="r")

        maxsize = hgvs.global_config.lru_cache.maxsize
        if "PYTEST_CURRENT_TEST" in os.environ:
            maxsize = None
            _logger.info("%s: Using unlimited cache size", __file__)

        self.data_version = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.data_version
        )
        self.schema_version = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.schema_version
        )
        self.get_acs_for_protein_seq = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_acs_for_protein_seq
        )
        self.get_gene_info = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_gene_info
        )
        self.get_pro_ac_for_tx_ac = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_pro_ac_for_tx_ac
        )
        self.get_seq = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(self.get_seq)
        self.get_similar_transcripts = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_similar_transcripts
        )
        self.get_tx_exons = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_exons
        )
        self.get_tx_for_gene = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_for_gene
        )
        self.get_tx_for_region = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_for_region
        )
        self.get_tx_identity_info = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_identity_info
        )
        self.get_tx_info = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_info
        )
        self.get_tx_mapping_options = lru_cache(maxsize=maxsize, mode=self.mode, cache=self.cache)(
            self.get_tx_mapping_options
        )

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

        msg = f"Incompatible versions: {type(self).__name__} requires schema version {self.required_version}, but {self.url} provides version {self.schema_version()}"
        raise RuntimeError(msg)

    # required_version: what version of the remote schema is required
    # by the subclass? This value is compared to the result of
    # schema_version, which must be implemented by the class.
    required_version: ClassVar[str | None] = None

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
    def get_tx_info(self, tx_ac: str, alt_ac: str, alt_aln_method: str) -> TranscriptInfo:
        pass

    @abc.abstractmethod
    def get_tx_mapping_options(self, tx_ac):
        pass


class _V4Bridge(Interface):
    """Presents an interface version 4 data provider as version 5.

    See :meth:`Interface.from_v4`.
    """

    def __init__(self, provider):
        # Interface.__init__ is deliberately not called: the wrapped
        # provider already caches its methods and has checked its schema
        # version, and the bridge is a thin translation layer over it.
        self._provider = provider

    def __getattr__(self, name):
        # attributes the interface does not define (url, seqfetcher, mode,
        # ...) are the wrapped provider's
        return getattr(self._provider, name)

    @property
    def required_version(self):
        return self._provider.required_version

    # -- unchanged contract: delegated as-is -------------------------------

    def data_version(self):
        return self._provider.data_version()

    def schema_version(self):
        return self._provider.schema_version()

    def get_acs_for_protein_seq(self, seq):
        return self._provider.get_acs_for_protein_seq(seq)

    def get_assembly_map(self, assembly_name):
        return self._provider.get_assembly_map(assembly_name)

    def get_gene_info(self, gene):
        return self._provider.get_gene_info(gene)

    def get_pro_ac_for_tx_ac(self, tx_ac):
        return self._provider.get_pro_ac_for_tx_ac(tx_ac)

    def get_seq(self, ac, start_i=None, end_i=None):
        return self._provider.get_seq(ac, start_i, end_i)

    def get_similar_transcripts(self, tx_ac):
        return self._provider.get_similar_transcripts(tx_ac)

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        return self._provider.get_tx_exons(tx_ac, alt_ac, alt_aln_method)

    def get_tx_for_gene(self, gene):
        return self._provider.get_tx_for_gene(gene)

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        return self._provider.get_tx_for_region(alt_ac, alt_aln_method, start_i, end_i)

    def get_tx_identity_info(self, tx_ac):
        return self._provider.get_tx_identity_info(tx_ac)

    def get_tx_mapping_options(self, tx_ac):
        return self._provider.get_tx_mapping_options(tx_ac)

    # -- changed contract: dict results become models ------------------------

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        result = self._provider.get_tx_info(tx_ac, alt_ac, alt_aln_method)
        # version 4 allowed None for "not found"; callers check for it
        return None if result is None else TranscriptInfo(**result)


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
