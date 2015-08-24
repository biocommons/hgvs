#!/usr/bin/env python

"""Starts UTA REST dataprovider on localhost:15032"""

# TODO: #262 -- make REST port configurable

port = 15032

from hgvs.dataproviders.utarest.server import Server
Server.run(port=port)
