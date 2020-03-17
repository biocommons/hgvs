import os
import re

import hgvs

# N.B. In sphinx parlance, version is supposed to be the undecorated
# version (e.g., 1.2 or 1.2.3).  release corresponds to __version__
# (e.g., 1.2.3.post3).
release = str(hgvs.__version__)
version = re.sub(r"\.post\d+$", "", release)

project = u'HGVS'
authors = project + ' Contributors'
copyright = u'2017, ' + authors

extlinks = {'issue': ('https://github.com/biocommons/hgvs/issues/%s', 'issue '), }

man_pages = [('index', project, project + u' Documentation', [project + u' Contributors'], 1)]

# even for local builds
html_theme = 'sphinx_rtd_theme'

#
# Boilerplate

autodoc_default_flags = ['members', 'undoc-members', 'show-inheritance']    #, 'inherited-members']
autoclass_content = 'both'
exclude_patterns = ['build', 'static', 'templates', 'themes']
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.extlinks',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
]
html_favicon = 'static/favicon.ico'
html_logo = 'static/hgvs-logo.png'
html_static_path = ['static']
html_title = '{project} {release}'.format(project=project, release=release)
intersphinx_mapping = {'http://docs.python.org/': None, }
master_doc = 'index'
pygments_style = 'sphinx'
source_suffix = '.rst'
templates_path = ['templates']


# rst_epilog is appended to all rst files
# it's a good place to define global aliases
# If ends in .rst, sphinx will append it to itself :-(
rst_epilog_fn = os.path.join(os.path.dirname(__file__), 'rst_epilog')
rst_epilog = open(rst_epilog_fn).read()

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
