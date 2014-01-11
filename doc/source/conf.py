import sys, os

import hgvs
version = hgvs.__version__
release = str(hgvs.__version__)

project = u'HGVS'
copyright = u'2014, InVitae'

source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['build','static','templates','themes']
pygments_style = 'sphinx'
templates_path = ['templates']

html_theme = 'invitae'
html_theme_path = ['../themes']
html_title = '{project} {release}'.format(project = project, release = release)
html_logo = '../static/hgvs-logo.svg'
html_favicon = '../static/favicon.ico'
html_static_path = ['../static']

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.intersphinx',
    'sphinx.ext.pngmath',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinxcontrib.fulltoc',
    ]

extlinks={
    'issue': ('https://bitbucket.org/invitae/hgvs/issue/%s' 'HGVS issue '),
    }

intersphinx_mapping = {
    'http://docs.python.org/': None,
    }

latex_documents = [
  ('index', 'HGVS.tex', u'HGVS Documentation', u'HGVS Contributors', 'manual'),
]

man_pages = [
    ('index', 'uta', u'HGVS Documentation', [u'HGVS Contributors'], 1)
]

texinfo_documents = [
  ('index', 'HGVS', u'HGVS Documentation',
   u'HGVS Contributors', 'HGVS', 'One line description of project.',
   'Miscellaneous'),
]
