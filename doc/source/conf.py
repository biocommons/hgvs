############################################################################
## Theme setup

html_theme = 'invitae'
#html_theme = 'f6'
#html_theme = 'bootstrap'
#html_theme = 'sphinx_rtd_theme'
#html_theme = 'pyramid'

html_theme_path = ['../themes']
if html_theme == 'sphinx_rtd_theme':
    import sphinx_rtd_theme
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
elif html_theme == 'bootstrap':
    import sphinx_bootstrap_theme
    html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

print(html_theme, html_theme_path)


############################################################################
## Project config

import hgvs
version = hgvs.__version__
release = str(hgvs.__version__)

project = u'HGVS'
authors = project + ' Contributors'
copyright = u'2014, ' + authors

extlinks={
    'issue': ('https://bitbucket.org/invitae/hgvs/issue/%s', 'HGVS issue '),
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


############################################################################
## Boilerplate

autodoc_default_flags = ['members', 'undoc-members', 'show-inheritance'] #, 'inherited-members']
exclude_patterns = ['build','static','templates','themes']
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
html_favicon = '../static/favicon.ico'
html_logo = '../static/logo.svg'
html_static_path = ['../static']
html_title = '{project} {release}'.format(project = project, release = release)
intersphinx_mapping = {
    'http://docs.python.org/': None,
    }
master_doc = 'index'
pygments_style = 'sphinx'
source_suffix = '.rst'
templates_path = ['templates']
