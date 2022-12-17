import os
from xml.etree import ElementTree

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'pgapack'
copyright = '1996-2022, David M. Levine, Philip L. Hallstrom, David M. Noelle, Brian P. Walenz, Dirk Eddelbuettel, Ralf Schlatterbeck'
author = 'David M. Levine, Philip L. Hallstrom, David M. Noelle, Brian P. Walenz, Dirk Eddelbuettel, Ralf Schlatterbeck'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'breathe',
    'exhale'
]

breathe_projects = {
    "pgapack": "./_doxygen/xml"
}
breathe_default_project = "pgapack"
breathe_implementation_filename_extensions = []
breathe_order_parameters_first = True

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "../..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    """
        INPUT = ../../include ../../source
        XML_PROGRAMLISTING     = NO
        OPTIMIZE_OUTPUT_FOR_C  = YES
        EXTRACT_ALL            = NO
        EXTRACT_LOCAL_CLASSES  = NO
        SHOW_INCLUDE_FILES     = NO
        SHOW_NAMESPACES        = NO
        JAVADOC_AUTOBRIEF      = YES
        SKIP_FUNCTION_MACROS   = YES
        GENERATE_MAN           = YES
        GENERATE_HTML          = YES
        EXPAND_ONLY_PREDEF     = YES
        EXCLUDE = ../../include/pgapackf.h
    """
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'c'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'c'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


man_pages = []
dir = '_doxygen/xml'
authors = ''

try:
    ldir = os.listdir (dir)
except FileNotFoundError:
    ldir = []
for fn in ldir:
    if fn.endswith ('8c.xml'):
        tree = ElementTree.parse (os.path.join (dir, fn))
        root = tree.getroot ()
        for memb in root.findall ('.//memberdef'):
            if memb.get ('kind') != 'function':
                continue
            id    = memb.get ('id')
            name  = memb.find ('name').text
            kind  = memb.get ('kind')
            try:
                brief = memb.find ('briefdescription') [0].text
            except IndexError:
                brief = ''
            man_pages.append \
                (('api/%s_%s' % (kind, id), name, brief, authors, '3'))

man_make_section_directory = True
