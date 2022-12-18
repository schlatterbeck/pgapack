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
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'PGAPack'
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
    "PGAPack": "./_doxygen/xml"
}
breathe_default_project = "PGAPack"
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
        EXTRACT_STATIC         = NO
        SHOW_INCLUDE_FILES     = NO
        SHOW_NAMESPACES        = NO
        QT_AUTOBRIEF           = YES
        JAVADOC_AUTOBRIEF      = YES
        SKIP_FUNCTION_MACROS   = YES
        GENERATE_MAN           = YES
        GENERATE_HTML          = YES
        EXPAND_ONLY_PREDEF     = YES
        EXCLUDE = ../../include/pgapackf.h ../../source/mpi_stub.c
        PREDEFINED += STATIC=static
        PREDEFINED += PGADebugEntered(a)
        PREDEFINED += PGADebugExited(a)
        PREDEFINED += PGACheckDataType(a,b)
        PREDEFINED += INDEX(a,b,c,d)
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
html_theme = 'sphinx_rtd_theme'

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
group_remove = {}
# Compute functions in group internal and explicit
# For these we do not create manpages
if ldir:
    for g in 'internal', 'explicit':
        tree = ElementTree.parse (os.path.join (dir, 'group__%s.xml' % g))
        root = tree.getroot ()
        for memb in root.findall ('.//memberdef'):
            if memb.get ('kind') != 'function':
                continue
            group_remove [memb.get ('id')] = True
for fn in ldir:
    if fn.endswith ('8c.xml'):
        tree = ElementTree.parse (os.path.join (dir, fn))
        root = tree.getroot ()
        for memb in root.findall ('.//memberdef'):
            if memb.get ('kind') != 'function':
                continue
            if memb.get ('static') == 'yes':
                continue
            id = memb.get ('id')
            if id in group_remove:
                continue
            name  = memb.find ('name').text
            kind  = memb.get ('kind')
            try:
                brief = memb.find ('briefdescription') [0].text
            except IndexError:
                brief = ''
            # Only create manpage if brief description exists
            #if not brief:
            #    continue
            man_pages.append \
                (('api/%s_%s' % (kind, id), name, brief, authors, '3a'))

man_make_section_directory = True
