import os
from xml.etree import ElementTree
import subprocess
# Stuff for monkey-patching
import inspect
import ast
import exhale
import exhale.graph
import exhale.utils
import exhale.configs
import exhale.parse
import breathe
import breathe.directives
import breathe.renderer
import breathe.directives.content_block
import breathe.renderer.sphinxrenderer
from docutils.parsers.rst.directives import flag
from docutils.nodes import Node
from docutils import nodes
from typing import List, cast

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
copyright = '1996-2025, David M. Levine, Philip L. Hallstrom, David M. Noelle, Brian P. Walenz, Dirk Eddelbuettel, Ralf Schlatterbeck'
author = 'David M. Levine, Philip L. Hallstrom, David M. Noelle, Brian P. Walenz, Dirk Eddelbuettel, Ralf Schlatterbeck'

# On readthedocs we need to run doxygen first

read_the_docs_build = os.environ.get ('READTHEDOCS', None) == 'True'
if read_the_docs_build:
    subprocess.call ('doxygen', shell=True)

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'breathe',
    'exhale',
    'sphinxcontrib.inkscapeconverter',
    #'sphinxfortran.fortran_domain'
]

breathe_projects = {
    "PGAPack": "./_doxygen/xml"
}
breathe_default_project = "PGAPack"
breathe_implementation_filename_extensions = []
breathe_order_parameters_first = True
breathe_domain_by_extension = dict (c = 'c', h = 'c')

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
    "exhaleExecutesDoxygen": False,
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
        ALIASES  = "rst=\\verbatim embed:rst"
        ALIASES += "endrst=\\endverbatim"
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
# Compute functions in some groups for which we do not create manpages
if ldir:
    for g in ('internal', 'linalg', 'utility', 'notimplemented', 'fun-bit'):
        tree = ElementTree.parse (os.path.join (dir, 'group__%s.xml' % g))
        root = tree.getroot ()
        for memb in root.findall ('.//memberdef'):
            if memb.get ('kind') != 'function':
                continue
            group_remove [memb.get ('id')] = True
for fn in ldir:
    if fn.endswith ('8c.xml') or fn.endswith ('8h.xml'):
        tree = ElementTree.parse (os.path.join (dir, fn))
        root = tree.getroot ()
        for memb in root.findall ('.//memberdef'):
            if memb.get ('kind') != 'function':
                continue
            name  = memb.find ('name').text
            # Don't allow static, but allow static inline
            if memb.get ('static') == 'yes' and memb.get ('inline') == 'no':
                continue
            id = memb.get ('id')
            if id in group_remove:
                continue
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

def src (obj):
    s = inspect.getsource (obj).split ('\n')
    indent = len (s [0]) - len (s [0].strip ())
    return '\n'.join (i [indent:] for i in s)
# end def src

def monkey_patch ():
    cls = exhale.graph.ExhaleRoot
    d   = dict \
        ( configs  = exhale.configs
        , utils    = exhale.utils
        , StringIO = exhale.graph.StringIO
        , textwrap = exhale.graph.textwrap
        , codecs   = exhale.graph.codecs
        , parse    = exhale.parse
        )

    n   = 'generateViewHierarchies'
    fun = getattr (cls, n)
    mod = ast.parse (src (fun))
    assert mod.body [0].body [4].value.args [0].keys [3].value == 'file_title'
    assert (  mod.body [0].body [4].value.args [0].values [3].value
           == 'Class Hierarchy'
           )
    mod.body [0].body [4].value.args [0].values [3].value = 'Data Structures'
    exec (compile (mod, '<string>', 'exec'), d)
    setattr (cls, n, d [n])

    n   = 'generateNamespaceChildrenString'
    fun = getattr (cls, n)
    mod = ast.parse (src (fun))
    assert mod.body [0].body [11].value.args [1].value == 'Classes'
    mod.body [0].body [11].value.args [1].value = 'Structs'
    exec (compile (mod, '<string>', 'exec'), d)
    setattr (cls, n, d [n])

    n   = 'generateFileNodeDocuments'
    fun = getattr (cls, n)
    mod = ast.parse (src (fun))
    assert mod.body[0].body [2].body [15].value.args [1].value == 'Classes'
    mod.body[0].body [2].body [15].value.args [1].value = 'Structs'
    exec (compile (mod, '<string>', 'exec'), d)
    setattr (cls, n, d [n])

    n   = 'generateUnabridgedAPI'
    fun = getattr (cls, n)
    mod = ast.parse (src (fun))
    assert (  mod.body [0].body [3].body [5].value.elts [1].elts [0].value
           == 'Classes and Structs'
           )
    mod.body [0].body [3].body [5].value.elts [1].elts [0].value = 'Structs'
    exec (compile (mod, '<string>', 'exec'), d)
    setattr (cls, n, d [n])

    cls = breathe.directives.content_block.DoxygenGroupDirective
    cls.option_spec.update (sort = flag)

    cls = breathe.renderer.sphinxrenderer.SphinxRenderer
    d   = dict \
        ( List = List, Node = Node, cast = cast, nodes = nodes
        , RenderContext = breathe.renderer.RenderContext
        )
    n   = 'visit_sectiondef'
    fun = getattr (cls, n)
    txt = src (fun).split ('\n')
    assert txt [6].lstrip ().startswith ('# Get all the memberdef info')
    txt.insert (7, "    if 'sort' in options:")
    txt.insert (8, "        node.memberdef.sort(key=lambda x: x.name)")
    mod = ast.parse ('\n'.join (txt))
    exec (compile (mod, '<string>', 'exec'), d)
    setattr (cls, n, d [n])
    cls.methods ['sectiondef'] = d [n]
# end def monkey_patch

monkey_patch ()

def get_popen (cmd):
    """ Return stripped output from external command
    """
    fd = os.popen (cmd)
    ret = fd.read ().strip ()
    fd.close ()
    return ret
# end def get_popen

today = get_popen ('git show -s --format=%cd --date=format:%Y-%m-%d')
