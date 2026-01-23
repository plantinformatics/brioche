# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Brioche'
copyright = "2026, James O'Dwyer, Gabriel Keeble-Gagnere, Kerrie L Forrest"
author = "James O'Dwyer, Gabriel Keeble-Gagnere, Kerrie L Forrest"
release = '2.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']


html_theme_options = {
    "github_user": "plantinformatics",
    "github_repo": "brioche",
    "github_banner": True,
    "show_powered_by": False,
}
