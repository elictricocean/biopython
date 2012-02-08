#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A package to work with the Gene Ontology.

"""

__author__ = 'Chris Lasher'
__email__ = 'chris DOT lasher <AT> gmail DOT com'

import Bio.GO.Parsers.oboparser


def read(handle, format="obo"):
    
    if format=="obo":
        return Bio.GO.Parsers.oboparser.Parser(handle)