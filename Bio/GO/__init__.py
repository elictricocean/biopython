#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A package to work with the Gene Ontology.

"""

__author__ = 'Chris Lasher'
__email__ = 'chris DOT lasher <AT> gmail DOT com'

import Bio.GO.Parsers.oboparser
import Bio.GO.Parsers.annotation

def read(handle, format="obo"):
    
    if format=="obo":
        return Bio.GO.Parsers.oboparser.Parser(handle)
    
    if format=="goa":
        return Bio.GO.Parsers.annotation.Parser(handle)