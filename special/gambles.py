#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 15:53:56 2018

@author: benavoli
"""

from sympy import Piecewise
from sympy.parsing.sympy_parser import parse_expr

def indicator(condition):
    """Returns 1 if condition 0 otherwise"""
    expr=parse_expr(condition)
    return Piecewise((1,expr),(0,True)) 


def conditionalGamble(g,condition):
    """Returns a conditional gamble"""
    return g*indicator(condition)



def conditionalGambleList(G,condition):
    """Returns list of conditional gambles"""
    Gnew=[]
    for i in range(0,len(G)):
        Gnew.append(conditionalGamble(G[i],condition))        
    return Gnew