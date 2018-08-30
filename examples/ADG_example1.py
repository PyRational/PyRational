#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 09:15:36 2018

@author: benavoli
"""

from BeliefModel import  ADG
from sympy import symbols,  Interval, FiniteSet, Piecewise, Eq




x1 = symbols('x1');
domain_x1=Interval(-5,5)

x2 = symbols('x2');
domain_x2=FiniteSet(0,1)


ListSymbols=[x1]
ListDomains=[domain_x1]
#
model = ADG(ListSymbols,ListDomains)
model.add_gambleList([x1-1,1-x1,x1**2-1.1,x1**4-1])

model.buildModel()




print(model)

optimoptions={'method_LISP': 'Cutting_plane', 
                   'SolverLP':'cplex',
                   'LP_acc_constraints':1e-8,
                   'SolverNLP':'linprog',
                   'NLP_acc_constraints':1e-3,
                   'num_support_points': 5,
                   'verbose':True}

#res=model.check_avs(options=optimoptions)

NewG=model.updating(sympy.Piecewise((1.0,True)),options=optimoptions)

#c=0.5
#f=Piecewise((1,x1>=c),(0,x1<c))
#f_range=(0,1)
#res=model.lower_prevision(f,f_range, options=optimoptions)

#print(res)
#c=0.5
#g1 = Piecewise((1,x1>=c),(0,x1<c))-0.5
#g2=x1-0.5
#model.add_gambleList([g1,g2])
#print(model)

