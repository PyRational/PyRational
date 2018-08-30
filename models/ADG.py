#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADG

Implements the belief model ADG: ALmost Desirable Gambles.

Created on Thu Aug 16 08:34:59 2018

@author: benavoli
"""
import sympy as sympy
from sympy import ProductSet
import numpy as np
from . import AbstractModels
import pyDOE as pyDOE
from scipy.optimize import linprog,differential_evolution




class ADG(AbstractModels.BeliefModel):
    """
    ADG Belief Model
    
    :param Symbols: List of Symbols (es. ['x1','x2',...])
    :param Domains: List of Domains
    :rtype: model object

    """
    _max_number_of_discrete_elements=10000 #max number of elements for powerset of discrete variables
    optimoptions={'method_LISP': 'Cutting_pane', #'Cutting_plane', 'discretise'
                       'SolverLP':'linprog', #'linprog', 'cplex'
                       'LP_acc_constraints':1e-8,
                       'SolverNLP':'differential_evolution',
                       'NLP_acc_constraints':1e-8,
                       'num_support_points': 30,
                       'verbose':False}
    
    def __init__(self, _Symbols, _Domains):
        if type(_Symbols) is not list:
            raise TypeError("Symbols must be a list of sympy symbols")
        if type(_Domains) is not list:
            raise TypeError("Domains must be a list of sympy domains")
        if len(_Symbols)!=len(_Domains):
            raise TypeError("Symbols and Domains must have same length")
            
        self.is_continuous=False
        self.is_mixed=False
        self.is_discrete=False
        
        self.is_built=False
        self.is_AVS=None
        self.is_Coherent=None
        
        self.Symbols=_Symbols#this removes duplicates
        self._Domains=_Domains
        self._TypeSymbols=[]
        for el in _Domains:
            if el.is_FiniteSet==True:
                 self._TypeSymbols.append('d')
            if el.is_Interval==True:
                 self._TypeSymbols.append('c')
                 

    
        ind=np.argsort(self._TypeSymbols)
        self._TypeSymbols=list(np.array(self._TypeSymbols)[ind])
        self.Symbols=list(np.array(self.Symbols)[ind])
        self._Domains=list(np.array(self._Domains)[ind])
        self.Domain=ProductSet(self._Domains)
        if ('c' in self._TypeSymbols and 'd' not in self._TypeSymbols):
            self.is_continuous=True
        elif ('c' not in self._TypeSymbols and 'd' in self._TypeSymbols):
            self.is_discrete=True
        else:
            self.is_mixed=True
            ind,=np.where(np.array(self._TypeSymbols)=='c')
            self.last_continuous=ind[-1]
            
        self.Gambles=[] #list of gambles
        self.SupportingPoints=[] #Supporting point of the discretisation
        self.DiscretePowerSet=[] #power set of discrete points 
        self.lambda_bounds=[] #bounds for lambda_i in LP
        self.type_relation=[] #type of relation: '>=' inequality or '==' equality 
        

    def buildModel(self):
        """Builds a BeliefModel"""
        G=self.Gambles.copy()
        NewG=[]
        BND=[]
        TypeR=[]
        while len(G)>0:
            g=G[0]
            delL=[0]
            bnd=[0,None]
            type_relation='>='
            for i in range(1,len(G)):
                val=sympy.Eq(g,G[i])
                if val==True: #removing redundant gambles
                    delL.append(i) 
                val=sympy.Eq(g,-G[i])
                if val==True: #removing g==-g gambles by imposing lambda in R
                    delL.append(i)
                    bnd[0]=None
                    type_relation='=='
            NewG.append(g)
            BND.append(bnd)
            TypeR.append(type_relation)
            G=[G[i] for i in set(np.arange(0,len(G)))-set(delL)]
        self.Gambles=NewG
        self.lambda_bounds=BND
        self.type_relation=TypeR
        self.is_built=True
    

    def coherence(self,discretisation=np.array([]),options=[], tolerance=1e-4):
         """Checks coherence (consistency) of the BeliefModel"""
         
         res=self.check_avs(discretisation=discretisation,options=options, tolerance=tolerance)
         
         return res
    

    def inference(self,f,f_range,h=sympy.Piecewise((1.0,True)),discretisation=np.array([]),options=[]):
         """Returns the result of an inference on a BeliefModel"""
         res=self.natural_extension(f, f_range,h=h,discretisation=discretisation,options=options)
         return res
         

    def updating(self,h,G=[],G_range=[],discretisation=np.array([]),options=[]):
        """Updates a BeliefModel"""
        if len(G)==len(G_range):
            print("G and G_range must have same length")
        NewG=[]
        if G==[]:
            G=self.Gambles
            G_range=[  [None] * 2  ]*len(G)
            type_relation=self.type_relation             
        else:
             type_relation=['>=']*len(G)
             

        for i in range(0,len(G)):
            g=G[i]
            g_range=G_range[i]
            
             
            sol=self.lower_prevision(g, g_range ,h=h, discretisation=discretisation,options=options)
            NewG.append(g-sol)
            if type_relation[i]=='==':
                sol=self.lower_prevision(-g, [-g_range[1] if g_range[1]!=None else None,-g_range[0] if g_range[0]!=None else None] ,h=h, discretisation=discretisation,options=options)
                NewG.append(-g-sol)
                 
        return NewG
             
             
    def revision(self,h):
        raise NotImplementedError("To be implemented")
        
    def __str__(self):
        """
        Prints BeliefModel
        """
        stri='List of Symbols: '+str(self.Symbols)+'\n'
        stri=stri+'Domain: '+str(self.Domain)+'\n'
        stri=stri+'List of desirable gambles: '+str(self.Gambles)+'\n'
        stri=stri+'Avoiding sure loss? '+str(self.is_AVS)+'\n'
        #stri=stri+'Coherence? '+str(self.is_Coherent)+'\n'
        return stri
    
    def _repr_html_(self):
        stri = "<h4>ADG model</h4>"
        
        stri = stri + "<table width=\"100%\" border=\"3\"  >"
        stri = stri + "<tr> <th bgcolor=\"#FFCC33\" width=\"30%\" style=\"text-align: left;\"> List of Symbols </th>"
        stri = stri + "<td bgcolor=\"#f1f1f1\" style=\"text-align: left;\">" + str(self.Symbols)+ "</th>"
        stri = stri + "</tr>"
        stri = stri + "<tr> <th bgcolor=\"#FFDD33\" width=\"30%\" style=\"text-align: left;\"> Domain </th>"
        stri = stri + "<td bgcolor=\"#f6f6f6\" style=\"text-align: left;\">" + "&Omega;="+str(self.Domain)+ "</th>"
        stri = stri + "</tr>"
        stri = stri + "<tr> <th bgcolor=\"#FFCC33\" width=\"30%\" style=\"text-align: left;\"> List of desirable gambles </th>"  
        if self.is_built==True:
            if len(self.Gambles)==0:
                stri = stri + "<td bgcolor=\"#f1f1f1\" style=\"text-align: left;\">" + "G=&Lscr;(&Omega;)<sup>+</sup>"+"</th>"
            elif len(self.Gambles)>0:
                stri = stri + "<td bgcolor=\"#f1f1f1\" style=\"text-align: left;\">" + "G=posi(&Lscr;(&Omega;)<sup>+</sup>      &cup;    " +str(self.Gambles)+ ")</th>"
        else:
            stri = stri + "<td bgcolor=\"#f1f1f1\" style=\"text-align: left;\">" + "G=" +str(self.Gambles)+ "</th>"
        stri = stri + "</tr>"
        stri = stri + "<tr> <th bgcolor=\"#FFDD33\" width=\"30%\" style=\"text-align: left;\"> Avoiding sure loss? </th>"
        stri = stri + "<td bgcolor=\"#f6f6f6\" style=\"text-align: left;\">" + str(self.is_AVS)+ "</th>"
        stri = stri + "</tr>"
        stri = stri + "</table>"
        return stri
        
        
    def to_dict(self, save_data=True):
        raise NotImplementedError
        
        
    def add_gambleList(self,_G):
        """
        Adds a list of gambles to BeliefModel
        """
        self.is_built=False
        self.is_AVS=False
        if type(_G) is not list:
            raise TypeError("Gambles must be a list of sympy functions")
        for g in _G:
            if not isinstance(g,sympy.Expr):
                print(g)
                TypeError("gamble must be a sympy expression")
            #if len(g.free_symbols-set(self.Symbols))>0:
            #    raise TypeError("unknown Symbol in gamble", g)
            self.Gambles.append(g)
            self.lambda_bounds.append((0,None))
            self.type_relation.append('>=')
            
    def __make_DiscretePoints(self,discretdomain):
        '''
        makes discretisation points for discrete variables
        '''        
        T=[]
        if len(discretdomain)<self._max_number_of_discrete_elements:                
            for d in discretdomain:
                if isinstance(d,tuple):
                    T.append(list(d))
                else:
                    T.append([d])
        else:
            raise TypeError("PowerSet of the domains of discrete variables has too many elements. You can try to increase '_max_number_of_discrete_elements'")
        return T
    
    def __make_ContinuousPoints(self,nvar,npoints, criterion='centermaximin'):
        '''
        makes discretisation points for discrete variables using design of experiment strategies
        '''
        if  criterion=='random':
            doe=np.random.rand(npoints,nvar)
        else:
            doe=pyDOE.lhs(nvar, samples=npoints, criterion=criterion)
            doe_n=doe.copy()
        T=[]
        for i in range(0,npoints):
            t=[]
            for j in range(0,nvar):
                d=self._Domains[j]    
                t.append(doe[i,j]*(d.sup-d.inf)+d.inf)
            T.append(t)                    
        return T
                
    def make_SupportingPoints(self,n, criterion='centermaximin'):
        '''
        makes discretisation points for continuous and discrete vars using design of experiment strategy + exahustive on discrete vars
        '''
        if self.is_discrete==True:
            T=self.__make_DiscretePoints(self.Domain)
        elif self.is_continuous==True:
            T=self.__make_ContinuousPoints(len(self.Symbols), n, criterion)
        elif self.is_mixed==True:
            Tc=self.__make_ContinuousPoints(self.last_continuous+1, n, criterion)
            discretedomain=ProductSet(self._Domains[self.last_continuous+1:])
            Td=self.__make_DiscretePoints(discretedomain)
            self.DiscretePowerSet=Td
            T=[]
            for i in range(0,len(Tc)):
                for j in range(0,len(Td)):
                    T.append(Tc[i]+Td[j])
        #T = np.unique(T,axis=0) #initial discretisation
        return np.array(T)
    
    ############## Inference ##############
    def check_avs(self, discretisation=np.array([]),options=[], tolerance=1e-4):
        """
        Checks if the set of desirable gambles avoids sure loss
        """
        c=-0.5
        f = sympy.Piecewise((c,True))
        f_range=[c,1]
        
        lambda0=self.lower_prevision(f, f_range , discretisation=discretisation,options=options)        
        
        if lambda0>=1-tolerance:
            print("Belief Model incurs in a sure  loss")
            self.is_AVS=False
        elif lambda0<=c+tolerance:
            print("Belief Model avoids sure  loss")
            self.is_AVS=True
        else:
            print("Result is "+str(lambda0)+"is inside ["+str(tolerance)+",1-"+str(tolerance)+"]. No decision")

        return self.is_AVS
    
    
    def natural_extension(self, f, f_range,h=sympy.Piecewise((1.0,True)), discretisation=np.array([]),options=[]):
        """
        Computes natural extension for f
        """
        lambda0=self.lower_prevision( f, f_range, h=h, discretisation=discretisation,options=options)
        res=False
        if lambda0>=0:
            res=True
        return res

        
    def lower_prevision(self, f, f_range ,h=sympy.Piecewise((1.0,True)), discretisation=np.array([]),options=[]):
        """
        Computes lower prevision of f
        """
        if options==[]:
            options=self.optimoptions
            
        if len(discretisation)==0 and len(self.SupportingPoints)==0: 
            self.SupportingPoints=self.make_SupportingPoints(options['num_support_points'])            
        elif len(discretisation)>0:
            self.SupportingPoints=discretisation
        
    
        if options['verbose']==True:
            print('Method for LSIP = ',options['method_LISP'],', Linear Programming solver = ',options['SolverLP'], 'NonLinear Programming solver = ',options['SolverNLP'])
        
        sol,T=self.__optimize(self.Symbols,f,f_range,self.Gambles, self.SupportingPoints,h=h,verbose=options['verbose'],method=options['method_LISP'],solver=options['SolverLP'],epsilonLP=options['LP_acc_constraints'],epsilonNLP=options['NLP_acc_constraints'])
        self.SupportingPoints=T
        
        
        return sol['fun']
        
    def upper_prevision(self, f, f_range, h=sympy.Piecewise((1.0,True)), discretisation=np.array([]),options=[]):
        """
        Computes upper prevision of f
        """
        lambda0=self.lower_prevision( -f, (-f_range[1],-f_range[0]),h=h, discretisation=discretisation,options=options)
        return -lambda0
    
    
    ############## Optimisation ##############
    
    def __optimize(self,x,f,f_range,G,T,h=sympy.Piecewise((1.0,True)),verbose=True,method='Cutting_plane',solver='linprog',epsilonLP=1e-8,epsilonNLP=1e-8): 
        if self.is_discrete==True:
            method='discretise'
            if verbose==True:
                print("In this case LSIP reduces to LP")
        if method=='discretise':
            sol=self.LP(x,f,f_range,G, T,h=h,verbose=verbose,solver=solver,epsilon=epsilonLP)
        elif method=='Cutting_plane': #cutting plane
            sol,T=self.Cutting_plane(x,f,f_range,G,T,h=h,verbose=verbose,solver=solver,epsilonLP=epsilonLP,epsilonNLP=epsilonNLP)
        if verbose:
            print('Optimisation Terminated')
        return sol, T

    def LP(self,x,f,f_range,G,T,h=sympy.Piecewise((1.0,True)),verbose=True,solver='linprog',epsilon=1e-8):    
        if solver=='linprog':
            res=self._linprog(x,f,f_range,G,T,h=h,verbose=verbose,epsilon=epsilon)
        elif solver=='cplex':
            res=self._cplex(x,f,f_range,G,T,h=h,verbose=verbose,epsilon=epsilon)
        return res
    
    def populatebycolumn(self,x,f,f_range,G,T,prob,h=sympy.Piecewise((1.0,True)),verbose=False,epsilon=1e-8):
        import cplex as cplex
        nvar=len(G)+1
        ncons=T.shape[0]
        
        
        
        #this is the weighting function (partition of unity)        
        hval=np.array([float(h.subs(list(zip(x,T[j,:])))) for j in range(0,T.shape[0])])
        col=np.multiply(np.ones(ncons),hval)
        
        nonzeroind = np.nonzero(col)[0] 
        #print(np.ndarray.tolist(col[nonzeroind]))
    
        #define matrix A of the linear constraints
        ind_cons=np.arange(0,ncons,1)
        A=[]
        A.append([np.ndarray.tolist(ind_cons[nonzeroind]),np.ndarray.tolist(col[nonzeroind])])#this multiplies lambda_0
        for i in range(0,len(G)):
            col=np.array([float(G[i].subs(list(zip(x,T[j,:])))) for j in range(0,ncons)]) #this multiplies lambda_i
            if np.max(col)<0:
                print('Warning: gamble '+str(i)+' is always negative in T')
            nonzeroind = np.nonzero(col)[0] 
            if len(nonzeroind)==0:
                print('Warning: constraint '+str(i)+' is zeros for all discretised values')
            else:
                #print([np.ndarray.tolist(ind_cons[nonzeroind]),np.ndarray.tolist(col[nonzeroind])])
                A.append([np.ndarray.tolist(ind_cons[nonzeroind]),np.ndarray.tolist(col[nonzeroind])])
                    
        #b vector    
        b=[float(hval[j]*f.subs(list(zip(x,T[j,:])))+epsilon) for j in range(0,ncons)]        
               
        #c vector
        c=np.zeros(nvar)    
        c[0]=-1.0 #we minimize -lambda_0
        #bounds
        if len(self.lambda_bounds)>0:
            bnd=np.array(self.lambda_bounds.copy())
            bnd=np.vstack([[f_range[0],f_range[1]],bnd])                
        else:
            bnd=np.vstack([[f_range[0],f_range[1]]])
        bnd[:,0]=[-cplex.infinity if v is None else v for v in bnd[:,0]]
        bnd[:,1]=[ cplex.infinity if v is None else v for v in bnd[:,1]]
        Lower=bnd[:,0]
        Upper=bnd[:,1]
        #Lower=np.hstack([f_range[0], np.zeros(nvar-1)])
        #Upper=[cplex.infinity]*nvar
        #Upper[0]=f_range[1]      
        
        prob.objective.set_sense(prob.objective.sense.minimize)
        #print(b)
        prob.linear_constraints.add(rhs = b, senses = ["L"]*(ncons))    
        prob.variables.add(obj = c, lb=Lower,  ub = Upper, columns=A)
    
    def _cplex(self,x,f,f_range,G,T,h=sympy.Piecewise((1.0,True)),verbose=True,epsilon=1e-8):
        import cplex as cplex
        my_prob = cplex.Cplex()
        handle=self.populatebycolumn(x,f,f_range,G,T,my_prob,h=h,epsilon=epsilon)
        my_prob.solve()
        if verbose:
            print(my_prob.solution.status[my_prob.solution.get_status()])
            print(my_prob.solution.get_objective_value())
        sol={'x':my_prob.solution.get_values(),
             'fun':-my_prob.solution.get_objective_value(),
             'LP max constraint violation': "To BE DONE for CPLEX"}
        return sol

    

    def _linprog(self,x,f,f_range, G,T,h=sympy.Piecewise((1.0,True)),verbose=True,epsilon=1e-8):
        #print(h)
        #this is the weighting function (partition of unity)        
        hval=np.array([[float(h.subs(list(zip(x,T[j,:])))) for j in range(0,T.shape[0])]]).T
        #print(hval)
        #define matrix A of the linear constraints
        A=np.multiply(np.ones((T.shape[0],1)),hval)#this multiplies lambda_0
        for i in range(0,len(G)):
            col=[float(G[i].subs(list(zip(x,T[j,:])))) for j in range(0,T.shape[0])] #this multiplies lambda_i
            #print(col)
            if np.max(col)<0:
                print('Warning: gamble '+str(i)+' is always negative in T')
            A=np.hstack([A,np.array([col]).T])
        
        #b vector 
        
        col=[float(hval[j]*f.subs(list(zip(x,T[j,:])))) for j in range(0,T.shape[0])]
        b=np.array([col]).T+epsilon
        
        #c vector
        c=np.zeros(A.shape[1])    
        c[0]=-1#we minimize -lambda_0
        
        #bounds
        if len(self.lambda_bounds)>0:
            bnd=np.array(self.lambda_bounds.copy())
            bnd=np.vstack([[f_range[0],f_range[1]],bnd])
        else:
            bnd=np.vstack([[f_range[0],f_range[1]]])
        #Lower=np.hstack([f_range[0], np.zeros(A.shape[1]-1)])
        #Upper=[None]*A.shape[1]
        #Upper[0]=f_range[1];
        #bnd=tuple(zip(Lower, Upper))

        res=linprog(c, A_ub=A, b_ub=b, A_eq=None, b_eq=None, bounds=bnd, method='interior-point')
        xopt=np.array([res.x]).T
        
        sol={'x':res.x,
             'fun':-res.fun,
             'LP max constraint violation': np.max(np.dot(A,xopt)-b)}
        
        if verbose:
            print('LP max constraint violation:',np.max(np.dot(A,xopt)-b))
        return sol

        

    def Cutting_plane(self,x,f,f_range,G,T,h=sympy.Piecewise((1.0,True)),verbose=True,solver='linprog',epsilonLP=1e-8,epsilonNLP=1e-3):
        loop = True
        eps=epsilonNLP
        niter=0
        
        while loop==True: 
            niter=niter+1

            resLP=self.LP(x,f,f_range,G,T,h=h,verbose=False,solver=solver,epsilon=epsilonLP)
            
            if verbose:
                print('LP iteration '+str(niter)+' fun=', resLP['fun'], ' number of support points='+str(T.shape[0]))
                #print('LP status:'+statusLP[resLP.status])
            lam = resLP['x'] 
            condition=(f-lam[0]-np.sum([lam[i+1]*G[i] for i in range(0,len(G))])) #the constraint we want to violate
            #D_condition=sp.diff(condition) #its jacobian

            def myfun(_z):
                _z=np.ndarray.tolist(_z)
                if self.is_mixed==True:
                    minv=np.inf
                    for d in self.DiscretePowerSet:
                        val=_z+d
                        cond=condition.subs(list(zip(x,val))) 
                        if cond<minv:
                            minv=cond
                            min_point=np.array(val)                    
                else:
                    val=_z
                    minv=condition.subs(list(zip(x,val)))
                    min_point=np.array(val) 
                c=np.zeros(1)
                c[0]=minv
                return c,min_point
            
            def myfun0(_z):
                c,none=myfun(_z)
                return c
    
            bnd=[]
            for ib in range(0,len(self.Symbols)):                
                if self._TypeSymbols[ib]=='c':
                    bnd.append([self._Domains[ib].boundary.inf,self._Domains[ib].boundary.sup]) 
            #bnd=np.array(bnd)
            #_z0 =self.make_SupportingPoints(1, criterion='random')

            resNP=differential_evolution(myfun0,bnd,seed=T.shape[0])

            min_fun,min_point=myfun(resNP.x)
            if verbose: 
                print('Violation of positivity of the constraint: ',min_fun)
            
            if min_fun<=-eps:
                T=np.vstack([T,min_point])
                if verbose: 
                    print('New support point:',resNP.x)
            else:    
                if verbose: 
                    print('Violation of positivity of the constraint is below the threshold ',eps)
                    loop=False      
    
        
        return resLP,T
    