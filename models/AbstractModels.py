#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main Gamble Module

Defines abstract BeliefModel

Created on Thu Aug 31 08:34:59 2018

@author: benavoli
"""
from abc import ABCMeta, abstractmethod


class BeliefModel(object, metaclass=ABCMeta):
    #__metaclass__  = abc.ABCMeta
    
    @abstractmethod
    def buildModel(self):
         """Builds a BeliefModel"""
    
    @abstractmethod
    def coherence(self):
         """Checks coherence (consistency) of the BeliefModel"""
    
    @abstractmethod
    def inference(self,f):
         """Returns the result of an inference on a BeliefModel"""
         
    @abstractmethod
    def updating(self,h):
         """Updates a BeliefModel"""
            
    @abstractmethod
    def revision(self,h):
         """Revises a BeliefModel"""