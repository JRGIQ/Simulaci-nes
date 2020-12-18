#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 18:00:10 2018

@author: jheison
"""

from control import*
from control.matlab import*
import math as mt
from scipy import* 
import matplotlib.pyplot as plt
g=tf(1,[1,1,1])
print(g)
t,y=step_response(g)
plt.plot(t,y)
plt.grid()
plt.xlabel('tiempo')
plt.ylabel('respuesta escalon')