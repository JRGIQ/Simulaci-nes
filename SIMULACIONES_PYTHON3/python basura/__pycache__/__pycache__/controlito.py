#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 18:23:57 2018

@author: jheison
"""


from control import*
import math as mt
import sympy as sy
import control.matlab as cm
import matplotlib.pyplot as plt
g=tf(1,[1,1,1])
print(g)
t,y=impulse_response(g) 
plt.plot(t,y)
t,y=step_response(g)
plt.plot(t,y)
plt.grid()
plt.xlabel('tiempo')
plt.ylabel('impulse y step')
####################
x=mt.exp(1)
f=cm.pade(x,2)
#fg=tf([-1.0, 0.7357588823428847], [1.0, 0.7357588823428847])
print(f)


