#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 15:29:50 2018

@author: jheison
""






from control.matlab import*
import matplotlib.pyplot as plt
g=tf(1,[1,1,1])
print(g)
t,y=step_response(g)
plt.plot(t,y)
plt.grid()
plt.xlabel('tiempo')
plt.ylabel('respuesta escalon')