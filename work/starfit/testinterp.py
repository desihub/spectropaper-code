#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from desispec.fluxcalibration import _compute_coef

# nodes
x=np.tile(np.arange(4),(4,1))
y=x.copy().T
nodes=np.array([x.ravel(),y.ravel()]).T
print(nodes.shape)

# test points
cx=np.tile(np.linspace(0,np.max(x),30),(30,1))
cy=cx.copy().T
cx=cx.ravel()
cy=cy.ravel()

# compute iterpolation coefs
coefs=[]
for x,y in zip(cx,cy) :
    coefs.append( _compute_coef([x,y],nodes) )
coefs=np.array(coefs)

# plot nodes and some coeffs
nn=np.arange(nodes.shape[0],dtype=int)[::-1]
for n in nn :
    color=c="C{}".format(n%10)
    plt.plot(nodes[n,0],nodes[n,1],"X",color=color,markersize=12)
    if n==6 or n==1 :
        for i in np.where(coefs[:,n]>0)[0] :
            plt.plot(cx[i],cy[i],"o",color=color,alpha=coefs[i,n])
plt.show()
