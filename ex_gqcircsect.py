#compute the integral of f in a circular sector for the interval (-omega,omega)
#prima area
n=5
omega=np.pi/6
f=lambda x,y: x-y**3+x**7*y
xyw=gqcircsect(n,omega,1,3)
ris=sum(xyw[:, 2]*f(xyw[:, 0],xyw[:, 1]))
