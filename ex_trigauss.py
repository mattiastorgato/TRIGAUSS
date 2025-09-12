#compute the integral of f in (a,b) 
#INTERVALLO 1

nv=5
a=np.pi/6
b=np.pi/4
f=lambda x: 5+(1/2)*np.sin(17*x)-6*np.cos(14*x)

tw=trigauss(nv,a,b)
int=sum(tw[:,1]*f(tw[:,0]))
