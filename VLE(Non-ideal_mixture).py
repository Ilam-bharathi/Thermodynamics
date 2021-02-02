import numpy as np
import matplotlib.pyplot as plt
#given
x1=0.581 #x_EA
xea=x1
x2=1-x1
paz=1.165
p=paz*np.ones(2)
T=(273.15+76)*np.ones(2)

def coeff(p,pi,x1,x2):#Finding alpha and beta
    g=p/pi #since x=y at the azeotrope in raoult's law
    alp=((1+((x2)*np.log(g[1]))/((x1)*np.log(g[0])))**2)*np.log(g[0])
    bet=((1+((x1)*np.log(g[0]))/((x2)*np.log(g[1])))**2)*np.log(g[1])
    return(alp,bet) #rewritten equations
#1    
A=[9.6830,9.3171]
B=[2842.2,2810.5]
C=[56.3209,51.2586]
pi=np.exp(A-(B/(T-C)))
alp,bet=coeff(p,pi,x1,x2) #Antoine equation
g=p/pi
print("The constants are found to be alpha =",alp,", beta=",bet)

#2
x1=np.arange(0.001,1,0.01) #Liquid mole fraction in limits of 0 to 1, in steps of 0.1
x1[-1]=0.9999
x2=1-x1
g1=np.exp(alp*(1+(alp*x1)/(bet*x2))**-2) #given van-laar
g2=np.exp(bet*(1+(bet*x2)/(alp*x1))**-2)
p=pi[0]*g1*x1+pi[1]*g2*x2 #Raoult's law
y1=pi[0]*g1*x1/p #generalized raoult's
plt.figure(figsize=(8,3.5),dpi=100)
plt.plot(x1,p,'k')
plt.plot(xea,paz,'ro')
plt.ylabel(r'$P$'+ " (bar)")
plt.xlabel(r'$X_{EA}$'+' (Liquid mole fraction)')
plt.text(0.1, 1.04, "Liquid", ha='left', rotation=45)
plt.legend(['Bubble curve at 76'+r'$^{\circ}$'+'C','Azeotrope'],loc='lower right',fontsize='small')
plt.title(r'$P$'+' vs '+r'$X_{EA}$'+' (Ethyl Acetate)')
plt.show()

#combined
plt.figure(figsize=(8,3.5),dpi=100)
plt.plot(x1,p,'k')
plt.plot(y1,p,'b')
plt.plot(xea,paz,'ro')
plt.ylabel(r'$P$'+ " (bar)")
plt.xlabel(r'$X_{EA}$'+' (Liquid mole fraction)')
plt.legend(['Bubble curve at 76'+r'$^{\circ}$'+'C','Dew curve at 76'+r'$^{\circ}$'+'C','Azeotrope'],loc='lower right',fontsize='small')
plt.title('Phase curves at 76'+r'$^{\circ}$'+'C')
plt.text(0.1, 1.04, "Liquid", ha='left', rotation=45)
plt.text(0.2, 0.94, "Vapor", ha='left', rotation=45)
plt.text(0.1, 0.96, "Co-existence region", ha='left', rotation=45)
plt.show()


#3
plt.figure(figsize=(8,3.5),dpi=100)
y1=pi[0]*g1*x1/p
y2=pi[1]*g2*x2/p
plt.plot(x1,y1,'k')
plt.plot(x1,x1,'g--')
plt.plot(0.581,0.581,'ro')
plt.ylabel(r'$Y_{EA}$'+' (Vapor mole fraction)')
plt.xlabel(r'$X_{EA}$'+' (Liquid mole fraction)')
plt.legend(['X vs Y',"X=Y",'Azeotrope'],loc='lower right',fontsize='small')
plt.title(r'$Y$'+' vs '+r'$X$'+ ' at 76'+r'$^{\circ}$'+'C'+' (Ethyl acetate)')
plt.show()

#4
plt.figure(figsize=(8,3.5),dpi=100)
p1=1.013
#Assume a temperature
T=(273.15+70)*np.ones(2)
pa=[]
pb=[]
To=[]
ka=1
w=0
for i in range(len(x1)):
    while ka >10**-3:
        pi=np.exp(A-(B/(T-C)))
        p=pi[0]*g1[i]*x1[i]+pi[1]*g2[i]*x2[i]
        sumy=p/p1
        k=sumy-1
        ka=abs(k)
        w=T[0]
        if k < 0:
            T=T+0.0001
        if k > 0:
            T=T-0.0001
    ka=1
    pa=np.append(pa,pi[0])
    pb=np.append(pb,pi[1])
    To=np.append(To,w)
plt.plot(x1,To-273.15,color='black')
q=np.sort(abs(ya-x1))
zx=np.array(x1)
c=zx[abs(ya-x1)==q[0]]
print("The Azeotrope for 1.013 bar is: ",c)
plt.plot(c,To[abs(ya-x1)==q[0]]-273.15,'ro')
plt.ylabel(r'$T$'+r'$^{\circ}$'+'C')
plt.text(0.1, 74, "Liquid", ha='left', rotation=-45)
plt.xlabel(r'$X_{EA}$'+' (Liquid mole fraction)')
plt.legend(['Bubble curve at 1.013 bar','Azeotrope'],loc='upper right',fontsize='small')
plt.title(r'$T$'+' vs '+r'$X_{EA}$'+' (Ethyl Acetate)')
plt.show() 

#5
plt.figure(figsize=(8,3.5),dpi=100)
ya=pa*g1*x1/p1   
plt.plot(x1,ya,color='black') 
plt.plot(x1,x1,'g--')
plt.plot(c,c,'ro')
plt.ylabel(r'$Y_{EA}$'+' (Vapor mole fraction)')
plt.xlabel(r'$X_{EA}$'+' (Liquid mole fraction)')
plt.legend(['X vs Y ','X=Y','Azeotrope'],loc='lower right',fontsize='small')
plt.title(r'$Y$'+' vs '+r'$X$'+' at 1.013 bar (Ethyl acetate)' )
plt.show() 

#combined
plt.figure(figsize=(8,3.5),dpi=100)
plt.plot(x1,To-273.15,color='r') 
plt.plot(y1,To-273.15,'b')
plt.plot(c,To[x1==c]-273.15,'ro')
plt.text(0.1, 74, "Co-existence region", ha='left', rotation=-45)
plt.text(0.2, 76, "Vapor", ha='left', rotation=-45)
plt.text(0.1, 74, "Liquid", ha='left', rotation=-45)
plt.ylabel(r'$T$'+r'$^{\circ}$'+'C')
plt.xlabel(r'$X_{EA}$'+' (Liquid mole fraction)')
plt.legend(['Bubble curve at  1.013 bar','Dew curve at 1.013 bar','Azeotrope'],loc='upper right',fontsize='small')
plt.title('Phase curves at 1.013 bar')
plt.show() 
