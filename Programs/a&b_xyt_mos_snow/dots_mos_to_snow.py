import numpy as np
import math
import random
import matplotlib.pyplot as plt
import collections
import numpy as geek

# #_____________________________________

H=500

x=np.arange(0,500,10)
y=[0]*len(x)

#POINT ZERO x=0
Rr=948 #радиус кривизны зеркала  #  mm
rr=528 #радиус кривизны мозаики  #  mm
#R_mir=630 
#R_mir=Rr*math.sin(52.29)  #  mm
R_mir=750  #mm
#R_mos=rr*math.sin(28.36)  #  mm
R_mos=250   #mm

random.seed(11)    # seed the uniform random generator

# zmoscc=8 #мм
z0=-H*1000  #  mm #Высота 
#z1=rr*math.cos(28.36)
z1=rr*math.cos(0.495)
#z2=Rr*math.cos(52.29)
z2=Rr*math.cos(0.91267)

rhood=465  #  mm ok
zhood=-300  #  mm ok
tableau_x1=[]
tableau_y1=[]
w=[]
# for i in range(101):
  # print("")
  # print(" **** i=",i)
for i in range(0,len(x)):
    x01=x[i]*1000  #  mm
    y01=y[i]*1000  #  mm
  # print("  x0=",x0," y0=",y0)
    for j in range(1000):
        #  for j in range(1):
        # rop = rhood*math.sqrt(1)        !!!!
        rop = rhood*math.sqrt(random.random())  #  generate random radius-vector length
        xi = random.random()*2.0*np.pi          #  generate random azimuthal angle
        #    print("  rop=",rop," xi=",xi)
        # xop = rop*math.sqrt(math.cos(2.0*np.pi*2))     !!!
        xop = rop*math.cos(xi)
        # yop = rop*math.sin(2.0*np.pi*2)          !!!!
        yop = rop*math.sin(xi)
        #    xop=18.4172466455
        #    yop=423.55360337886
        #    print("    xop=",xop," yop=",yop)
        #ph0 =0
        th01=math.atan(math.sqrt((math.pow((xop-x01),2)+math.pow(yop-y01,2)))/(zhood-z0))
        if xop>x01:
          ph01=math.atan((yop-y01)/(xop-x01))
        if xop<x01:
          ph01=np.pi + math.atan((yop-y01)/(xop-x01))
        if xop==x01 and yop>=y01:
          ph01 = 0.5*np.pi
        if xop==x01 and yop<y01:
          ph01=1.5*np.pi
        if ph01<0:
          ph01=ph01+2.0*np.pi
        #    print("  th0=",th0," ph0=",ph0)
        #      break    
        # incident ray at mirror level:  (x2,y2,z2)
        t11 =(z2-z0)/math.cos(th01)
        x11 = x01 + t11 *math.sin(th01)*math.cos(ph01)
        y11 = y01 + t11*math.sin(th01)*math.sin(ph01)
        #    print("    x2=",x2," y2=",y2)
        # проверить, попадает ли точка (x2, y2) в круг зеркала радиуса
        if x11*x11+y11*y11<=R_mos*R_mos:
          continue
        #      break    
        # incident ray at mirror level:  (x2,y2,z2)
        t21 =(z2-z0)/math.cos(th01)
        x21 = x01 + t11 *math.sin(th01)*math.cos(ph01)
        y21 = y01 + t11*math.sin(th01)*math.sin(ph01)
        if x21*x21+y21*y21>R_mir*R_mir:
          continue
        # Определение точки P (точка на зеркале-> пересечение луча с зеркалом)
        A1 = x01*math.sin(th01)*math.cos(ph01) + y01*math.sin(th01)*math.sin(ph01) + z0*math.cos(th01)
        B1=x01*x01+y01*y01+z0*z0-Rr*Rr
        D1=A1*A1-B1
        #    print("     A=",A," B=",B,"  D=",D)
        if D1<0:
          continue
        #       break
        tP1=math.sqrt(D1)-A1
        xP1= x01 + tP1 *math.sin(th01)*math.cos(ph01)
        yP1= y01 + tP1*math.sin(th01)*math.sin(ph01)
        zP1 = z0 + tP1*math.cos(th01)
        #    print("     tP=",tP," xP=",xP," yP=",yP," zP=",zP)
        #Определение отраженного луча
        tA1 =(0.9*Rr*Rr -xP1*x01-yP1*y01-zP1*z0)/(xP1 * math.sin(th01)*math.cos(ph01)+ yP1*math.sin(th01)*math.sin(ph01) + zP1*math.cos(th01))
        xA1= x01 + tA1 * math.sin(th01)*math.cos(ph01)
        yA1 = y01 + tA1*math.sin(th01)*math.sin(ph01)
        zA1 = z0 + tA1*math.cos(th01)
        #    print("    xA=",xA," yA=",yA," zA=",zA)
        xC1=0.9*xP1
        yC1=0.9*yP1
        zC1=0.9*zP1
        #    print("    xC=",xC," yC=",yC," zC=",zC)
        xB1 = 2*xC1-xA1
        yB1 = 2*yC1-yA1
        zB1 = 2*zC1-zA1
        #    print("    xB=",xB," yB=",yB," zC=",zB)
        # reflected ray at mosaic level:  (x1,y1,z1)
        t11=(z1-zP1)/(zB1-zP1)
        x11 = xP1 + t11*(xB1-xP1)
        y11 = yP1 + t11*(yB1-yP1)
        #    print("    t1=",t1," x1=",x1," y1=",y1)
        if x11*x11+y11*y11>R_mos*R_mos: #неравенство гарантирующее прохождение луча мимо мозаики
          continue
        S1=math.pow((xB1-xP1),2) + math.pow((yB1-yP1),2)+math.pow((zB1-zP1),2)
        W1= xP1*xB1 + yP1*yB1 + zP1*zB1-Rr*Rr
        # U = R_mir*R_mir-r*r                        !!!  r->R_mos ,  see Fig.2 in image_formation.pdf
        U = Rr*Rr-rr*rr
    
        #    print("   U=",U," W=",W," S=",S)
        tauQ1 =(-W1-math.sqrt(W1*W1-S1*U))/S1
        xQ1 = xP1 + tauQ1*(xB1-xP1)
        yQ1 = yP1 + tauQ1 *(yB1-yP1)
        zQ1 = zP1 + tauQ1 * (zB1-zP1)
        # print("  >>>>>  point Q found!  i=",i," j=",j," tauQ=",tauQ," xQ=",xQ," yQ=",yQ," zQ=",zQ)
        # print(" ")
        tableau_x1.append(xQ1)
        tableau_y1.append(yQ1)
        w.append(i)
        x_y1=np.array([tableau_x1,tableau_y1,w])
        x_y1=x_y1.T

m=collections.Counter(w) 
dict = m 
result = dict.items() 
data = list(result) 
tt = np.array(data) 
t = geek.cumsum(tt[:,1])

t=np.insert(t,0,0)

decoupage_x=[tableau_x1[t[i]:t[i+1]] for i in range(0,len(t)-1)]
decoupage_y=[tableau_y1[t[i]:t[i+1]] for i in range(0,len(t)-1)]

moy_x=[]
moy_y=[]
for i in range(0,len(decoupage_x)):
    moy_x.append(np.mean(decoupage_x[i]))
    moy_y.append(np.mean(decoupage_y[i]))
    
d_mos=[np.sqrt(moy_x[i]**2+moy_y[i]**2) for i in range(0,len(moy_x))] #distance du centre des photons sur la mosaique 
d_snow=x

a=[((d_mos[i]*H*1000)/422)/1000 for i in range(0,len(d_mos))] #метод сферического зеркала

#Approximation x snow -> x mos
cut=len(moy_x)-5
a,b,c=np.polyfit(-np.array(moy_x)[:cut],x[:cut],2)
approx_onlyx_X=np.arange(0,max(-np.array(moy_x)[:cut]),1)
approx_onlyx_Y=a*approx_onlyx_X**2+b*approx_onlyx_X



#Graphiques      
params = {"xtick.direction": "in", "ytick.direction": "in"}
plt.rcParams.update(params)

# plt.scatter(x,y,s=1)
# plt.show()
# plt.scatter(tableau_x1,tableau_y1,s=1)
# plt.xlim(-400,400)
# plt.ylim(-400,400)
# plt.show()

# plt.scatter(d_mos,d_snow,label='по данным телескопа')
# plt.scatter(d_mos,a, label='по формуле сферического зеркала')
# plt.xlabel('r, мм')
# plt.ylabel('R, м')
# plt.title('Зависимость точек на снегу от точек на мозаике')
# plt.legend()
# plt.xlim(0,250)
# plt.ylim(0,510)
# plt.grid(True)
# plt.show()

plt.scatter(-np.array(moy_x)[:cut]/10,x[:cut])
plt.plot(approx_onlyx_X/10,approx_onlyx_Y,c='red',label='a={:.6f}, b={:.3f}'.format(a/10,b/10))
plt.ylabel('R, м')
plt.xlabel('r, см')
plt.legend(title='H={} м \n R=ar^2+br'.format(H))
plt.grid(True)
plt.xlim(0,32)
plt.ylim(0,260)
plt.show()

# #r(R)
# plt.scatter(-np.array(moy_x)[:cut],x[:cut])
# plt.plot(approx_onlyx_X,approx_onlyx_Y,c='red',label='a={:.5f}, b={:.2f}'.format(a,b))
# plt.ylabel('R, m')
# plt.xlabel('r, mm')
# plt.legend(title='H={} m'.format(H))
# plt.grid(True)
# plt.xlim(0,None)
# plt.ylim(0,None)
# plt.show()



