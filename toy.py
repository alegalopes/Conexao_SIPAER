# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 08:14:55 2018

@author: alexandre lopes
"""
import numpy as np
import matplotlib.pyplot as plt
#%% Dados do problema
vc=117
ve=38
m=11610
p0=101325
p=101000
T0=288.15
r0=1.225
mu0=17.894*(10**-6)
gt=-0.0065
T=T0+18
g=9.81
R=287
S=42.4
b=22.81
c=S/b
e=0.8
e1=0.85
RA=(b**2)/S
W=m*g
#%% Funções de conversão
#Equação de estado
def rho(p,T):
    ro=p/(R*T)
    return(ro)
#Equação de Sutherland
def mu(T):
    mi=mu0*((T/T0)**(3/2))*((T0+110)/(T+110))
    return(mi)
#Variáveis corrigidas
r=rho(p,T)
mu=mu(T)
#%% Função para obtenção do número de Reynolds
def re(r,v,d,mi):
    rn=r*v*d/mi
    return(rn)
#Nr Reynolds para velocidades de cruzeiro e de estol
rec=re(r,vc,c,mu)
ree=re(r,ve,c,mu)
#Obtenção da velocidade para o Reynolds escolhido
rei=9*(10**6) #Reynolds ideal
reref=0
v=ve
while (rei-reref)>=(10**3):
    v=v+1
    reref=re(r,v,c,mu)
v=v-1
#%% Correção da inclinação devida à razão de aspecto elevada
#Ângulos em graus
alp1=-1.8
cl1=0
alp2=11
cl2=1.3
a0=(cl2-cl1)/(alp2-alp1)
a=a0/(1+(57.3*a0/(np.pi*e1*RA)))
#cl para Alpha=0, para asa infinita
cli=0.18
#cl para Alpha=0, para asa finita
clf=-a*alp1+cl1
#Geração do gráfico
ax = plt.subplot(111)
degx = np.arange(alp1, alp2, 0.01)
cly = a0*degx+cli
cl0y = a*degx+clf
plt.plot(degx,cly,label="a0=%.3f"%(a0,))
plt.plot(degx,cl0y,label="a=%.3f"%(a,))
leg = plt.legend(loc='best', ncol=2, mode="expand", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.grid(1)
plt.title('Aerofólio NACA 63-210 - Região Linear')
plt.xlabel('Ângulo de ataque, em graus')
plt.ylabel('Coeficiente de sustentação')
plt.show()
#%% Cálculo do CL
CL=W/(0.5*r*(v**2)*S)
#%% Cálculo do novo CL, com variação de 5°
alf=5.8+5
CLn=a0*(alf-alp1)
#%% Nova sustentação
L=CLn*0.5*r*(v**2)*S
av=(L-W)/m
ng=av/g
#%% Cálculo do CD
cd0=0.005
CD=cd0+(CLn**2)/(np.pi*e*RA)
D=CD*0.5*r*(v**2)*S
#%% Cálculo da tração
Tr=W/(CL/0.0075)
#%% Cálculo da aceleração horizontal
ah=(Tr-D)/m
#%% Aceleeração resultante
at=np.sqrt(av**2+ah**2)
