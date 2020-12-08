
import numpy as np
from tqdm import tqdm
from pylab import plot, show, xlabel, ylabel, legend, title
import math

#modifying functions
def alpha_n(V):
	return 0.01*(V-10)/(1-math.exp(-(V-10)/10.0))

def beta_n(V):
	return 0.125*math.exp(-0.0125*(V))

def alpha_m(V):
	return 0.1*(V-25)/(1-math.exp(-(V-25)/10.0))

def beta_m(V):
	return 4*math.exp(-0.0556*(V))

def beta_h(V):
	return 1.0/(1+(math.exp(-(V-30)/10.0)))

def alpha_h(V):
	return 0.07*math.exp(-0.05*(V))

def del_n(an, bn, n):
	#print("del_n: ", alpha_n*(1-n) - beta_n*n)
	return an*(1-n) - bn*n

def del_m(am, bm, m):
	#print(alpha_m*(1-m) - beta_m*m)
	return am*(1-m) - bm*m

def del_h(ah, bh, h):
	#print("del_h ", alpha_h*(1-h) - beta_h*h)
	return ah*(1-h) - bh*h

#constants
Cm = 1.0 #nF for 0.1 mm^2
gNa = 1.2 #mS/mm^2
gK = 0.36 #mS/mm^2
gL = 0.003 #mS/mm^2
Vk = -77.0 #mV
VNa = 50.0 #mV
VL = -54.387 #mV

#initial conditions
m0 = 0.0
h0 = 1.0
n0 = 0.0
I = 100.0
V0 = -60.0 #V initial equi potential

#Calculating current
def del_Ik(n, m, h, V):
	K = gK*n**4*(V-Vk)
	Na = gNa*m**3*h*(V-VNa)
	L = gL*((V-VL))
	return (K+Na+L)

#Voltage updation
def del_v(Ik, I):
	return (I-Ik)/Cm

#Perfrming integration for time = b seconds
N = 100000
b = 200
time = np.linspace(0, b, N)
delh = 10.0/N
V = np.linspace(-100,100, N)
delv = float(b)/N
Voltages = []
current = []
n_s = []
m_s = []
h_s = []
an_s=[]
bn_s=[]
xn_s=[]
xh_s=[]
xm_s=[]

for t in tqdm(V):
    	
	# if int(int(t))%250==0 and t>2 :
	# 	I=200.0
	# else:
	# 	I=0.0
	# if t>100 and t<101:
	# 	I=200.0
	# else:
	# 	I=0.0
	n_s.append(n0)
	m_s.append(m0)
	h_s.append(h0)
	Voltages.append(V0)

	n_alpha = alpha_n(V0)
	h_alpha = alpha_h(V0)
	m_alpha = alpha_m(V0)
	#an_s.append(n_alpha)
	n_beta = beta_n(V0)
	h_beta = beta_h(V0)
	m_beta = beta_m(V0)
	#bn_s.append(n_beta)

	xn_s.append(1.0/(n_alpha+n_beta)*n_alpha)
	xm_s.append(1.0/(m_alpha+m_beta)*m_alpha)
	xh_s.append(1.0/(h_alpha+h_beta)*h_alpha)
	
	n1 = n0 + delv*del_n(n_alpha, n_beta, n0)
	n0 = n1
	m1 = m0 + delv*del_m(m_alpha, m_beta, m0)
	m0 = m1
	h1 = h0 + delv*del_h(h_alpha, h_beta, h0)
	h0 = h1
	#current.append(del_Ik(n0,m0,h0,V0))
	V1 = V0 + delv*del_v(del_Ik(n0,m0,h0,V0),I)
	V0 = V1


plot(Voltages, xn_s, label="n Infinity")
plot(Voltages, xm_s, label="m Infinity")
plot(Voltages, xh_s, label="h Infinity")
#plot(time[50:], Voltages[50:], label="Voltages")
#plot(time, m_s, label="m")
#plot(time, n_s, label="n", color="green")
#plot(time, h_s, label="h", color="black")
xlabel("Voltage")
ylabel('Probabilities')
legend()
title("Limiting Values for varying Voltages")
show()
