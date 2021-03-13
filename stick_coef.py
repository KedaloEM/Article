import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

nu1 = 10 ** (11)
nu_des = 10 ** (14)
nu_r = 10 ** (13)
kb = 1.38 * 10 ** (-23)
Ea = 1.6 * 1.6 * 10 ** (-19)
Edes = 1 * 1.6 * 10 ** (-19)
Ea_step = 0.8 * 1.6 * 10 ** (-19)
Edes_step = 0.9 * 1.6 * 10 ** (-19)
E_inv1 = 1.5*1.6*10**(-19)
E_inv2 = 1.6*10**(-19)
alpha = 0.5
NL = 10 ** (7)
a_max = 7.5 * 10 ** (18)
k = -0.0064
b = 11.41
P = 0.57*10**2



def flux(T, t):
    Ng = P / (T)
    m = 1.6 * 10 ** (-27) * 28
    V = (3 * T/ m) ** (1 / 2)
    s = 0.3
    flux = Ng * V * s * t / 4
    return flux


def fun(y, t, T,step_rate):
    z = -np.exp(-Edes / T) * y[0] - nu_r / nu_des * np.exp(-Ea / T) * y[0] * (1 - y[0]/(a_max*(1-step_rate)) - y[1]/(a_max*(1-step_rate))) + (1-step_rate)*flux(T, 1) / (nu_des) * (1 - y[0]/(a_max*(1-step_rate)) - y[1]/(a_max*(1-step_rate))) + y[1]**2*10**(-6)*np.exp(-E_inv1/T)/nu_des
    x = 2 * nu_r / nu_des * np.exp(-Ea / T) * y[0] * (1 - y[0]/(a_max*(1-step_rate)) - y[1]/(a_max*(1-step_rate))) - 2*y[1]**2*10**(-6)*np.exp(-E_inv1/T)/nu_des
    if step_rate > 0:
        q = -np.exp(-Edes_step / T) * y[2] - nu_r / nu_des * np.exp(-Ea_step / T) * y[2] * (1 - y[2]/(a_max*step_rate) - y[3]/(a_max*step_rate)) + step_rate*flux(T, 1) / (nu_des) * (1 - y[2]/(a_max*step_rate) - y[3]/(a_max*step_rate)) + y[3]**2*10**(-6)*np.exp(-E_inv1/T)/nu_des
        r = 2 * nu_r / nu_des * np.exp(-Ea_step / T) * y[2] * (1 - y[2]/(a_max*step_rate) - y[3]/(a_max*step_rate)) - 2*y[3]**2*10**(-6)*np.exp(-E_inv1/T)/nu_des
    else:
        q = 0
        r = 0
    return [z, x, q, r]


def solver(T, step_rate,total_time):
    time1 = np.arange(0, total_time, total_time/NL)
    return integrate.odeint(fun, [0, 0, 0, 0], time1, args = (T,step_rate,))


T_range = np.arange(300, 1200, 50)
sticks = np.zeros(len(T_range))
step_rate = np.array([0,0.001,0.01,0.05,0.1])
fig,ax = plt.subplots()
for j in range(len(step_rate)):
    for i in range(len(T_range)):
        total_time = nu_des*10**((k*T_range[i]+b))*1.33*10**(-6)/(P*10**(-2))
        #if i > (4/5)*len(T_range):
        #    P = 6.3*10**(-5)
        ui = solver(T_range[i]*kb,step_rate[j],total_time)
        sticks[i] = ((1 - step_rate[j]) * ui[NL - 1, 1] + step_rate[j] * ui[NL - 1, 3])*nu_des/flux(T_range[i]*kb,total_time)
    line, = ax.plot(T_range,np.log10(sticks))
    if step_rate[j] == 0:
        line.set_label('terrace')
    else:
        line.set_label(r'$\theta = %.03f$'% step_rate[j])

x_exp = np.array([300,350,400,424,448])
y_exp = np.array([-11.6,-10.8,-10.1,-9.8,-9.6])
x1_exp = np.array([1160])
y1_exp = np.array([-6.2])
ax.scatter(x_exp,y_exp,s = 120, c = 'black',marker = 'o')
ax.scatter(x1_exp,y1_exp,s = 120, c = 'black',marker = 'x')
ax.legend(loc = 4, fontsize = 18)
ax.set_xlabel('T(K)', fontsize = 20)
ax.set_ylabel('log(S)', fontsize = 20)
ax.tick_params(labelsize=18)
plt.show()
