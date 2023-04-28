import numpy as np 
import matplotlib.pyplot as plt

Cr_025 = np.loadtxt('correlation_0.25.txt')
Cr_05 = np.loadtxt('correlation_0.5.txt')

# Cr_025 = np.round(Cr_025) #round up to avoid too much precision
# Cr_05 = np.round(Cr_05)

print(len(Cr_025))
print(len(Cr_05))
r = np.linspace(0, len(Cr_025) - 1, len(Cr_025))


def Cr_anal(L, r, T):
    J =  1
    beta = 1/T
    a = np.exp(beta*J)
    ledd1 = (a+2)**r * (a-1)**(L-r) + (a-1)**r * (a+2)**(L-r) + (a-1)**L
    ledd2 = 2 * (a-1)**L + (a+2)**L
    return ledd1/ledd2

Cr_anal_025 = Cr_anal(16, r, 0.25)
Cr_anal_05 = Cr_anal(16, r, 0.5)

fig = plt.figure()
plt.scatter(r, Cr_025, label = 'MC, T = 0.25 J')
plt.scatter(r, Cr_05, label = 'MC, T = 0.5 J')
plt.plot(r, Cr_anal_025, label = 'Analytical, T = 0.25 J')
plt.plot(r, Cr_anal_05, label = 'Analytical, T = 0.5 J')
plt.xlabel('r')
plt.ylabel('C(r)')
plt.legend()
plt.title('Correlation function for L = 16')
# plt.ylim(min(correlation), max(correlation))
plt.savefig('correlation_L16.pdf')
plt.show()

#plotting the real part of magnetisation as a func of temperature 

magnetic_moments_L16_T0 = np.loadtxt('magnetic_moments_T0.000000.txt')
magnetic_moments_L16_T025 = np.loadtxt('magnetic_moments_T0.250000.txt')
magnetic_moments_L16_T05 = np.loadtxt('magnetic_moments_T0.500000.txt')
magnetic_moments_L16_T075 = np.loadtxt('magnetic_moments_T0.750000.txt')
magnetic_moments_L16_T1 = np.loadtxt('magnetic_moments_T1.000000.txt')
magnetic_moments_L16_T125 = np.loadtxt('magnetic_moments_T1.250000.txt')
magnetic_moments_L16_T15 = np.loadtxt('magnetic_moments_T1.500000.txt')
magnetic_moments_L16_T175 = np.loadtxt('magnetic_moments_T1.750000.txt')
magnetic_moments_L16_T2 = np.loadtxt('magnetic_moments_T2.000000.txt')


#do the same for files with L=8
magnetic_moments_L8_T0 = np.loadtxt('magnetic_moments_L8_T0.000000.txt')
magnetic_moments_L8_T025 = np.loadtxt('magnetic_moments_L8_T0.250000.txt')
magnetic_moments_L8_T05 = np.loadtxt('magnetic_moments_L8_T0.500000.txt')
magnetic_moments_L8_T075 = np.loadtxt('magnetic_moments_L8_T0.750000.txt')
magnetic_moments_L8_T1 = np.loadtxt('magnetic_moments_L8_T1.000000.txt')
magnetic_moments_L8_T125 = np.loadtxt('magnetic_moments_L8_T1.250000.txt')
magnetic_moments_L8_T15 = np.loadtxt('magnetic_moments_L8_T1.500000.txt')
magentic_moments_L8_T175 = np.loadtxt('magnetic_moments_L8_T1.750000.txt')
magnetic_moments_L8_T2 = np.loadtxt('magnetic_moments_L8_T2.000000.txt')

#do the same for files with L=32 
magnetic_moments_L32_T0 = np.loadtxt('magnetic_moments_L32_T0.000000.txt')
magnetic_moments_L32_T025 = np.loadtxt('magnetic_moments_L32_T0.250000.txt')
magnetic_moments_L32_T05 = np.loadtxt('magnetic_moments_L32_T0.500000.txt')
magnetic_moments_L32_T075 = np.loadtxt('magnetic_moments_L32_T0.750000.txt')
magnetic_moments_L32_T1 = np.loadtxt('magnetic_moments_L32_T1.000000.txt')
magnetic_moments_L32_T125 = np.loadtxt('magnetic_moments_L32_T1.250000.txt')
magnetic_moments_L32_T15 = np.loadtxt('magnetic_moments_L32_T1.500000.txt') 
magnetic_moments_L32_T175 = np.loadtxt('magnetic_moments_L32_T1.750000.txt')
magnetic_moments_L32_T2 = np.loadtxt('magnetic_moments_L32_T2.000000.txt')

#extract the last element of the first column of the array - that is real part of the average magnetisation

#make a function to extract the real part of the average magnetisation for each temperature
def m_avg(array_name):
    magnetic_moments = array_name
    m_avg = magnetic_moments[-1, 0]
    return m_avg

list_of_m_avg_L16 = [m_avg(magnetic_moments_L16_T0), m_avg(magnetic_moments_L16_T025), m_avg(magnetic_moments_L16_T05), m_avg(magnetic_moments_L16_T075), m_avg(magnetic_moments_L16_T1), m_avg(magnetic_moments_L16_T125), m_avg(magnetic_moments_L16_T15), m_avg(magnetic_moments_L16_T175), m_avg(magnetic_moments_L16_T2)]
list_of_m_avg_L16 = np.round(list_of_m_avg_L16, 2)
list_of_T = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]

"""
fig = plt.figure()
plt.scatter(list_of_T, list_of_m_avg_L16)   
plt.plot(list_of_T, list_of_m_avg_L16)
plt.xlabel('T [J]')
plt.ylabel('<m>')
plt.title('Average magnetisation as a function of temperature for a 2D lattice with L=16')
fig.savefig('m_avg_T_2D.pdf')
plt.show()
"""


def m2_avg(array_name):
    magnetic_moments = array_name
    m2_avg = magnetic_moments[-1, 2]
    return m2_avg

list_of_m2_avg_L16 = [m2_avg(magnetic_moments_L16_T0), m2_avg(magnetic_moments_L16_T025), m2_avg(magnetic_moments_L16_T05), m2_avg(magnetic_moments_L16_T075), m2_avg(magnetic_moments_L16_T1), m2_avg(magnetic_moments_L16_T125), m2_avg(magnetic_moments_L16_T15), m2_avg(magnetic_moments_L16_T175), m2_avg(magnetic_moments_L16_T2)]
print(list_of_m2_avg_L16)


"""
fig = plt.figure()
plt.scatter(list_of_T, list_of_m2_avg_L16)
plt.xlabel('T [J]')
plt.ylabel('$<m^2>$')
plt.title('Average magnetisation squared as a function of temperature for a 2D lattice with L = 16 ')
plt.savefig('m2_avg_T_2D.pdf')
plt.show()
"""



#compute the forth moment divided by the second moment squared
def m4_avg(array_name):
    magnetic_moments = array_name
    m4_avg = magnetic_moments[-1, 3]
    return m4_avg

list_of_m4_avg_L16 = [m4_avg(magnetic_moments_L16_T0), m4_avg(magnetic_moments_L16_T025), m4_avg(magnetic_moments_L16_T05), m4_avg(magnetic_moments_L16_T075), m4_avg(magnetic_moments_L16_T1), m4_avg(magnetic_moments_L16_T125), m4_avg(magnetic_moments_L16_T15), m4_avg(magnetic_moments_L16_T175), m4_avg(magnetic_moments_L16_T2)]

#extract the 2nd and 4th magnetic moment for L = 8 system 
list_of_m2_avg_L8 = [m2_avg(magnetic_moments_L8_T0), m2_avg(magnetic_moments_L8_T025), m2_avg(magnetic_moments_L8_T05), m2_avg(magnetic_moments_L8_T075), m2_avg(magnetic_moments_L8_T1), m2_avg(magnetic_moments_L8_T125), m2_avg(magnetic_moments_L8_T15), m2_avg(magentic_moments_L8_T175), m2_avg(magnetic_moments_L8_T2)]
list_of_m4_avg_L8 = [m4_avg(magnetic_moments_L8_T0), m4_avg(magnetic_moments_L8_T025), m4_avg(magnetic_moments_L8_T05), m4_avg(magnetic_moments_L8_T075), m4_avg(magnetic_moments_L8_T1), m4_avg(magnetic_moments_L8_T125), m4_avg(magnetic_moments_L8_T15), m4_avg(magentic_moments_L8_T175), m4_avg(magnetic_moments_L8_T2)]

#do the same for L = 32 system

list_of_m2_avg_L32 = [m2_avg(magnetic_moments_L32_T0), m2_avg(magnetic_moments_L32_T025), m2_avg(magnetic_moments_L32_T05), m2_avg(magnetic_moments_L32_T075), m2_avg(magnetic_moments_L32_T1), m2_avg(magnetic_moments_L32_T125), m2_avg(magnetic_moments_L32_T15), m2_avg(magnetic_moments_L32_T175), m2_avg(magnetic_moments_L32_T2)]
list_of_m4_avg_L32 = [m4_avg(magnetic_moments_L32_T0), m4_avg(magnetic_moments_L32_T025), m4_avg(magnetic_moments_L32_T05), m4_avg(magnetic_moments_L32_T075), m4_avg(magnetic_moments_L32_T1), m4_avg(magnetic_moments_L32_T125), m4_avg(magnetic_moments_L32_T15), m4_avg(magnetic_moments_L32_T175), m4_avg(magnetic_moments_L32_T2)]

gamma_L16 = np.divide(list_of_m4_avg_L16, np.square(list_of_m2_avg_L16))
gamma_L8 = np.divide(list_of_m4_avg_L8, np.square(list_of_m2_avg_L8))
gamma_L32 = np.divide(list_of_m4_avg_L32, np.square(list_of_m2_avg_L32))

"""
fig = plt.figure()
plt.scatter(list_of_T, gamma_L16, label = 'L = 16', color = 'red')
plt.scatter(list_of_T, gamma_L8, label = 'L = 8', color = 'blue')
plt.scatter(list_of_T, gamma_L32, label = 'L = 32', color = 'green')
plt.plot(list_of_T, gamma_L16, color = 'red')
plt.plot(list_of_T, gamma_L8, color = 'blue')
plt.plot(list_of_T, gamma_L32, color = 'green')
plt.legend()
plt.title("Ratio of dimensionless magnetic moments $\Gamma$ for different system sizes")
plt.savefig('gamma_L.pdf' )
plt.show()
"""




