import numpy as np
from scipy import optimize

'''
Modelo de potencia de Morrell
Relaciona parámetros cacterísticos del molino: Largo (L), diámetro (D)
'''

def funcion_morrell(fecha,rpm,Pbruta,p_solidos,diccionario_constantes):
    L = float(diccionario_constantes['L'])*0.3048
    d = float(diccionario_constantes['d'])*0.3048
    rho_0 = float(diccionario_constantes['rho_0'])
    rho_b = float(diccionario_constantes['rho_b'])
    g = float(diccionario_constantes['g'])
    E = float(diccionario_constantes['E'])
    U = float(diccionario_constantes['U'])
    rho_p = float(diccionario_constantes['rho_p'])
    Pv = float(diccionario_constantes['Pv'])
    k = float(diccionario_constantes['k'])
    P61 = float(diccionario_constantes['P61'])

    def Nc(J):
        return 42.3/np.sqrt(d)

    #Ecuación 19: velocidad critica teórica
    def phi(J):
        return rpm/Nc(J)

    #Ecuación 17: nivel de llenado de bolas
    def j_b(J):
        jb1 = Pbruta/(np.power(d,2.5)*L)
        jb2 = 0.0348*P61
        jb3 = 4/p_solidos
        jb4 = 6.7*phi(J)
        return (jb1 - jb2 - jb3 + jb4)/15

    #Ecuación de velocidad angular del molino en m/s
    def N_m(J):   
        return rpm/60

    def N(J):
        return N_m(J)/2

    #Ecuación de velocidad crítica
    def phi_c(J):
        return 0.35*(3.364 - J)

    #Ecuación 2:  angulo de pie
    def theta_T(J):
        tt1 = -19.42*(phi_c(J) - phi(J))
        tt2 = 1 - np.exp(tt1)
        tt3 = 2.5307*(1.2796 - J)
        return tt2*tt3 + np.pi/2

    #Ecuación 3:  angulo de hombro
    def theta_S(J):
        ts1 = theta_T(J) - np.pi/2
        ts2 = 0.3386 + 0.1041*phi(J)
        ts3 = (1.54 - 2.5673*phi(J))*J
        return np.pi/2 - ts1*(ts2 + ts3) 

    #Ecuación 14: posición radial media 
    def r(J):
        r1 = (2*np.pi*J)/(2*np.pi + theta_S(J) - theta_T(J))
        r2 = np.sqrt(1 - r1)
        r3 = 1 + r2
        return (d/2)/2*r3

    #Ecuación 16: densidad de carga
    def rho_c(J):
        rc1 = J*rho_0*(1 - E + E*U*rho_p)
        rc2 = j_b(J)*(rho_b - rho_0)*(1 - E)
        rc3 = J*E*U*(1 - rho_p)
        return (rc1 + rc2 + rc3)/J

    def z(J):
        return np.power(1 - J,0.4532)

    #Ecuación 11: tiempo de vuelo
    def t_c(J):
        return (2*np.pi - theta_T(J) + theta_S(J)) / (2*np.pi*N(J))

    #Ecuación 12 : tiempo de free fall
    def t_f(J):
        tf1 = np.sin(theta_S(J)) - np.sin(theta_T(J))
        tf2 = 2*r(J)*tf1/g
        return np.sqrt(tf2)

    #Ecuación 10: parámetro beta
    def beta(J):
        return t_c(J)/(t_f(J) + t_c(J))

    #Ecuación 9: radio interno
    def r_i(J):
        ri1 = (2*np.pi*beta(J)*J)/(2*np.pi + theta_S(J) - theta_T(J))
        ri2 = 1 - ri1
        return (d/2)*np.sqrt(ri2)

    #Ecuación de theta_0 relacionada con el tipo de descarga del molino
    def theta_T0(J):
        return np.pi + np.arcsin(r_i(J)/(d/2) )

    #Ecuación 5: Potencia del molino debido a las energías potencial y cinética
    def P(J):
        a1 = (np.pi*g*L*(d/2) *N_m(J))/(3*(d/2)  - 3*z(J)*r_i(J))
        a2 = 2*np.power((d/2) ,3) - 3*z(J)*np.power((d/2) ,2)*r_i(J) + np.power(r_i(J),3)*(3*z(J) - 2)
        a3 = rho_c(J)*(np.sin(theta_S(J)) - np.sin(theta_T(J)))
        a4 = (N_m(J)*(d/2) *np.pi)/((d/2)  - z(J)*r_i(J))
        a5 = L*rho_c(J)*np.power(a4,3)
        a6 = np.power((d/2)  - z(J)*r_i(J),4) - np.power(r_i(J),4)*np.power(z(J) - 1,4)
        a7 = rho_p*(np.sin(theta_T(J)) - np.sin(theta_T0(J)))
        return a1*a2*(a3+a7)+a5*a6
        #return a1*a2*a3+a5*a6

   #Ecuación 21: Ecuación de Morrell     

    def P_t(J):
        return (Pbruta - Pv)/k 

    def equation_to_solve(J):
        return P(J) - P_t(J) 


    try:
        J_calculado = optimize.bisect(equation_to_solve, 0.1, 0.37, maxiter=100)
        results = {
        'timestamp': fecha,
        'rpm': rpm,
        'potenciaconsumida': int(Pbruta),
        'J': J_calculado*100,
        'Jb': j_b(J_calculado)*100,
        'theta_T': theta_T(J_calculado)*180/np.pi,
        'theta_s' : theta_S(J_calculado)*180/np.pi,
        'vel_critica' : phi(J_calculado)*100,
        'porcentajesolidos': p_solidos,
        'solver_error_exception': None}

        return results

    except RuntimeError:

        results = {
            'timestamp': fecha,
            'rpm': rpm,
            'potenciaconsumida': int(Pbruta),
            'J': np.nan,
            'Jb': np.nan,
            'theta_T': np.nan,
            'theta_s' : np.nan,
            'vel_critica' : np.nan,
            'porcentajesolidos': p_solidos,
            'solver_error_exception': "error"}
        
        return results

