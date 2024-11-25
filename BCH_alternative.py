#Baker-Campbell-Hausdorff python Implementation for the Simple Harmonic.expand(ommutator=True) Oscillator

#IMPORTS
from sympy import *
from sympy.series.series import series
from sympy.physics.quantum import *
from sympy.printing.latex import LatexPrinter, print_latex
init_printing()
#------------------------------------------------------
#early definitions
t = symbols(r't',real = True)
x = Operator(r'\hat{x}')
p = Operator(r'\hat{p}')
omega = symbols(r'\omega')
H = (p**2+omega**2*x**2)/2
alpha = symbols(r'\alpha')
commutator_subs = {
    Commutator(x, p): I,
    Commutator(p, x): -I,
    Commutator(p, p): 0,
    Commutator(x, x): 0}
#--------------------------------------------------------
func_list_1 = [sin, lambda x: -sin(x), cos, lambda x: -cos(x), 
               lambda x: sin(x)/omega, lambda x: -sin(x)/omega]
#----------------------------------------------------------
def conmutador(A, B):
    result = Commutator(A, B).expand(commutator=True).expand(commutator=True).subs(commutator_subs)
    return result
#-----------------------------------------------------------
def com_order(A, B, n):
    if n <= 1:
        return conmutador(A, B)
    result = com_order(A, conmutador(A, B), n - 1)
    return result
#---------------------------------------------------------
def BCH(A, B, alpha, n):
    factorials = [factorial(i) for i in range(n+1)]
    a = 0
    for i in range(1, n+1):
        alpha_pow_i = alpha**i  # Calculate alpha^i explicitly
        a += (alpha_pow_i / factorials[i]) * com_order(A, B, i)
    return B + a  
#---------------------------------------------------------
def extract_coeff(expression,Operator):
    return expression.coeff(Operator)
#------------------------------------
def find_and_replace(expression,operator_list,function_list):
    coeff_list = []
    offset = [-1,0,1]
    for i in operator_list:
        coeff_list.append(expression.coeff(i))
    for i in range(len(coeff_list)):
        for j in range(len(function_list)):
            for k in offset:
                if coeff_list[i] == series(function_list[j](omega*t),omega*t,x0 = 0,n = 9 + k).removeO():
                    coeff_list[i] = function_list[j](omega*t)
    sumf = 0
    for j in range(len(coeff_list)):
        sumf += coeff_list[j]*operator_list[j]
    return sumf

        
#test zone
BCH_test = BCH(H,x,alpha,9).subs({alpha:I*t})
#test = extract_coeff(BCH(H,x,alpha,9).subs({alpha:I*t}),x)
#test = com_order(H,x,8)
#test = series(cos(omega*t),omega*t,x0 = 0, n = 9).removeO()
#test = extract_coeff(BCH(H,x,alpha,9).subs({alpha:I*t}),x)
test = find_and_replace(BCH_test,[x,p],func_list_1)
"""a = BCH_test.coeff(x)
print(a)
lol = (a == series(cos(omega*t),omega*t,x0 = 0,n = 9).removeO())
if lol == True:
    a = cos(omega*t)
print(a)"""
print(print_latex(test))
