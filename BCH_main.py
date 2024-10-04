from operator import indexOf
import sympy
from sympy import symbols
from sympy.plotting import plot
sympy.init_printing()
from sympy import*
from sympy.physics.quantum import *
from sympy.physics.quantum import Operator, Commutator, AntiCommutator
from sympy.physics.quantum import Dagger, Commutator
from sympy.physics.quantum.boson import BosonOp
from math import ceil
from sympy.printing.latex import LatexPrinter, print_latex
#-----------------------------
#-----------------------------

def collect_coeff(expresion,A): #factoriza el operador A de una expresión "expresion"
  vt = 0
  for k in range(len(expand(expresion).args)):
    if (expand(expresion).args[k]).coeff(A)!=0:
      vt+= (expand(expresion).args[k].coeff(A))
  return (vt)*A
#--------------------------------------
def collect_coeff2(expresion,A): #factoriza el operador A de una expresión "expresion" NOTA AL DR. RAMOS: MODIFIQUÉ LA EXPRESIÓN FINAL: AHORA REGRESA (vt), antes regresaba (vt)*A
  vt = 0
  for k in range(len(expand(expresion).args)):
    if (expand(expresion).args[k]).coeff(A)!=0:
      vt+= (expand(expresion).args[k].coeff(A))
  return (vt)
#-----------------------------
def conmutador(A,B):
  return Commutator(A,B).expand(commutator=True).expand(commutator=True).subs({Commutator(x,p):I, Commutator(p,x):-I,Commutator(p,p):0,Commutator(x,x):0,

                                                                               Commutator(a,a_dag):1,Commutator(a_dag,a):-1,Commutator(a_dag,a_dag):0,Commutator(a,a):0,

                                                                               Commutator(sigmap,sigmam):sigmaz,Commutator(sigmam,sigmap):-sigmaz,
                                                                               Commutator(sigmaz,sigmap):2*sigmap,
                                                                               Commutator(sigmap,sigmaz):-2*sigmap,Commutator(sigmam,sigmaz):2*sigmam,
                                                                               Commutator(sigmaz,sigmam):-2*sigmam,
                                                                               Commutator(sigmap,sigmap):0,Commutator(sigmam,sigmam):0,Commutator(sigmaz,sigmaz):0})
def com_order(A,B,n): #regresa el conmutador anidado [A,[A,..[A,B]] n veces.
    if n<=1:
        return conmutador(A,B)

    else:
        return com_order(A,conmutador(A,B),n-1)
#-----------------------------
#AQUI IBA LAS FUNCIONES RELACIONADAS CON LATEX
#------------------------------
def BCH(A, B, alpha, n): # Baker_Campbell_Hausdorff formula
    a=0
    for i in range(1,n+1):
        a += (alpha**i/factorial(i))*com_order(A,B,i)
    return (B + a)
#-------------------------------
def collect_all(expression,operator_list): # itera "len(operator list)" veces la función "collect_coef" y suma todas sus iteraciones.
  a=0
  for i in range(len(operator_list)):
    a += collect_coeff(expression,operator_list[i])
  return a
#-------------------------------
def collect_and_save_all(expression,operator_list):
  coeff_list=[] #lista donde vamos a meter los coeficientes de cada operador
  for z in range(len(operator_list)): #iterador para cada operador
    a = collect_coeff2(expression,operator_list[z])
    coeff_list.append(a)
  return coeff_list
#-------------------------------
def compare_all(coeff_list,operator_list,n): #compara los coeficientes y los reemplaza.
  list_1= coeff_list.copy()
  functions = [cos, sin , exp] #lista de funciones para comparar
  for z in range(len(operator_list)):
    for i in functions: #iterador para cada funcion
      for y in (-2,1,0,1,2): #iterador para el orden de las expansiones
        if coeff_list[z] == series(i(omega*t),omega*t,0,n+y).removeO(): #funciones en 'functions' positivas
          list_1[z] = i(omega*t)
          break
        elif coeff_list[z] == series(-i(omega*t),omega*t,0,n+y).removeO():# funciones en 'functions' negativas
          list_1[z] = -i(omega*t)
          break
        elif coeff_list[z] == fraction(series(i(omega*t)/omega,omega,0,n+y).removeO())[0]: #funciones en 'functions' cardinales positivas
          list_1[z] = i(omega*t)/omega
          break
        elif coeff_list[z] == fraction(series(-i(omega*t)/omega,omega,0,n+y).removeO())[0]:#funciones en 'functions' cardinales negativas
          list_1[z] = -i(omega*t)/omega
          break
        elif coeff_list[z] == series(i(omega*t)*omega,omega,0,n+y).removeO(): #funciones en 'functions' cardinales positivas
          list_1[z] = i(omega*t)*omega
          break
        elif coeff_list[z] == series(-i(omega*t)*omega,omega,0,n+y).removeO():#funciones en 'functions' cardinales negativas
          list_1[z] = -i(omega*t)*omega
          break

  return list_1

#-------------------------------
def sum_collect(coef_list,operator_list):
  a=0
  for i in range(len(operator_list)):
    a += coef_list[i]*(operator_list[i])
  return a
#-------------------------------

def differentiate(T,t):
  return T.diff(t)
 #-------------------------------
def diff_term(T,t):
  return - I*Dagger(T)*T.diff(t)
#-------------------------------

  


x     = Operator(r'\hat{x}')
p     = Operator(r'\hat{p}')
omega = symbols(r'\omega', real = True)
alpha = symbols(r'\alpha')
t     = symbols(r't', real = True)
H     = Rational(1,2)*(p**2+omega**2*x**2)
l =[x,p]
sigmap = Operator(r'\hat{\sigma}_+')
sigmam = Operator(r'\hat{\sigma}_-')
sigmaz = Operator(r'\hat{\sigma}_z')

H_2 = sigmaz + sigmap + sigmam
l_2 = [sigmaz,sigmap,sigmam]
#---------------------------------------
a_dag = Operator(r'\hat{a}^{\dagger}')
a = Operator(r'\hat{a}')
H_3 = omega*(a_dag*a+1/2)
l_3 = [a,a_dag]
#---------------------------------------
orden = 8
