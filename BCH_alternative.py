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
rho = Function(r'\rho')(t)
commutator_subs = {
    Commutator(x, p): I,
    Commutator(p, x): -I,
    Commutator(p, p): 0,
    Commutator(x, x): 0}
#function Lists
func_list_1 = [sin, lambda x: -sin(x), cos, lambda x: -cos(x), 
               lambda x: sin(x)/omega, lambda x: -sin(x)/omega,
               lambda x: -omega*sin(x)]
#Transformation definitions
T_1 = exp((I*rho.diff(t)/(2*rho))*x**2)
T_2 = exp(-I*log(rho)/2*(x*p + p*x))

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
                    break
            else:
                continue
            break
    sumf = 0
    for j in range(len(coeff_list)):
        sumf += coeff_list[j]*operator_list[j]
    return sumf
#--------------------------------------------
def extract_operator_and_coeff(expression):
    """
    Extracts the operator and coefficient from an expression of the form exp(coeff * operator).

    Args:
        expression: The expression to extract from.

    Returns:
        A tuple containing the operator and the coefficient.
    """

    if not expression.func == exp:
        raise ValueError("Expression must be an exponential function")

    arg = expression.args[0]

    # Handle simple cases:
    if arg.is_Mul:
        operator = arg.args[-1]
        coeff = arg / operator
        return operator, coeff
    elif arg.is_Mul or arg.is_Add:
        # Handle more complex expressions:
        # Use pattern matching or other techniques to identify operators and coefficients
        # ... (Implement more complex pattern matching or symbolic manipulation here)
        raise NotImplementedError("Complex expressions not yet supported")

    return arg, 1  # If no operator is found, return the entire argument as the operator
#-----------------------------------------------------        
def operadores_inmiscuidos(S, x):
    """
    Calculates the sandwich <S|x|S> using the BCH formula.

    Args:
        S: The transformation operator.
        x: The operator to be sandwiched.

    Returns:
        The BCH expression for the sandwiched operator.
    """

    operator, coeff = extract_operator_and_coeff(S)
    bch_expression = BCH(-operator, x, coeff, n=9)

    # Simplify the BCH expression
    simplified_expr = bch_expression.simplify()

    return simplified_expr.coeff(x)
#----------------------------------------------
def single_find(expression,func_list,var_list):
    for i in range(len(func_list)):
        for j in range(len(var_list)):
            for k in [-1,0,1]:
                if expression == series(func_list[i](t),t,x0 = 0,n = 9 + k).removeO().subs({t:var_list[j]}):
                    expression = func_list[i](t).subs({t:var_list[j]})
                    break
            else:
                continue
            break
        else:
            continue
        break
    return expression
#---------------------------------------------
def substitute_operators(expr, substitutions):
    """
    Substitutes operators in an expression with their corresponding functions.

    Args:
        expr: The expression to substitute in.
        substitutions: A dictionary mapping operators to their substitutions.

    Returns:
        The expression with substituted operators.
    """

    # Recursively substitute operators in the expression
    def _substitute_recursive(expr, substitutions):
        if isinstance(expr, Add):
            return Add(*[_substitute_recursive(arg, substitutions) for arg in expr.args])
        elif isinstance(expr, Mul):
            return Mul(*[_substitute_recursive(arg, substitutions) for arg in expr.args])
        elif isinstance(expr, Pow):
            base, exponent = expr.args
            return Pow(_substitute_recursive(base, substitutions), exponent)
        elif expr in substitutions:
            return substitutions[expr]
        else:
            return expr

    return _substitute_recursive(expr, substitutions)
#------------------------------------------------------
def time_dep_exp(BCH,transformation):
  return BCH -I*transformation.args[0].diff(t)
def calculate_time_dep_exp(exponential, hamiltonian, n=9):
    """
    Calculates the time-dependent exponential of a Hamiltonian.

    Args:
    exponential: The exponential expression.
    hamiltonian: The Hamiltonian.
    n: The order of the BCH expansion.

    Returns:
    The time-dependent exponential expression.
    """
    operator, coeff = extract_operator_and_coeff(exponential)
    bch_expression = BCH(-operator, hamiltonian, coeff, n)
    time_dep_expression = expand(time_dep_exp(bch_expression, exponential))
    return time_dep_expression
        
#test zone
#BCH_test = BCH(H,x,alpha,9).subs({alpha:I*t})
#BCH_test = BCH(H,p,alpha,9).subs({alpha:I*t})
#test = extract_coeff(BCH(H,x,alpha,9).subs({alpha:I*t}),x)
#test = com_order(H,x,8)
#test = series(cos(omega*t),omega*t,x0 = 0, n = 9).removeO()
#test = extract_coeff(BCH(H,x,alpha,9).subs({alpha:I*t}),x)
#test = find_and_replace(BCH_test,[x,p],func_list_1)
#test = extract_operator_and_coeff(T_1)
#test = extract_operator_and_coeff(T_2)
#test2 = operadores_inmiscuidos(T_2,x)
#test = series(exp(t),t,x0 = 0,n = 10).removeO().subs({t:log(rho)})
#test = single_find(test2,[exp],[log(rho)])
#print(print_latex(test))
test = calculate_time_dep_exp(T_1,H)
#-------------------------------------------------
"""H = (p**2 + omega**2*x**2 + x*p + p*x)/2
substitutions = {x: cos(omega*t), p: -I*omega*sin(omega*t)}
new_H = substitute_operators(H, substitutions)"""
#------------------------------------------------
print(print_latex(test))
