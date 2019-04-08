import mpmath
from sympy import *	
x, y, z, t = symbols("x y z t")

menuLoop = True

def debug(method, y0, t0, h, expression):
	print("metodo =" + method + "\ny0 = " + str(y0) + "\nt0 = " + str(t0) + "\nh = " + str(h) + "\nexpression = " + str(expression))

class methods:
	def euler(): # Pronto
		y0 = float(input("Insira o valor de y0: "))
		t0 = float(input("Insira o valor de t0: "))
		h = float(input("Insira o valor de h: "))
		steps = int(input("Insira o numero de passos: "))

		stringExpression = input("Insira a expressao: ")
		expression = sympify(stringExpression)

		Y = y0
		T = t0

		print("0 " + str(y0))
		for i in range(1, (steps + 1)):
			F = h*expression.subs([(t, T), (y, Y)]) # Substituir o y da equacao pela variavel Y ( Assim para o t e T )

			Y += F
			T += h

			print(str(i), Y)

	def eulerInverse(): # Pronto
		y0 = float(input("Insira o valor de y0: "))
		t0 = float(input("Insira o valor de t0: "))
		h = float(input("Insira o valor de h: "))
		steps = int(input("Insira o numero de passos: "))

		stringExpression = input("Insira a expressao: ")
		expression = sympify(stringExpression)
		
		Y = y0
		T = t0

		print("0 " + str(y0))
		for i in range(1, (steps + 1)):
			# Yn+1 = Yn + h *F(tn+1, yn+1)
			# Previsoes 
			Tn1 = T + h # Equivalente ao tn+1
			Yn1 = Y + h*expression.subs([(t, T), (y, Y)]) # Equivalente ao yn+1

			Y += h*expression.subs([(t, Tn1), (y, Yn1)])
			T += h

			print(str(i), Y)

	def eulerImproved(): # Pronto
		print("Metodo: euler aprimorado")
		y0 = float(input("Insira o valor de y0: "))
		t0 = float(input("Insira o valor de t0: "))
		h = float(input("Insira o valor de h: "))
		steps = int(input("Insira o numero de passos: "))

		stringExpression = input("Insira a expressao: ")
		expression = sympify(stringExpression)

		Y = y0
		T = t0

		print("0 " + str(y0))
		for i in range(1, (steps + 1)):
			# Y = Y * h((F1+F2)/2)

			F1 = expression.subs([(t, T), (y, Y)])
			Yh = Y + h*F1 # Utilizando o metodo de Euler pra calcular o Yn+1
			F2 = expression.subs([(t, T+h), (y, Yh)]) 

			Y += h*((F1+F2)/2)
			T += h

			print(str(i), Y)

	def rungeKutta(): # Pronto
		print("Metodo: euler aprimorado")
		y0 = float(input("Insira o valor de y0: "))
		t0 = float(input("Insira o valor de t0: "))
		h = float(input("Insira o valor de h: "))
		steps = int(input("Insira o numero de passos: "))

		stringExpression = input("Insira a expressao: ")
		expression = sympify(stringExpression)

		Y = y0
		T = t0

		print("0 " + str(y0))
		for i in range(1, (steps + 1)):
			# Y = Y * h((F1+F2)/2)

			F1 = expression.subs([(t, T), (y, Y)])
			Yh = Y + h*F1 # Utilizando o metodo de Euler pra calcular o Yn+1
			F2 = expression.subs([(t, T+h), (y, Yh)]) 

			Y += h*((F1+F2)/2)
			T += h

			print(str(i), Y)
		y0 = float(input("Insira o valor de y0: "))
		t0 = float(input("Insira o valor de t0: "))
		h = float(input("Insira o valor de h: "))
		steps = int(input("Insira o numero de passos: "))

		stringExpression = input("Insira a expressao: ")
		expression = sympify(stringExpression)

		Y = y0
		T = t0

		print("0 " + str(y0))
		for i in range(1, (steps + 1)):
			# Y += Y + ((K1+ 2*K2 + 2*K3 + K4)/6)*h 

			K1 = expression.subs([(t, T), (y, Y)])

			T2 = T + (h/2)
			Y2 = Y+(h/2)*K1
			K2 = expression.subs([(t, T2), (y, Y2)])

			T3 = T + (h/2)
			Y3 = Y+(h/2)*K2
			K3 = expression.subs([(t, T3), (y, Y3)])
			
			T4 = T + h
			Y4 = Y + h*K3
			K4 = expression.subs([(t, T4), (y, Y4)])

			K = (K1+(2*K2)+(2*K3)+K4)/6

			Y += K*h
			T += h

			print(str(i), Y)
		
	def adamBashforth():
		
	#def adamMulton():
		
	#def inverseFormula():
		

def menuOptions():
	print("Escolha o metodo a ser utilizado a baixo pelo seu codigo:\n\nEuler            = euler\nEuler Inverso    = euler_inverso\nEuler Aprimorado = euler_aprimorado\n\nRunge-Kutta      = runge_kutta\n")
	print("Adam-Bashforth   = adam_bashforth\n                 = adam_bashforth_by_euler\n                 = adam_bashforth_by_euler_inverso\n                 = adam_bashforth_by_euler_aprimorado\n                 = adam_bashforth_by_runge_kutta\n\n")
	print("Adam-Multon      = adam_multon\n                 = adam_multon_by_euler\n                 = adam_multon_by_euler_inverso\n                 = adam_multon_by_euler_aprimorado\n                 = adam_multon_by_runge_kutta\n\n")
	print("Formula Inversa  = formula_inversa\n                 = formula_inversa\n                 = formula_inversa_by_euler\n                 = formula_inversa_by_euler_inverso\n                 = formula_inversa_by_euler_aprimorado\n                 = formula_inversa_by_runge_kutta\n\n")
	print("Sair do menu     = sair")
	return input("Opcao: ")

def switchCase(string):
	if string == "euler":
		methods.euler()
	elif string == "euler_inverso":
		methods.eulerInverse()
	elif string == "euler_aprimorado":
		methods.eulerImproved()
	elif string == "runge_kutta":
		methods.rungeKutta()
	elif string == "adam_bashforth"
		methods.adam_bashforth()
	elif string == "sair":
		exit()

def main():
	while menuLoop == True:
		methodOption = menuOptions()
		switchCase(methodOption)
main()
