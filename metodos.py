import mpmath
from sympy import *	
x, y, z, t = symbols("x y z t")

class methods:
	def euler(entrie): # Pronto
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Euler\n")
		output.write("y( 0.0 ) = " + str(y0) + "\n")
		output.write("h = " + str(h) + "\n")

		Y = y0
		T = t0

		for i in range(1, (steps + 1)):
			# Yn+1 = Yn + h*fn
			F = h*expression.subs([(t, T), (y, Y)]) # Substituir o y da equacao pela variavel Y ( Assim para o t e T )

			Y += F
			T += h

			output.write(str(i) + " " + str(Y) + "\n")
		output.write("\n")
		output.close()

	def eulerInverse(entrie): # Pronto
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Euler Inverso\n")
		output.write("y( 0.0 ) = " + str(y0) + "\n")
		output.write("h = " + str(h) + "\n")
		
		Y = y0
		T = t0

		for i in range(1, (steps + 1)):
			# Yn+1 = Yn + h *F(tn+1, yn+1)
			# Previsoes 
			Tn1 = T + h # Equivalente ao tn+1
			Yn1 = Y + h*expression.subs([(t, T), (y, Y)]) # Equivalente ao yn+1

			Y += h*expression.subs([(t, Tn1), (y, Yn1)])
			T += h

			output.write(str(i) + " " + str(Y) + "\n")
		output.write("\n")
		output.close()

	def eulerImproved(entrie): # Pronto
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Euler Aprimorado\n")
		output.write("y( 0.0 ) = " + str(y0) + "\n")
		output.write("h = " + str(h) + "\n")
		

		Y = y0
		T = t0

		for i in range(1, (steps + 1)):
			# Y = Y * h((F1+F2)/2)

			F1 = expression.subs([(t, T), (y, Y)])
			Yh = Y + h*F1 # Utilizando o metodo de Euler pra calcular o Yn+1
			F2 = expression.subs([(t, T+h), (y, Yh)]) 

			Y += h*((F1+F2)/2)
			T += h

			output.write(str(i) + " " + str(Y) + "\n")
		output.write("\n")
		output.close()

	def rungeKutta(entrie): # Pronto
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Runge-Kutta\n")
		output.write("y( 0.0 ) = " + str(y0) + "\n")
		output.write("h = " + str(h) + "\n")

		Y = y0
		T = t0

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

			output.write(str(i) + " " + str(Y) + "\n")
		output.write("\n")
		output.close()
		
	def adamBashforth(entrie):
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		
		degree = entrie.pop()
		degree = int(degree)
		stringExpression = entrie.pop()
		X = entrie[3:degree+2]

		expression = sympify(stringExpression)
		steps = int(entrie.pop())
		h = float(entrie.pop())

		#output = open('saida.txt', 'a')
		#output.write("Metodo de Adam-Bashforth\n")
		#output.write("y( 0.0 ) = " + str(y0) + "\n")
		#output.write("h = " + str(h) + "\n")

		#Y = y0
		#T = t0
		#F[Y] = answer
		#F[0][0] = 0

		#for i in range(1, (steps + 1)):
			# Yn+1 = Yn + 3/2 * h * fn - 1/2 * h * fn-1
			#F0 = (3/2)*h*expression.subs([(t, T), (y, Y)]) 
			#F1 = (-1)*(1/2)*h*expression.subs([(t, T), (y, Y)]) 
			

			#output.write(str(i) + " " + str(Y) + "\n")
		#output.write("\n")
		#output.close()

	#def adamMulton():
		
	#def inverseFormula():
		
def switchCase(strct):
	if strct[0] == "euler":
		methods.euler(strct)
	elif strct[0] == "euler_inverso":
		methods.eulerInverse(strct)
	elif strct[0] == "euler_aprimorado":
		methods.eulerImproved(strct)
	elif strct[0] == "runge_kutta":
		methods.rungeKutta(strct)
	elif strct[0] == "adam_bashforth":
		methods.adamBashforth(strct)

def main():
	output = open('saida.txt', 'w') # Criar o arquivo vazio inicialmente
	output.close()

	Input = open('entradas.txt', 'r')
	data = Input.readlines()

	for line in data:
		array = line.split()
		switchCase(array)

	Input.close()
main()
