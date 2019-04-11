import mpmath
from sympy.plotting import plot as symplot
from sympy import *	
x, y, z, t = symbols("x y z t")

class methods:
	def euler(entrie):
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Euler\n")
		output.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
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

	def eulerInverse(entrie): 
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Euler Inverso\n")
		output.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
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

	def eulerImproved(entrie):
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Euler Aprimorado\n")
		output.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
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

	def rungeKutta(entrie):
		y0 = float(entrie[1])
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])

		stringExpression = entrie[5]
		expression = sympify(stringExpression)

		output = open('saida.txt', 'a')
		output.write("Metodo de Runge-Kutta\n")
		output.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
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
		degree = int(entrie.pop())
		Y = []
		
		t0 = float(entrie[degree + 1])
		h = float(entrie[degree + 2])
		steps = int(entrie[degree + 3])
		expression = entrie[degree + 4]

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Bashforth de ordem " + str(degree) + "\n")
		output.write("y( 0.0 ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		
		for i in range(1, degree + 1):
			Y.append(float(entrie[i]))
			output.write(str(i-1) + " " + entrie[i] + "\n")
			
		expression = sympify(expression)

		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps):
			if degree == 2:
				F1 = (3/2)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/2)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F1 + F0)
			elif degree == 3:
				F2 = (23/12)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-4/3)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (5/12)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F2 + F1 + F0)
			elif degree == 4:
				F3 = (55/24)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-59/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (37/24)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (-3/8)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
			elif degree == 5:
				F4 = (1901/720)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-1387/360)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (109/30)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-637/360)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (251/720)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
			elif degree == 6:
				F5 = (4277/1440)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-2641/480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (4991/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-3649/720)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (959/480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-95/288)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 7:
				F6 = (198721/60480)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-18637/2520)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (235183/20160)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-10754/945)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (135713/20160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (19087/60480)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)

			elif degree == 8:
				F7 = (16083/4480)*expression.subs([(t, T), (y, Y[i])])
				F6 = (-1152169/120960)*expression.subs([(t, T-h), (y, Y[i-1])])
				F5 = (242653/13440)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F4 = (-296053/13440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F3 = (2102243/120960)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F2 = (-115747/13440)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F1 = (32863/13440)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])
				F0 = (-5257/17280)*expression.subs([(t, T-(7*h)), (y, Y[i-7])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h

		output.write("\n")
		output.close()

	def adamBashforthEuler(entrie): 
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Bashforth por Euler de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yEuler = y0
		tEuler = t0

		for i in range(1, degree):
			F = h*expression.subs([(t, tEuler), (y, yEuler)])

			yEuler += F
			tEuler += h
			Y.append(yEuler)
			output.write(str(i) + " " + str(yEuler) + "\n")

		T = t0 + (degree - 1) * h				

		for i in range(degree - 1, steps):
			if degree == 2:
				F1 = (3/2)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/2)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F1 + F0)
			elif degree == 3:
				F2 = (23/12)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-4/3)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (5/12)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F2 + F1 + F0)
			elif degree == 4:
				F3 = (55/24)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-59/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (37/24)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (-3/8)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
			elif degree == 5:
				F4 = (1901/720)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-1387/360)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (109/30)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-637/360)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (251/720)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
			elif degree == 6:
				F5 = (4277/1440)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-2641/480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (4991/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-3649/720)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (959/480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-95/288)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 7:
				F6 = (198721/60480)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-18637/2520)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (235183/20160)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-10754/945)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (135713/20160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (19087/60480)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def adamBashforthEulerInverse(entrie): 
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Bashforth por Euler Inverso de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yEulerInverse = y0
		tEulerInverse = t0

		for i in range(1, degree):
			Tn1 = tEulerInverse + h # Equivalente ao tn+1
			Yn1 = yEulerInverse + h*expression.subs([(t, tEulerInverse), (y, yEulerInverse)]) # Equivalente ao yn+1

			yEulerInverse += h*expression.subs([(t, Tn1), (y, Yn1)])
			tEulerInverse += h
			Y.append(yEulerInverse)
			output.write(str(i) + " " + str(yEulerInverse) + "\n")

		T = t0 + (degree - 1) * h				

		for i in range(degree - 1, steps):
			if degree == 2:
				F1 = (3/2)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/2)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F1 + F0)
			elif degree == 3:
				F2 = (23/12)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-4/3)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (5/12)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F2 + F1 + F0)
			elif degree == 4:
				F3 = (55/24)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-59/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (37/24)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (-3/8)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
			elif degree == 5:
				F4 = (1901/720)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-1387/360)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (109/30)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-637/360)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (251/720)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
			elif degree == 6:
				F5 = (4277/1440)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-2641/480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (4991/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-3649/720)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (959/480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-95/288)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 7:
				F6 = (198721/60480)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-18637/2520)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (235183/20160)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-10754/945)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (135713/20160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (19087/60480)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 8:
				F7 = (16083/4480)*expression.subs([(t, T), (y, Y[i])])
				F6 = (-1152169/120960)*expression.subs([(t, T-h), (y, Y[i-1])])
				F5 = (242653/13440)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F4 = (-296053/13440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F3 = (2102243/120960)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F2 = (-115747/13440)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F1 = (32863/13440)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])
				F0 = (-5257/17280)*expression.subs([(t, T-(7*h)), (y, Y[i-7])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h# Pronto
		output.write("\n")
		output.close()

	def adamBashforthEulerImproved(entrie): 
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Bashforth por Euler Aprimorado de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yEulerImproved = y0
		tEulerImproved = t0

		for i in range(1, degree):
			F1_improved = expression.subs([(t, tEulerImproved), (y, yEulerImproved)])
			Yh = yEulerImproved + h*F1_improved # Utilizando o metodo de Euler pra calcular o Yn+1
			F2_improved = expression.subs([(t, tEulerImproved+h), (y, Yh)]) 

			yEulerImproved += h*((F1_improved+F2_improved)/2)
			tEulerImproved += h
			Y.append(yEulerImproved)
			output.write(str(i) + " " + str(yEulerImproved) + "\n")

		T = t0 + (degree - 1) * h				

		for i in range(degree - 1, steps):
			if degree == 2:
				F1 = (3/2)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/2)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F1 + F0)
			elif degree == 3:
				F2 = (23/12)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-4/3)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (5/12)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F2 + F1 + F0)
			elif degree == 4:
				F3 = (55/24)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-59/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (37/24)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (-3/8)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
			elif degree == 5:
				F4 = (1901/720)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-1387/360)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (109/30)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-637/360)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (251/720)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
			elif degree == 6:
				F5 = (4277/1440)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-2641/480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (4991/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-3649/720)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (959/480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-95/288)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 7:
				F6 = (198721/60480)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-18637/2520)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (235183/20160)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-10754/945)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (135713/20160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (19087/60480)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 8:
				F7 = (16083/4480)*expression.subs([(t, T), (y, Y[i])])
				F6 = (-1152169/120960)*expression.subs([(t, T-h), (y, Y[i-1])])
				F5 = (242653/13440)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F4 = (-296053/13440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F3 = (2102243/120960)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F2 = (-115747/13440)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F1 = (32863/13440)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])
				F0 = (-5257/17280)*expression.subs([(t, T-(7*h)), (y, Y[i-7])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def adamBashforthRungeKutta(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Bashforth por Runge-Kutta de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yRungeKutta = y0
		tRungeKutta = t0

		for i in range(1, degree):
			K1 = expression.subs([(t, tRungeKutta), (y, yRungeKutta)])

			T2 = tRungeKutta + (h/2)
			Y2 = yRungeKutta+(h/2)*K1
			K2 = expression.subs([(t, T2), (y, Y2)])

			T3 = tRungeKutta + (h/2)
			Y3 = yRungeKutta+(h/2)*K2
			K3 = expression.subs([(t, T3), (y, Y3)])
			
			T4 = tRungeKutta + h
			Y4 = yRungeKutta + h*K3
			K4 = expression.subs([(t, T4), (y, Y4)])

			K = (K1+(2*K2)+(2*K3)+K4)/6

			yRungeKutta += K*h
			tRungeKutta += h
			Y.append(yRungeKutta)
			output.write(str(i) + " " + str(yRungeKutta) + "\n")

		T = t0 + (degree - 1) * h				

		for i in range(degree - 1, steps):
			if degree == 2:
				F1 = (3/2)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/2)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F1 + F0)
			elif degree == 3:
				F2 = (23/12)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-4/3)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (5/12)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F2 + F1 + F0)
			elif degree == 4:
				F3 = (55/24)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-59/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (37/24)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (-3/8)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
			elif degree == 5:
				F4 = (1901/720)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-1387/360)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (109/30)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-637/360)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (251/720)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
			elif degree == 6:
				F5 = (4277/1440)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-2641/480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (4991/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-3649/720)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (959/480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-95/288)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 7:
				F6 = (198721/60480)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-18637/2520)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (235183/20160)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-10754/945)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (135713/20160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (19087/60480)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
			elif degree == 8:
				F7 = (16083/4480)*expression.subs([(t, T), (y, Y[i])])
				F6 = (-1152169/120960)*expression.subs([(t, T-h), (y, Y[i-1])])
				F5 = (242653/13440)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F4 = (-296053/13440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F3 = (2102243/120960)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F2 = (-115747/13440)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F1 = (32863/13440)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])
				F0 = (-5257/17280)*expression.subs([(t, T-(7*h)), (y, Y[i-7])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def adamMulton(entrie):
		degree = int(entrie.pop())
		Y = []
		
		t0 = float(entrie[degree])
		h = float(entrie[degree + 1])
		steps = int(entrie[degree + 2])
		expression = entrie[degree + 3]

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Multon de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		
		for i in range(1, degree):
			Y.append(float(entrie[i]))
			output.write(str(i-1) + " " + entrie[i] + "\n")
		expression = sympify(expression)

		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps): 
			if degree == 1:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F1 = (1/2)*expression.subs([(t, Tn1), (y, Yn1)])
				F0 = (1/2)*expression.subs([(t, T), (y, Y[i])])

				yni = Y[i] + h*(F1 + F0)
				yni += (-1/12)
			elif degree == 2:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
				F1 = (2/3)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/12)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F2 + F1 + F0)
				yni += (1/24)
			elif degree == 3:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
				F2 = (19/24)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-5/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 4:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
				F3 = (323/360)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-11/30)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])

				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 5:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F5 = (95/288)*expression.subs([(t, Tn1), (y, Yn1)])
				F4 = (1427/1440)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-133/240)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-863/60480)
			elif degree == 6:
				Tn1 = T + h 
				Yn1 = Y[i-1] + h*expression.subs([(t, T), (y, Y[i-1])]) 

				F6 = (19087/60480)*expression.subs([(t, Tn1), (y, Yn1)])
				F5 = (2713/2520)*expression.subs([(t, T), (y, Y[i-1])])
				F4 = (-15487/20160)*expression.subs([(t, T-h), (y, Y[i-2])])
				F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, Y[i-3])])
				F2 = (-6737/20160)*expression.subs([(t, T-(3*h)), (y, Y[i-4])])
				F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, Y[i-5])])
				F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, Y[i-6])])
				
				yni = Y[i-1] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (275/24192)
			elif degree == 7:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F7 = (5257/17280)*expression.subs([(t, Tn1), (y, Yn1)])
				F6 = (139849/120960)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-4511/4480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-88547/120960)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-33953/3628800)
			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close() # Pronto# Pronto \

	def adamMultonEuler(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Multon por Euler de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yEuler = y0
		tEuler = t0

		for i in range(1, degree):
			F = h*expression.subs([(t, tEuler), (y, yEuler)])

			yEuler += F
			tEuler += h
			Y.append(yEuler)
			output.write(str(i) + " " + str(yEuler) + "\n")

		T = t0 + (degree - 1) * h			

		for i in range(degree - 1, steps): 
			if degree == 1:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F1 = (1/2)*expression.subs([(t, Tn1), (y, Yn1)])
				F0 = (1/2)*expression.subs([(t, T), (y, Y[i])])

				yni = Y[i] + h*(F1 + F0)
				yni += (-1/12)
			elif degree == 2:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
				F1 = (2/3)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/12)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F2 + F1 + F0)
				yni += (1/24)
			elif degree == 3:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
				F2 = (19/24)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-5/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 4:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
				F3 = (323/360)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-11/30)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])

				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 5:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F5 = (95/288)*expression.subs([(t, Tn1), (y, Yn1)])
				F4 = (1427/1440)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-133/240)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-863/60480)
			elif degree == 6:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F6 = (19087/60480)*expression.subs([(t, Tn1), (y, Yn1)])
				F5 = (2713/2520)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-15487/20160)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-6737/20160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (275/24192)
			elif degree == 7:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F7 = (5257/17280)*expression.subs([(t, Tn1), (y, Yn1)])
				F6 = (139849/120960)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-4511/4480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-88547/120960)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-33953/3628800)

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()# Pronto \

	def adamMultonEulerInverse(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Multon por Euler Inverso de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yEulerInverse = y0
		tEulerInverse = t0

		for i in range(1, degree):
			Tn1 = tEulerInverse + h # Equivalente ao tn+1
			Yn1 = yEulerInverse + h*expression.subs([(t, tEulerInverse), (y, yEulerInverse)]) # Equivalente ao yn+1

			yEulerInverse += h*expression.subs([(t, Tn1), (y, Yn1)])
			tEulerInverse += h
			Y.append(yEulerInverse)
			output.write(str(i) + " " + str(yEulerInverse) + "\n")

		T = t0 + (degree - 1) * h			

		for i in range(degree - 1, steps): 
			if degree == 1:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F1 = (1/2)*expression.subs([(t, Tn1), (y, Yn1)])
				F0 = (1/2)*expression.subs([(t, T), (y, Y[i])])

				yni = Y[i] + h*(F1 + F0)
				yni += (-1/12)
			elif degree == 2:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
				F1 = (2/3)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/12)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F2 + F1 + F0)
				yni += (1/24)
			elif degree == 3:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
				F2 = (19/24)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-5/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 4:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
				F3 = (323/360)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-11/30)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])

				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 5:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F5 = (95/288)*expression.subs([(t, Tn1), (y, Yn1)])
				F4 = (1427/1440)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-133/240)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-863/60480)
			elif degree == 6:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F6 = (19087/60480)*expression.subs([(t, Tn1), (y, Yn1)])
				F5 = (2713/2520)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-15487/20160)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-6737/20160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (275/24192)
			elif degree == 7:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F7 = (5257/17280)*expression.subs([(t, Tn1), (y, Yn1)])
				F6 = (139849/120960)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-4511/4480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-88547/120960)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-33953/3628800)

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()# Pronto \

	def adamMultonEulerImproved(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Multon por Euler Aprimorado de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yEulerImproved = y0
		tEulerImproved = t0

		for i in range(1, degree):
			F1_improved = expression.subs([(t, tEulerImproved), (y, yEulerImproved)])
			Yh = yEulerImproved + h*F1_improved # Utilizando o metodo de Euler pra calcular o Yn+1
			F2_improved = expression.subs([(t, tEulerImproved+h), (y, Yh)]) 

			yEulerImproved += h*((F1_improved+F2_improved)/2)
			tEulerImproved += h
			Y.append(yEulerImproved)
			output.write(str(i) + " " + str(yEulerImproved) + "\n")

		T = t0 + (degree - 1) * h			

		for i in range(degree - 1, steps): 
			if degree == 1:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F1 = (1/2)*expression.subs([(t, Tn1), (y, Yn1)])
				F0 = (1/2)*expression.subs([(t, T), (y, Y[i])])

				yni = Y[i] + h*(F1 + F0)
				yni += (-1/12)
			elif degree == 2:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
				F1 = (2/3)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/12)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F2 + F1 + F0)
				yni += (1/24)
			elif degree == 3:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
				F2 = (19/24)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-5/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 4:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
				F3 = (323/360)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-11/30)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])

				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 5:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F5 = (95/288)*expression.subs([(t, Tn1), (y, Yn1)])
				F4 = (1427/1440)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-133/240)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-863/60480)
			elif degree == 6:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F6 = (19087/60480)*expression.subs([(t, Tn1), (y, Yn1)])
				F5 = (2713/2520)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-15487/20160)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-6737/20160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (275/24192)
			elif degree == 7:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F7 = (5257/17280)*expression.subs([(t, Tn1), (y, Yn1)])
				F6 = (139849/120960)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-4511/4480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-88547/120960)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-33953/3628800)

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()# Pronto \

	def adamMultonRungeKutta(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo de Adam-Multon por Runge Kutta de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")

		yRungeKutta = y0
		tRungeKutta = t0

		for i in range(1, degree):
			K1 = expression.subs([(t, tRungeKutta), (y, yRungeKutta)])

			T2 = tRungeKutta + (h/2)
			Y2 = yRungeKutta+(h/2)*K1
			K2 = expression.subs([(t, T2), (y, Y2)])

			T3 = tRungeKutta + (h/2)
			Y3 = yRungeKutta+(h/2)*K2
			K3 = expression.subs([(t, T3), (y, Y3)])
			
			T4 = tRungeKutta + h
			Y4 = yRungeKutta + h*K3
			K4 = expression.subs([(t, T4), (y, Y4)])

			K = (K1+(2*K2)+(2*K3)+K4)/6

			yRungeKutta += K*h
			tRungeKutta += h
			Y.append(yRungeKutta)
			output.write(str(i) + " " + str(yRungeKutta) + "\n")

		T = t0 + (degree - 1) * h			

		for i in range(degree - 1, steps): 
			if degree == 1:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F1 = (1/2)*expression.subs([(t, Tn1), (y, Yn1)])
				F0 = (1/2)*expression.subs([(t, T), (y, Y[i])])

				yni = Y[i] + h*(F1 + F0)
				yni += (-1/12)
			elif degree == 2:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
				F1 = (2/3)*expression.subs([(t, T), (y, Y[i])])
				F0 = (-1/12)*expression.subs([(t, T-h), (y, Y[i-1])])

				yni = Y[i] + h*(F2 + F1 + F0)
				yni += (1/24)
			elif degree == 3:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
				F2 = (19/24)*expression.subs([(t, T), (y, Y[i])])
				F1 = (-5/24)*expression.subs([(t, T-h), (y, Y[i-1])])
				F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])

				yni = Y[i] + h*(F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 4:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
				F3 = (323/360)*expression.subs([(t, T), (y, Y[i])])
				F2 = (-11/30)*expression.subs([(t, T-h), (y, Y[i-1])])
				F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])

				yni = Y[i] + h*(F4 + F3 + F2 + F1 + F0)
				yni += (-19/720)
			elif degree == 5:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F5 = (95/288)*expression.subs([(t, Tn1), (y, Yn1)])
				F4 = (1427/1440)*expression.subs([(t, T), (y, Y[i])])
				F3 = (-133/240)*expression.subs([(t, T-h), (y, Y[i-1])])
				F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				
				yni = Y[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-863/60480)
			elif degree == 6:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F6 = (19087/60480)*expression.subs([(t, Tn1), (y, Yn1)])
				F5 = (2713/2520)*expression.subs([(t, T), (y, Y[i])])
				F4 = (-15487/20160)*expression.subs([(t, T-h), (y, Y[i-1])])
				F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F2 = (-6737/20160)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				
				yni = Y[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (275/24192)
			elif degree == 7:
				Tn1 = T + h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				F7 = (5257/17280)*expression.subs([(t, Tn1), (y, Yn1)])
				F6 = (139849/120960)*expression.subs([(t, T), (y, Y[i])])
				F5 = (-4511/4480)*expression.subs([(t, T-h), (y, Y[i-1])])
				F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, Y[i-2])])
				F3 = (-88547/120960)*expression.subs([(t, T-(3*h)), (y, Y[i-3])])
				F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, Y[i-4])])
				F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, Y[i-5])])
				F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, Y[i-6])])

				yni = Y[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
				yni += (-33953/3628800)

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()
	
	def inverseFormula(entrie):
		degree = int(entrie.pop())
		Y = []
		
		t0 = float(entrie[degree])
		h = float(entrie[degree + 1])
		steps = int(entrie[degree + 2])
		expression = entrie[degree + 3]

		output = open('saida.txt', 'a')
		output.write("Metodo da Formula Inversa de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		
		for i in range(1, degree):
			Y.append(float(entrie[i]))
			output.write(str(i-1) + " " + entrie[i] + "\n")
		expression = sympify(expression)

		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps): 
			if degree == 2:
				Tn1 = T + 2*h # Equivalente ao tn+1
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) # Equivalente ao yn+1

				yni = (4/3)*Y[i] + (-1/3)*Y[i-1] + (2/3)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 3:
				Tn1 = T + 3*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (18/11)*Y[i] + (-9/11)*Y[i-1] + (2/11)*Y[i-2] + (6/11)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 4:
				Tn1 = T + 4*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (48/25)*Y[i] + (-36/25)*Y[i-1] + (16/25)*Y[i-2] + (-3/25)*Y[i-3] + (12/25)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 5:
				Tn1 = T + 5*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (300/137)*Y[i] + (-300/137)*Y[i-1] + (200/137)*Y[i-2] + (-75/137)*Y[i-3] + (-75/137)*Y[i-4] + (12/137)*Y[i-5] + (60/137)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 6:
				Tn1 = T + 6*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (300/147)*Y[i] + (-450/147)*Y[i-1] + (400/147)*Y[i-2] + (-225/147)*Y[i-3] + (72/147)*Y[i-4] + (-10/147)*Y[i-5] + (60/147)*h*expression.subs([(t, Tn1), (y, Yn1)])

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def inverseFormulaEuler(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo da Formula Inversa por Euler de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")
		
		yEuler = y0
		tEuler = t0

		for i in range(1, degree):
			F = h*expression.subs([(t, tEuler), (y, yEuler)])

			yEuler += F
			tEuler += h
			Y.append(yEuler)
			output.write(str(i) + " " + str(yEuler) + "\n")


		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps): 
			if degree == 2:
				Tn1 = T + 2*h # Equivalente ao tn+1
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) # Equivalente ao yn+1

				yni = (4/3)*Y[i] + (-1/3)*Y[i-1] + (2/3)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 3:
				Tn1 = T + 3*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (18/11)*Y[i] + (-9/11)*Y[i-1] + (2/11)*Y[i-2] + (6/11)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 4:
				Tn1 = T + 4*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (48/25)*Y[i] + (-36/25)*Y[i-1] + (16/25)*Y[i-2] + (-3/25)*Y[i-3] + (12/25)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 5:
				Tn1 = T + 5*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (300/137)*Y[i] + (-300/137)*Y[i-1] + (200/137)*Y[i-2] + (-75/137)*Y[i-3] + (-75/137)*Y[i-4] + (12/137)*Y[i-5] + (60/137)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 6:
				Tn1 = T + 6*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (300/147)*Y[i] + (-450/147)*Y[i-1] + (400/147)*Y[i-2] + (-225/147)*Y[i-3] + (72/147)*Y[i-4] + (-10/147)*Y[i-5] + (60/147)*h*expression.subs([(t, Tn1), (y, Yn1)])

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def inverseFormulaEulerInverse(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo da Formula Inversa por Euler Inverso de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")
		
		yEulerInverse = y0
		tEulerInverse = t0

		for i in range(1, degree):
			Tn1 = tEulerInverse + h # Equivalente ao tn+1
			Yn1 = yEulerInverse + h*expression.subs([(t, tEulerInverse), (y, yEulerInverse)]) # Equivalente ao yn+1

			yEulerInverse += h*expression.subs([(t, Tn1), (y, Yn1)])
			tEulerInverse += h
			Y.append(yEulerInverse)
			output.write(str(i) + " " + str(yEulerInverse) + "\n")

		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps): 
			if degree == 2:
				Tn1 = T + 2*h # Equivalente ao tn+1
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) # Equivalente ao yn+1

				yni = (4/3)*Y[i] + (-1/3)*Y[i-1] + (2/3)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 3:
				Tn1 = T + 3*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (18/11)*Y[i] + (-9/11)*Y[i-1] + (2/11)*Y[i-2] + (6/11)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 4:
				Tn1 = T + 4*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (48/25)*Y[i] + (-36/25)*Y[i-1] + (16/25)*Y[i-2] + (-3/25)*Y[i-3] + (12/25)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 5:
				Tn1 = T + 5*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (300/137)*Y[i] + (-300/137)*Y[i-1] + (200/137)*Y[i-2] + (-75/137)*Y[i-3] + (-75/137)*Y[i-4] + (12/137)*Y[i-5] + (60/137)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 6:
				Tn1 = T + 6*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (300/147)*Y[i] + (-450/147)*Y[i-1] + (400/147)*Y[i-2] + (-225/147)*Y[i-3] + (72/147)*Y[i-4] + (-10/147)*Y[i-5] + (60/147)*h*expression.subs([(t, Tn1), (y, Yn1)])

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def inverseFormulaEulerImproved(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo da Formula Inversa por Euler Aprimorado de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")
		
		yEulerImproved = y0
		tEulerImproved = t0

		for i in range(1, degree):
			F1_improved = expression.subs([(t, tEulerImproved), (y, yEulerImproved)])
			Yh = yEulerImproved + h*F1_improved # Utilizando o metodo de Euler pra calcular o Yn+1
			F2_improved = expression.subs([(t, tEulerImproved+h), (y, Yh)]) 

			yEulerImproved += h*((F1_improved+F2_improved)/2)
			tEulerImproved += h
			Y.append(yEulerImproved)
			output.write(str(i) + " " + str(yEulerImproved) + "\n")

		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps): 
			if degree == 2:
				Tn1 = T + 2*h # Equivalente ao tn+1
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) # Equivalente ao yn+1

				yni = (4/3)*Y[i] + (-1/3)*Y[i-1] + (2/3)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 3:
				Tn1 = T + 3*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (18/11)*Y[i] + (-9/11)*Y[i-1] + (2/11)*Y[i-2] + (6/11)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 4:
				Tn1 = T + 4*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (48/25)*Y[i] + (-36/25)*Y[i-1] + (16/25)*Y[i-2] + (-3/25)*Y[i-3] + (12/25)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 5:
				Tn1 = T + 5*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (300/137)*Y[i] + (-300/137)*Y[i-1] + (200/137)*Y[i-2] + (-75/137)*Y[i-3] + (-75/137)*Y[i-4] + (12/137)*Y[i-5] + (60/137)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 6:
				Tn1 = T + 6*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (300/147)*Y[i] + (-450/147)*Y[i-1] + (400/147)*Y[i-2] + (-225/147)*Y[i-3] + (72/147)*Y[i-4] + (-10/147)*Y[i-5] + (60/147)*h*expression.subs([(t, Tn1), (y, Yn1)])

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

	def inverseFormulaRungeKutta(entrie):
		Y = []
		y0 = float(entrie[1])
		Y.append(y0)
		t0 = float(entrie[2])
		h = float(entrie[3])
		steps = int(entrie[4])
		expression = entrie[5]
		expression = sympify(expression)
		degree = int(entrie[6])

		output = open('saida.txt', 'a')
		output.write("Metodo da Formula Inversa por Runge-Kutta de ordem " + str(degree) + "\n")
		output.write("y( " + str(t0) + " ) = " + str(entrie[1]) + "\n")
		output.write("h = " + str(h) + "\n")
		output.write("0 " + str(y0) + "\n")
		
		yRungeKutta = y0
		tRungeKutta = t0

		for i in range(1, degree):
			K1 = expression.subs([(t, tRungeKutta), (y, yRungeKutta)])

			T2 = tRungeKutta + (h/2)
			Y2 = yRungeKutta+(h/2)*K1
			K2 = expression.subs([(t, T2), (y, Y2)])

			T3 = tRungeKutta + (h/2)
			Y3 = yRungeKutta+(h/2)*K2
			K3 = expression.subs([(t, T3), (y, Y3)])
			
			T4 = tRungeKutta + h
			Y4 = yRungeKutta + h*K3
			K4 = expression.subs([(t, T4), (y, Y4)])

			K = (K1+(2*K2)+(2*K3)+K4)/6

			yRungeKutta += K*h
			tRungeKutta += h
			Y.append(yRungeKutta)
			output.write(str(i) + " " + str(yRungeKutta) + "\n")

		T = t0 + (degree - 1) * h		

		for i in range(degree - 1, steps): 
			if degree == 2:
				Tn1 = T + 2*h # Equivalente ao tn+1
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) # Equivalente ao yn+1

				yni = (4/3)*Y[i] + (-1/3)*Y[i-1] + (2/3)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 3:
				Tn1 = T + 3*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (18/11)*Y[i] + (-9/11)*Y[i-1] + (2/11)*Y[i-2] + (6/11)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 4:
				Tn1 = T + 4*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (48/25)*Y[i] + (-36/25)*Y[i-1] + (16/25)*Y[i-2] + (-3/25)*Y[i-3] + (12/25)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 5:
				Tn1 = T + 5*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])]) 

				yni = (300/137)*Y[i] + (-300/137)*Y[i-1] + (200/137)*Y[i-2] + (-75/137)*Y[i-3] + (-75/137)*Y[i-4] + (12/137)*Y[i-5] + (60/137)*h*expression.subs([(t, Tn1), (y, Yn1)])
			elif degree == 6:
				Tn1 = T + 6*h 
				Yn1 = Y[i] + h*expression.subs([(t, T), (y, Y[i])])

				yni = (300/147)*Y[i] + (-450/147)*Y[i-1] + (400/147)*Y[i-2] + (-225/147)*Y[i-3] + (72/147)*Y[i-4] + (-10/147)*Y[i-5] + (60/147)*h*expression.subs([(t, Tn1), (y, Yn1)])

			output.write(str(i+1) + " " + str(yni) + "\n")
			Y.append(yni)
			T += h
		output.write("\n")
		output.close()

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
	elif strct[0] == "adam_bashforth_by_euler":
		methods.adamBashforthEuler(strct)
	elif strct[0] == "adam_bashforth_by_euler_aprimorado":
		methods.adamBashforthEulerImproved(strct)
	elif strct[0] == "adam_bashforth_by_euler_inverso":
		methods.adamBashforthEulerInverse(strct)
	elif strct[0] == "adam_bashforth_by_runge_kutta":
		methods.adamBashforthRungeKutta(strct)
	elif strct[0] == "adam_multon":
		methods.adamMulton(strct)
	elif strct[0] == "adam_multon_by_euler":
		methods.adamMultonEuler(strct)
	elif strct[0] == "adam_multon_by_euler_inverso":
		methods.adamMultonEulerInverse(strct)
	elif strct[0] == "adam_multon_by_euler_aprimorado":
		methods.adamMultonEulerImproved(strct)
	elif strct[0] == "adam_multon_by_runge_kutta":
		methods.adamMultonRungeKutta(strct)
	elif strct[0] == "formula_inversa":
		methods.inverseFormula(strct)
	elif strct[0] == "formula_inversa_by_euler":
		methods.inverseFormulaEuler(strct)
	elif strct[0] == "formula_inversa_by_euler_inverso":
		methods.inverseFormulaEulerInverse(strct)
	elif strct[0] == "formula_inversa_by_euler_aprimorado":
		methods.inverseFormulaEulerImproved(strct)
	elif strct[0] == "formula_inversa_by_runge_kutta":
		methods.inverseFormulaRungeKutta(strct)

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
