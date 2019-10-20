import mpmath
from sympy.plotting import plot as symplot
from sympy import *	
import matplotlib
import matplotlib.pyplot as plt
x, y, z, t = symbols("x y z t")

def euler(y0, t0, h, step, exp, writeFlag):
	xeuler = []
	yeuler = []

	if writeFlag:
		printHeader("Metodo de Euler simples", y0, t0, h)

	expression = sympify(exp)

	Y = y0
	T = t0

	yeuler.append(y0)
	xeuler.append(t0)

	for i in range(step + 1):
		# Yn+1 = Yn + h*fn
		Y = Y + h*expression.subs([(t, T), (y, Y)])
		T = T + h

		yeuler.append(Y)
		xeuler.append(T)

		if writeFlag:
			outputAppend(str(i) + " " + str(yeuler[i]) + "\n")
	if writeFlag:
		outputAppend("\n")

	return xeuler, yeuler

def eulerInverse(y0, t0, h, step, exp, writeFlag):
	xeulerInverso = []
	yeulerInverso = []

	if writeFlag:
		printHeader("Metodo de Euler Inverso", y0, t0, h)

	expression = sympify(exp)

	Y = y0
	T = t0

	yeulerInverso.append(y0)
	xeulerInverso.append(t0)

	for i in range(step+1):
		# Yn+1 = Yn + h *F(tn+1, yn+1)
		# Previsoes 
		Tn1 = T + h # Equivalente ao tn+1
		Yn1 = Y + h*expression.subs([(t, T), (y, Y)]) # Equivalente ao yn+1

		Y += h*expression.subs([(t, Tn1), (y, Yn1)])
		T += h

		yeulerInverso.append(Y)
		xeulerInverso.append(T)
		if writeFlag:
			outputAppend(str(i) + " " + str(yeulerInverso[i]) + "\n")
	if writeFlag:
		outputAppend("\n")

	return xeulerInverso, yeulerInverso

def eulerImproved(y0, t0, h, step, exp, writeFlag):
	xeulerImproved = []
	yeulerImproved = []

	if writeFlag:
		printHeader("Metodo de Euler Aprimorado", y0, t0, h)

	expression = sympify(exp)
	
	yeulerImproved.append(y0)
	xeulerImproved.append(t0)

	Y = y0
	T = t0

	for i in range(step):
		# Y = Y * h((F1+F2)/2)

		F1 = expression.subs([(t, T), (y, Y)])
		Yh = Y + h*F1 # Utilizando o metodo de Euler pra calcular o Yn+1
		F2 = expression.subs([(t, T+h), (y, Yh)]) 

		Y += h*((F1+F2)/2)
		T += h

		xeulerImproved.append(T)
		yeulerImproved.append(Y)

		if writeFlag:
			outputAppend(str(i+1) + " " + str(Y) + "\n")

	if writeFlag:
		outputAppend("\n")

	return xeulerImproved, yeulerImproved

def rungeKutta(y0, t0, h, step, exp, writeFlag):
	xrungeKutta = []
	yrungeKutta = []

	if writeFlag:
		printHeader("Metodo de Runge-Kutta", y0, t0, h)

	expression = sympify(exp)
	
	yrungeKutta.append(y0)
	xrungeKutta.append(t0)

	Y = y0
	T = t0

	for i in range(1, (step + 1)):
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

		yrungeKutta.append(Y)
		xrungeKutta.append(T)

		if writeFlag:
			outputAppend(str(i) + " " + str(Y) + "\n")
	if writeFlag:
		outputAppend("\n")

	return xrungeKutta, yrungeKutta

def adamBashforth(y0, t0, h, step, exp, degree, method, flag, flagPrint):
	xAdam = []
	yAdam = []

	if flagPrint:
		printHeader(method, y0[0], t0[0], h)

	expression = sympify(exp)
	T = t0[0]

	for i in range(degree):
		if flagPrint:
			outputAppend(str(i) + " " + str(y0[i]) + "\n")
		xAdam.append(T)
		yAdam.append(y0[i])
		if flag:
			T = t0[i]
		else:
			T = T + h
	T = t0[0] + (degree - 1) * h

	for i in range(degree - 1, step):
		if degree == 2: 
			F1 = (3/2)*expression.subs([(t, T), (y, yAdam[i])]) 
			F0 = (-1/2)*expression.subs([(t, T-h), (y, yAdam[i-1])]) 

			yni = yAdam[i] + h*(F1 + F0)
		elif degree == 3: 
			F2 = (23/12)*expression.subs([(t, T), (y, yAdam[i])])
			F1 = (-4/3)*expression.subs([(t, T-h), (y, yAdam[i-1])]) 
			F0 = (5/12)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])]) 

			yni = yAdam[i] + h*(F2 + F1 + F0)
		elif degree == 4: 
			F3 = (55/24)*expression.subs([(t, T), (y, yAdam[i])]) 
			F2 = (-59/24)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F1 = (37/24)*expression.subs([(t, T-(3*h)), (y, yAdam[i-2])]) 
			F0 = (-3/8)*expression.subs([(t, T-(4*h)), (y, yAdam[i-3])]) 

			yni = yAdam[i] + h*(F3 + F2 + F1 + F0)
		elif degree == 5:
			F4 = (1901/720)*expression.subs([(t, xAdam[i]), (y,yAdam[i])])
			F3 = (-1387/360)*expression.subs([(t, xAdam[i-1]), (y, yAdam[i-1])]) 
			F2 = (109/30)*expression.subs([(t, xAdam[i-2]), (y, yAdam[i-2])])
			F1 = (-637/360)*expression.subs([(t, xAdam[i-3]), (y, yAdam[i-3])]) 
			F0 = (251/720)*expression.subs([(t, xAdam[i-4]), (y, yAdam[i-4])]) 
			
			yni = yAdam[i] + h*(F4 + F3 + F2 + F1 + F0)
		elif degree == 6: 
			F5 = (4277/1440)*expression.subs([(t,T), (y, yAdam[i])]) 
			F4 = (-2641/480)*expression.subs([(t, T-h), (y, yAdam[i-1])]) 
			F3 = (4991/720)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])]) 
			F2 = (-3649/720)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F1 = (959/480)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])]) 
			F0 = (-95/288)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])]) 
			
			yni = yAdam[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
		elif degree == 7:
			F6 = (198721/60480)*expression.subs([(t, T), (y, yAdam[i])]) 
			F5 = (-18637/2520)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F4 = (235183/20160)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])  
			F3 = (-10754/945)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])]) 
			F2 = (135713/20160)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			F1 = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])]) 
			F0 = (19087/60480)*expression.subs([(t, T-(6*h)), (y, yAdam[i-6])]) 

			yni = yAdam[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
		elif degree == 8:
			F7 = (16083/4480) * expression.subs([(t,T), (y, yAdam[i])])
			F6 = (-1152169/120960) * expression.subs([(t, T-h), (y, yAdam[i-1])]) 
			F5 = (242653/13440) * expression.subs([(t, T-(2*h)), (y, yAdam[i-2])]) 
			F4 = (-296053/13440) * expression.subs([(t, T-(3*h)), (y, yAdam[i-3])]) 
			F3 = (2102243/120960) * expression.subs([(t, T-(4*h)), (y, yAdam[i-4])]) 
			F2 = (-115747/13440) * expression.subs([(t, T-(5*h)), (y, yAdam[i-5])])
			F1 = (32863/13440) * expression.subs([(t, T-(6*h)), (y, yAdam[i-6])]) 
			F0 = (-5257/17280) * expression.subs([(t, T-(7*h)), (y, yAdam[i-7])]) 

			yni = yAdam[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
		if flagPrint:
			outputAppend(str(i+1) + " " + str(yni) + "\n")
		yAdam.append(yni)
		T += h
		xAdam.append(T)
		
	if flagPrint:
		outputAppend("\n")
	return xAdam, yAdam

def adamMulton(y0, t0, h, step, exp, degree, method, flag, flagPrint):
	xAdam = []
	yAdam = []

	if flagPrint:
		printHeader(method, y0[0], t0[0], h)

	expression = sympify(exp)
	T = t0[0]
	
	for i in range(degree):
		if flagPrint:
			outputAppend(str(i) + " " + str(y0[i]) + "\n")
		xAdam.append(T)
		yAdam.append(y0[i])
		if flag:
			T = t0[i]
		else:
			T = T + h
	T = t0[0] + degree * h

	for i in range(degree - 1, step): 
		if degree == 2: #yn+1 = yn + (h/2)*(f(tn+1, yn+1) + f(tn, yn))

			F1B = (3/2)*expression.subs([(t, T), (y, yAdam[i])]) 
			F0B = (-1/2)*expression.subs([(t, T-h), (y, yAdam[i-1])]) 

			Yn1 = yAdam[i] + h*(F1B + F0B)
			# --------------------------------------------------------
			Tn1 = T + h 
			F1 = expression.subs([(t, Tn1), (y, Yn1)])
			F0 = expression.subs([(t, T), (y, yAdam[i])])

			yni = yAdam[i] + (h/2)*(F1 + F0)
		elif degree == 3: # yn+1 = yn + h*((5/12)*f(tn+1, yn+1) + (2/3)*f(tn, yn) - (1/12)*f(tn-1, yn-1))
			F2B = (23/12)*expression.subs([(t, T), (y, yAdam[i])])
			F1B = (-4/3)*expression.subs([(t, T-h), (y, yAdam[i-1])]) 
			F0B = (5/12)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])]) 

			Yn1 = yAdam[i] + h*(F2B + F1B + F0B)
			#--------------------------------------------------------
			Tn1 = T + h 
			F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
			F1 = (2/3)*expression.subs([(t, T), (y, yAdam[i])])
			F0 = (-1/12)*expression.subs([(t, T-h), (y, yAdam[i-1])])

			yni = yAdam[i] + h*(F2 + F1 + F0)
		elif degree == 4: # yn+1 = yn + h*((3/8)*f(tn+1, yn+1) + (19/24)*f(tn, yn) - (5/24)*f(tn-1, yn-1) + (1/24)*f(tn-2, yn-2))
			F3B = (55/24)*expression.subs([(t, T), (y, yAdam[i])]) 
			F2B = (-59/24)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F1B = (37/24)*expression.subs([(t, T-(3*h)), (y, yAdam[i-2])]) 
			F0B = (-3/8)*expression.subs([(t, T-(4*h)), (y, yAdam[i-3])]) 

			Yn1 = yAdam[i] + h*(F3B + F2B + F1B + F0B)
			#-------------------------------------------------------
			Tn1 = T + h 
			F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
			F2 = (19/24)*expression.subs([(t, T), (y, yAdam[i])])
			F1 = (-5/24)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])

			yni = yAdam[i] + h*(F3 + F2 + F1 + F0)
		elif degree == 5:# yn+1 = yn + h*((251/720)*f(tn+1, yn+1) + (323/360)*f(tn, yn) - (11/30)*f(tn-1, yn-1) + (53/360)*f(tn-2, yn-2) - (19/720)*f(tn-3, yn-3))
			F4B = (1901/720)*expression.subs([(t, xAdam[i]), (y,yAdam[i])])
			F3B = (-1387/360)*expression.subs([(t, xAdam[i-1]), (y, yAdam[i-1])]) 
			F2B = (109/30)*expression.subs([(t, xAdam[i-2]), (y, yAdam[i-2])])
			F1B = (-637/360)*expression.subs([(t, xAdam[i-3]), (y, yAdam[i-3])]) 
			F0B = (251/720)*expression.subs([(t, xAdam[i-4]), (y, yAdam[i-4])]) 
			
			Yn1 = yAdam[i] + h*(F4B + F3B + F2B + F1B + F0B)
			#------------------------------------------------------
			Tn1 = T + h 

			F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
			F3 = (323/360)*expression.subs([(t, T), (y, yAdam[i])])
			F2 = (-11/30)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])

			yni = yAdam[i] + h*(F4 + F3 + F2 + F1 + F0)
		elif degree == 6: # yn+1 = yn + h*((95/288)*f(tn+1, yn+1) + (1427/1440)*f(tn, yn) - (133/240)*f(tn-1, yn-1) + (241/720)*f(tn-2, yn-2) - (173/1440)*f(tn-3, yn-3) + (3/160)*f(tn-4, yn-4))
			F5B = (4277/1440)*expression.subs([(t,T), (y, yAdam[i])]) 
			F4B = (-2641/480)*expression.subs([(t, T-h), (y, yAdam[i-1])]) 
			F3B = (4991/720)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])]) 
			F2B = (-3649/720)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F1B = (959/480)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])]) 
			F0B = (-95/288)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])]) 
			
			Yn1 = yAdam[i] + h*(F5B + F4B + F3B + F2B + F1B + F0B)

			#--------------------------------------------------------

			F5 = (95/288)*expression.subs([(t, T+h), (y, Yn1)])
			F4 = (1427/1440)*expression.subs([(t, T), (y, yAdam[i])])
			F3 = (-133/240)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			
			yni = yAdam[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			
		elif degree == 7: # yn+1 = yn + h*((19087/60480)*f(tn+1, yn+1) + (2713/2520)*f(tn, yn) - (15487/20160)*f(tn-1, yn-1) + (586/945)*f(tn-2, yn-2) - (5737/20160)*f(tn-3, yn-3) + (263/2520)*f(tn-4, yn-4) - (863/60480)*f(tn-5, yn-5))
			F6B = (198721/60480)*expression.subs([(t, T), (y, yAdam[i])]) 
			F5B = (-18637/2520)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F4B = (235183/20160)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])  
			F3B = (-10754/945)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])]) 
			F2B = (135713/20160)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			F1B = (-5603/2520)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])]) 
			F0B = (19087/60480)*expression.subs([(t, T-(6*h)), (y, yAdam[i-6])]) 

			Yn1 = yAdam[i] + h*(F6B + F5B + F4B + F3B + F2B + F1B + F0B)
			#---------------------------------------------------------

			F6 = (19087/60480)*expression.subs([(t, T+h), (y, Yn1)])
			F5 = (2713/2520)*expression.subs([(t, T), (y, yAdam[i])])
			F4 = (-15487/20160)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F2 = (-5737/20160)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])])
			
			yni = yAdam[i] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
		elif degree == 8: # yn+1 = yn + h*((5257/17280)*f(tn+1, yn+1) + (139849/120960)*f(tn, yn) - (4511/4480)*f(tn-1, yn-1) + (123133/120960)*f(tn-2, yn-2) - (88574/120960)*f(tn-3, yn-3) + (1537/4480)*f(tn-4, yn-4) - (11351/120960)*f(tn-5, yn-5) + (275/24192)*f(tn-6, yn-6))
			F7B = (16083/4480) * expression.subs([(t,T), (y, yAdam[i])])
			F6B = (-1152169/120960) * expression.subs([(t, T-h), (y, yAdam[i-1])]) 
			F5B = (242653/13440) * expression.subs([(t, T-(2*h)), (y, yAdam[i-2])]) 
			F4B = (-296053/13440) * expression.subs([(t, T-(3*h)), (y, yAdam[i-3])]) 
			F3B = (2102243/120960) * expression.subs([(t, T-(4*h)), (y, yAdam[i-4])]) 
			F2B = (-115747/13440) * expression.subs([(t, T-(5*h)), (y, yAdam[i-5])])
			F1B = (32863/13440) * expression.subs([(t, T-(6*h)), (y, yAdam[i-6])]) 
			F0B = (-5257/17280) * expression.subs([(t, T-(7*h)), (y, yAdam[i-7])]) 

			Yn1 = yAdam[i] + h*(F7B + F6B + F5B + F4B + F3B + F2B + F1B + F0B)
			#-------------------------------------------------------
			F7 = (5257/17280)*expression.subs([(t, T+h), (y, Yn1)])
			F6 = (139849/120960)*expression.subs([(t, T), (y, yAdam[i])])
			F5 = (-4511/4480)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F3 = (-88574/120960)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])])
			F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, yAdam[i-6])])

			yni = yAdam[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
		if flagPrint:
			outputAppend(str(i+1) + " " + str(yni) + "\n")
		yAdam.append(yni)
		T += h
		xAdam.append(T)
	if flagPrint:
		outputAppend("\n")
	return xAdam, yAdam

def inverseFormula(y0, t0, h, step, exp, degree, method, flag, flagPrint):
	xInverse = []
	yInverse = []

	if flagPrint:
		printHeader(method, y0[0], t0[0], h)

	expression = sympify(exp)
	T = t0[0]

	for i in range(degree):
		if flagPrint:
			outputAppend(str(i) + " " + str(y0[i]) + "\n")
		xInverse.append(T)
		yInverse.append(y0[i])
		if flag:
			T = t0[i]
		else:
			T = T + h
	T = t0[0] + degree * h

	for i in range(degree-1, step):
		if(degree == 2): # yn+1 = (4/3)*yn - (1/3)*yn-1 + (2/3)*h*f(tn+1, yn+1)
			F1B = (3/2)*expression.subs([(t, T), (y, yInverse[i])]) 
			F0B = (-1/2)*expression.subs([(t, T-h), (y, yInverse[i-1])]) 

			Yn1 = yInverse[i] + h*(F1B + F0B)
			#---------------------------------------------------------------

			yni = (4/3)*yInverse[i-1] - (1/3)*yInverse[i-2] + (2/3)*h*expression.subs([(t, T+h), (y, Yn1)])

		if(degree == 3): # yn+1 = (18/11)*yn - (9/11)*yn-1 + (2/11)*yn-2 + (6/11)*h*f(tn+1, yn+1)
			F2B = (23/12)*expression.subs([(t, T), (y, yInverse[i])])
			F1B = (-4/3)*expression.subs([(t, T-h), (y, yInverse[i-1])]) 
			F0B = (5/12)*expression.subs([(t, T-(2*h)), (y, yInverse[i-2])]) 

			Yn1 = yInverse[i] + h*(F2B + F1B + F0B)
			#---------------------------------------------------------------

			yni = (18/11)*yInverse[i] - (9/11)*yInverse[i-1] + (2/11)*yInverse[i-2] + (6/11)*h*expression.subs([(t, T+h), (y, Yn1)])

		if(degree == 4): # yn+1 = (48/25)*yn - (36/25)*yn-1 + (16/25)*yn-2 - (3/25)*yn-3 + (12/25)*h*f(tn+1, yn+1)
			F3B = (55/24)*expression.subs([(t, T), (y, yInverse[i])]) 
			F2B = (-59/24)*expression.subs([(t, T-h), (y, yInverse[i-1])])
			F1B = (37/24)*expression.subs([(t, T-(3*h)), (y, yInverse[i-2])]) 
			F0B = (-3/8)*expression.subs([(t, T-(4*h)), (y, yInverse[i-3])]) 

			Yn1 = yInverse[i] + h*(F3B + F2B + F1B + F0B)			
			#---------------------------------------------------------------

			yni = (48/25)*yInverse[i] - (36/25)*yInverse[i-1] + (16/25)*yInverse[i-2] - (3/25)*yInverse[i-3] + (12/25)*h*expression.subs([(t,T+h), (y, Yn1)])

		if(degree == 5): # yn+1 = (300/137)*yn - (300/137)*yn-1 + (200/137)*yn-2 - (75/137)*yn-3 + (12/137)*yn-4 + (60/137)*h*f(tn+1, yn+1)
			F4B = (1901/720)*expression.subs([(t, xInverse[i]), (y,yInverse[i])])
			F3B = (-1387/360)*expression.subs([(t, xInverse[i-1]), (y, yInverse[i-1])]) 
			F2B = (109/30)*expression.subs([(t, xInverse[i-2]), (y, yInverse[i-2])])
			F1B = (-637/360)*expression.subs([(t, xInverse[i-3]), (y, yInverse[i-3])]) 
			F0B = (251/720)*expression.subs([(t, xInverse[i-4]), (y, yInverse[i-4])]) 
			
			Yn1 = yInverse[i] + h*(F4B + F3B + F2B + F1B + F0B)
			#---------------------------------------------------------------

			yni = (300/137)*yInverse[i] - (300/137)*yInverse[i-1] + (200/137)*yInverse[i-2] - (75/137)*yInverse[i-3] + (12/137)*yInverse[i-4] + (60/137)*h*expression.subs([(t, T+h), (y, Yn1)])

		if(degree == 6): # yn+1 = (360/147)*yn - (450/147)*yn-1 + (400/147)*yn-2 - (225/147)*yn-3 + (72/147)*yn-4 - (10/147)*yn-5 + (60/147)*h*f(tn+1, yn+1)
			F5B = (4277/1440)*expression.subs([(t,T), (y, yInverse[i])]) 
			F4B = (-2641/480)*expression.subs([(t, T-h), (y, yInverse[i-1])]) 
			F3B = (4991/720)*expression.subs([(t, T-(2*h)), (y, yInverse[i-2])]) 
			F2B = (-3649/720)*expression.subs([(t, T-(3*h)), (y, yInverse[i-3])])
			F1B = (959/480)*expression.subs([(t, T-(4*h)), (y, yInverse[i-4])]) 
			F0B = (-95/288)*expression.subs([(t, T-(5*h)), (y, yInverse[i-5])]) 
			
			Yn1 = yInverse[i] + h*(F5B + F4B + F3B + F2B + F1B + F0B)
			#---------------------------------------------------------------

			yni = (360/147)*yInverse[i] - (450/147)*yInverse[i-1] + (400/147)*yInverse[i-2] - (225/147)*yInverse[i-3] + (72/147)*yInverse[i-4] - (10/147)*yInverse[i-5] + (60/147)*h*expression.subs([(t, T+h), (y, Yn1)])
		if flagPrint:
			outputAppend(str(i+1) + " " + str(yni) + "\n")
		yInverse.append(yni)
		T += h
		xInverse.append(T)
	if flagPrint:
		outputAppend("\n")
	return xInverse, yInverse

def outputAppend(text):
	output = open('saida.txt', 'a')
	output.write(text)
	output.close()

def plotGraphic(xPlot, yPlot, Name, Savename):
	fig, pltQ = plt.subplots()
	line = pltQ.plot(xPlot, yPlot, 'k-o', label=Name)
	pltQ.set(xlabel='T', ylabel='Y', title=Name)
	pltQ.grid()
	pltQ.tick_params(labelcolor='k', labelsize='small', width=1)
	pltQ.legend()

	Savename += ".png"
	fig.savefig(Savename)
	return

def printHeader(name, y0, t0, h):
	output = open('saida.txt', 'a')
	output.write(name + "\n")
	output.write("y( " + str(t0) + " ) = " + str(y0) + "\n")
	output.write("h = " + str(h) + "\n")
	output.close()
	return

def switchCase(strct):
	if strct[0] == "euler": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		plotxEuler, plotyEuler = euler(y0, t0, h, step, expression, True)
		plotGraphic(plotxEuler, plotyEuler, "Euler", "Euler")

	elif strct[0] == "euler_inverso": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		plotxEulerInverse, plotyEulerInverse = eulerInverse(y0, t0, h, step, expression, True)
		plotGraphic(plotxEulerInverse, plotyEulerInverse, "Euler Inverso", "Euler_inverso")
	
	elif strct[0] == "euler_aprimorado": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		plotxEulerImproved, plotyEulerImproved = eulerImproved(y0, t0, h, step, expression, True)
		plotGraphic(plotxEulerImproved, plotyEulerImproved, "Euler Aprimorado", "euler_aprimorado")
	
	elif strct[0] == "runge_kutta": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		plotxRungeKutta, plotyRungeKutta = rungeKutta(y0, t0, h, step, expression, True)
		plotGraphic(plotxRungeKutta, plotyRungeKutta, "Runge-Kutta", "runge_kutta")

	elif strct[0] == "adam_bashforth": # Done
		y0 = []
		t0 = []
		degree = int(strct.pop())

		for i in range(1,degree):
			y0.append(float(strct[i]))

		t0.append(float(strct[degree]))
		h =  float(strct[degree + 1])
		step = int(strct[degree + 2])
		expression = strct[degree + 3]

		methodTitle = "Metodo de Adam-Bashforth de Ordem " + str(degree)
		fileImageTitle = "bashforth" + str(degree) + "th"

		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0, t0, h, step, expression, degree-1, methodTitle, False, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_euler": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Euler de Ordem " + str(degree)
		fileImageTitle = "bashforthEuler" + str(degree) + "th"
		x0Euler, y0Euler = euler(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0Euler, x0Euler, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_euler_aprimorado": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Euler aprimorado de Ordem " + str(degree)
		fileImageTitle = "bashforthEulerImproved" + str(degree) + "th"
		x0EulerImproved, y0EulerImproved = eulerImproved(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0EulerImproved, x0EulerImproved, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_euler_inverso": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Euler Inverso de Ordem " + str(degree)
		fileImageTitle = "bashforthEulerInverse" + str(degree) + "th"
		x0EulerInverse, y0EulerInverse = eulerInverse(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0EulerInverse, x0EulerInverse, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_runge_kutta": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Runge-Kutta de Ordem " + str(degree)
		fileImageTitle = "bashforthRungeKutta" + str(degree) + "th"
		x0RungeKutta, y0RungeKutta = rungeKutta(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0RungeKutta, x0RungeKutta, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon": # Done
		y0 = []
		t0 = []
		degree = int(strct.pop())

		for i in range(1,degree):
			y0.append(float(strct[i]))

		t0.append(float(strct[degree]))
		h =  float(strct[degree + 1])
		step = int(strct[degree + 2])
		expression = strct[degree + 3]

		methodTitle = "Metodo de Adam-Multon de Ordem " + str(degree)
		fileImageTitle = "multon" + str(degree) + "th"

		plotxMulton, plotyMulton = adamMulton(y0, t0, h, step, expression, degree-1, methodTitle, False, True)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)
 
	elif strct[0] == "adam_multon_by_euler": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Multon por Euler de Ordem " + str(degree)
		fileImageTitle = "multonEuler" + str(degree) + "th"
		x0Euler, y0Euler = euler(y0, t0, h, step, expression, False)
		plotxMulton, plotyMulton = adamMulton(y0Euler, x0Euler, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon_by_euler_inverso": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Multon por Euler Inverso de Ordem " + str(degree)
		fileImageTitle = "multonEulerInverse" + str(degree) + "th"
		x0EulerInverse, y0EulerInverse = eulerInverse(y0, t0, h, step, expression, False)
		plotxMulton, plotyMulton = adamMulton(y0EulerInverse, x0EulerInverse, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon_by_euler_aprimorado": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Multon por Euler Aprimorado de Ordem " + str(degree)
		fileImageTitle = "multonEulerImproved" + str(degree) + "th"
		x0EulerImproved, y0EulerImproved = eulerImproved(y0, t0, h, step, expression, False)
		plotxMulton, plotyMulton = adamMulton(y0EulerImproved, x0EulerImproved, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon_by_runge_kutta": # Done
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Multon por Runge-Kutta de Ordem " + str(degree)
		fileImageTitle = "multonRungeKutta" + str(degree) + "th"
		x0RungeKutta, y0RungeKutta = rungeKutta(y0, t0, h, step, expression, False)
		plotxMulton, plotyMulton = adamMulton(y0RungeKutta, x0RungeKutta, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)

	elif strct[0] == "formula_inversa":
		y0 = []
		t0 = []
		degree = int(strct.pop())

		for i in range(1,degree):
			y0.append(float(strct[i]))

		t0.append(float(strct[degree]))
		h =  float(strct[degree + 1])
		step = int(strct[degree + 2])
		expression = strct[degree + 3]

		methodTitle = "Metodo de Formula Inversa de Ordem " + str(degree)
		fileImageTitle = "formulaInversa" + str(degree) + "th"

		plotxInverseFormula, plotyInverseFormula = inverseFormula(y0, t0, h, step, expression, degree-1, methodTitle, False, True)
		plotGraphic(plotxInverseFormula,plotyInverseFormula, methodTitle, fileImageTitle)

	elif strct[0] == "formula_inversa_by_euler":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Formula Inversa por Euler de Ordem " + str(degree)
		fileImageTitle = "formulaInversaEuler" + str(degree) + "th"
		x0Euler, y0Euler = euler(y0, t0, h, step, expression, False)
		plotxInverseFormula, plotyInverseFormula = inverseFormula(y0Euler, x0Euler, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxInverseFormula,plotyInverseFormula, methodTitle, fileImageTitle)

	elif strct[0] == "formula_inversa_by_euler_inverso":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Formula Inversa por Euler Inverso de Ordem " + str(degree)
		fileImageTitle = "formulaInversaEulerInverse" + str(degree) + "th"
		x0EulerInverse, y0EulerInverse = eulerInverse(y0, t0, h, step, expression, False)
		plotxInverseFormula, plotyInverseFormula = inverseFormula(y0EulerInverse, x0EulerInverse, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxInverseFormula,plotyInverseFormula, methodTitle, fileImageTitle)

	elif strct[0] == "formula_inversa_by_euler_aprimorado":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Formula Inversa por Euler Aprimorado de Ordem " + str(degree)
		fileImageTitle = "formulaInversaEulerImproved" + str(degree) + "th"
		x0EulerImproved, y0EulerImproved = eulerImproved(y0, t0, h, step, expression, False)
		plotxInverseFormula, plotyInverseFormula = inverseFormula(y0EulerImproved, x0EulerImproved, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxInverseFormula,plotyInverseFormula, methodTitle, fileImageTitle)

	elif strct[0] == "formula_inversa_by_runge_kutta":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Formula Inversa por Runge-Kutta de Ordem " + str(degree)
		fileImageTitle = "formulaInversaRungeKutta" + str(degree) + "th"
		x0RungeKutta, y0RungeKutta = rungeKutta(y0, t0, h, step, expression, False)
		plotxInverseFormula, plotyInverseFormula = inverseFormula(y0RungeKutta, x0RungeKutta, h, step, expression, degree, methodTitle, True, True)
		plotGraphic(plotxInverseFormula,plotyInverseFormula, methodTitle, fileImageTitle)
	
	elif strct[0] == "all": 
		fig, pltQ = plt.subplots()

		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		plotxEuler, plotyEuler = euler(y0, t0, h, step, expression, True)
		lineEuler = pltQ.plot(plotxEuler, plotyEuler, 'k-o', label="Euler")
		
		plotxEulerInverse, plotyEulerInverse = eulerInverse(y0, t0, h, step, expression, True)
		lineEulerInverse = pltQ.plot(plotxEulerInverse, plotyEulerInverse, 'b-o', label="Euler Inverso")

		plotxEulerImproved, plotyEulerImproved = eulerImproved(y0, t0, h, step, expression, True)
		lineEulerImproved = pltQ.plot(plotxEulerImproved, plotyEulerImproved, 'c-o', label="Euler Aprimorado")

		plotxRungeKutta, plotyRungeKutta = rungeKutta(y0, t0, h, step, expression, True)
		lineRungeKutta = pltQ.plot(plotxRungeKutta, plotyRungeKutta, 'm-o', label="Runge-Kutta")

		methodTitle = "Metodo de Adam-Bashforth Euler de Ordem " + str(degree)
		plotxBashforthEuler, plotyBashforthEuler = adamBashforth(plotyEuler, plotxEuler, h, step, expression, degree, methodTitle, True, True)
		lineBashforthEuler = pltQ.plot(plotxBashforthEuler, plotyBashforthEuler, marker='o', color='xkcd:grey', label="Bashforth Euler")

		methodTitle = "Metodo de Adam-Bashforth Euler Inverso de Ordem " + str(degree)
		plotxBashforthEulerInverse, plotyBashforthEulerInverse = adamBashforth(plotyEulerInverse, plotxEulerInverse, h, step, expression, degree, methodTitle, True, True)
		lineBashforthEulerInverse = pltQ.plot(plotxBashforthEulerInverse, plotyBashforthEulerInverse, marker='o', color='xkcd:beige', label="Bashforth Euler Inverso")

		methodTitle = "Metodo de Adam-Bashforth Euler Aprimorado de Ordem " + str(degree)
		plotxBashforthEulerImproved, plotyBashforthEulerImproved = adamBashforth(plotyEulerImproved, plotxEulerImproved, h, step, expression, degree, methodTitle, True, True)
		lineBashforthEulerImproved = pltQ.plot(plotxBashforthEulerImproved, plotyBashforthEulerImproved, marker='o', color='xkcd:chartreuse', label="Bashforth Euler Aprimorado")

		methodTitle = "Metodo de Adam-Bashforth Runge-Kutta de Ordem " + str(degree)
		plotxBashforthRungeKutta, plotyBashforthRungeKutta = adamBashforth(plotyRungeKutta, plotxRungeKutta, h, step, expression, degree, methodTitle, True, True)
		lineBashforthRungeKutta = pltQ.plot(plotxBashforthRungeKutta, plotyBashforthRungeKutta, marker='o', color='xkcd:brown', label="Bashforth Runge-Kutta")

		methodTitle = "Metodo de Adam-Multon Euler de Ordem " + str(degree)
		plotxMultonEuler, plotyMultonEuler = adamMulton(plotyEuler, plotxEuler, h, step, expression, degree, methodTitle, True, True)
		lineMultonEuler = pltQ.plot(plotxMultonEuler, plotyMultonEuler,  marker='o', color='xkcd:coral', label="Moulton Euler")

		methodTitle = "Metodo de Adam-Multon Euler Inverso de Ordem " + str(degree)
		plotxMultonEulerInverse, plotyMultonEulerInverse = adamMulton(plotyEulerInverse, plotxEulerInverse, h, step, expression, degree, methodTitle, True, True)
		lineMultonEulerInverse = pltQ.plot(plotxMultonEulerInverse, plotyMultonEulerInverse,  marker='o', color='xkcd:fuchsia', label="Moulton Euler Inverso")
	
		methodTitle = "Metodo de Adam-Multon Euler Aprimorado de Ordem " + str(degree)
		plotxMultonEulerImproved, plotyMultonEulerImproved = adamMulton(plotyEulerImproved, plotxEulerImproved, h, step, expression, degree, methodTitle, True, True)
		lineMultonEulerImproved = pltQ.plot(plotxMultonEulerImproved, plotyMultonEulerImproved,  marker='o', color='xkcd:gold', label="Moulton Euler Aprimorado")

		methodTitle = "Metodo de Adam-Multon Runge-Kutta de Ordem " + str(degree)
		plotxMultonRungeKutta, plotyMultonRungeKutta = adamMulton(plotyRungeKutta, plotxRungeKutta, h, step, expression, degree, methodTitle, True, True)
		lineMultonRungeKutta = pltQ.plot(plotxMultonRungeKutta, plotyMultonRungeKutta,  marker='o', color='xkcd:green', label="Moulton Runge-Kutta")

		methodTitle = "Metodo de Formula Inversa Euler de Ordem " + str(degree)
		plotxInverseFormulaEuler, plotyInverseFormulaEuler = inverseFormula(plotyEuler, plotxEuler, h, step, expression, degree, methodTitle, True, True)
		lineInverseFormulaEuler = pltQ.plot(plotxInverseFormulaEuler, plotyInverseFormulaEuler, marker='o', color='xkcd:lavender', label="Formula Inversa Euler")

		methodTitle = "Metodo de Formula Inversa Euler Inverso de Ordem " + str(degree)
		plotxInverseFormulaEulerInverse, plotyInverseFormulaEulerInverse = inverseFormula(plotyEulerInverse, plotxEulerInverse, h, step, expression, degree, methodTitle, True, True)
		lineInverseFormulaEulerInverse = pltQ.plot(plotxInverseFormulaEulerInverse, plotyInverseFormulaEulerInverse, marker='o', color='xkcd:purple', label="Formula Inversa Euler Inverso")

		methodTitle = "Metodo de Formula Inversa Euler Aprimorado de Ordem " + str(degree)
		plotxInverseFormulaEulerImproved, plotyInverseFormulaEulerImproved = inverseFormula(plotyEulerImproved, plotxEulerImproved, h, step, expression, degree, methodTitle, True, True)
		lineInverseFormulaEulerImproved = pltQ.plot(plotxInverseFormulaEulerImproved, plotyInverseFormulaEulerImproved, marker='o', color='xkcd:orange', label="Formula Inversa Euler Aprimorado")

		methodTitle = "Metodo de Formula Inversa Runge-Kutta de Ordem " + str(degree)
		plotxInverseFormulaRungeKutta, plotyInverseFormulaRungeKutta = inverseFormula(plotyRungeKutta, plotxRungeKutta, h, step, expression, degree, methodTitle, True, True)
		lineInverseFormulaRungeKutta = pltQ.plot(plotxInverseFormulaRungeKutta, plotyInverseFormulaRungeKutta, marker='o', color='xkcd:olive', label="Formula Inversa Runge-Kutta")
		
		pltQ.set(xlabel='T', ylabel='Y', title='Graficos dos Metodos')
		pltQ.grid()
		pltQ.tick_params(labelcolor='k', labelsize='small', width=1)
		pltQ.legend()
		plt.show()
		fig.savefig("GraficoMetodos.png")

	return

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