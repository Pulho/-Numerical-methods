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

def adamBashforth(y0, t0, h, step, exp, degree, method, flag):
	xAdam = []
	yAdam = []

	printHeader(method, y0[0], t0[0], h)

	expression = sympify(exp)
	T = t0[0]

	for i in range(degree):
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
		outputAppend(str(i+1) + " " + str(yni) + "\n")
		yAdam.append(yni)
		T += h
		xAdam.append(T)
		
	outputAppend("\n")
	return xAdam, yAdam

def adamMulton(y0, t0, h, step, exp, degree, method, flag):
	xAdam = []
	yAdam = []

	printHeader(method, y0[0], t0[0], h)

	expression = sympify(exp)
	T = t0[0]
	
	for i in range(degree):
		outputAppend(str(i) + " " + str(y0[i]) + "\n")
		xAdam.append(T)
		yAdam.append(y0[i])
		if flag:
			T = t0[i]
		else:
			T = T + h
	T = t0[0] + (degree - 1) * h

	for i in range(degree - 1, step): 
		if degree == 1:
			Tn1 = T + h 
			Yn1 = yAdam[i] + h*expression.subs([(t, T), (y, yAdam[i])]) 

			F1 = (1/2)*expression.subs([(t, Tn1), (y, Yn1)])
			F0 = (1/2)*expression.subs([(t, T), (y, yAdam[i])])

			yni = yAdam[i] + h*(F1 + F0)
			yni += (-1/12)
		elif degree == 2:
			Tn1 = T + h 
			Yn1 = yAdam[i] + h*expression.subs([(t, T), (y, yAdam[i])]) 

			F2 = (5/12)*expression.subs([(t, Tn1), (y, Yn1)])
			F1 = (2/3)*expression.subs([(t, T), (y, yAdam[i])])
			F0 = (-1/12)*expression.subs([(t, T-h), (y, yAdam[i-1])])

			yni = yAdam[i] + h*(F2 + F1 + F0)
			yni += (1/24)
		elif degree == 3:
			Tn1 = T + h 
			Yn1 = yAdam[i] + h*expression.subs([(t, T), (y, yAdam[i])]) 

			F3 = (3/8)*expression.subs([(t, Tn1), (y, Yn1)])
			F2 = (19/24)*expression.subs([(t, T), (y, yAdam[i])])
			F1 = (-5/24)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F0 = (1/24)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])

			yni = yAdam[i] + h*(F3 + F2 + F1 + F0)
			yni += (-19/720)
		elif degree == 4:
			Tn1 = T + h 
			Yn1 = yAdam[i] + h*expression.subs([(t, T), (y, yAdam[i])]) 

			F4 = (251/720)*expression.subs([(t, Tn1), (y, Yn1)])
			F3 = (323/360)*expression.subs([(t, T), (y, yAdam[i])])
			F2 = (-11/30)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F1 = (53/360)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F0 = (3/160)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])

			yni = yAdam[i] + h*(F4 + F3 + F2 + F1 + F0)
			yni += (-19/720)
		elif degree == 5:
			Tn1 = T + h 
			Yn1 = yAdam[i] + h*expression.subs([(t, T), (y, yAdam[i])]) 

			F5 = (95/288)*expression.subs([(t, Tn1), (y, Yn1)])
			F4 = (1427/1440)*expression.subs([(t, T), (y, yAdam[i])])
			F3 = (-133/240)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F2 = (241/720)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F1 = (-173/1440)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F0 = (3/160)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			
			yni = yAdam[i] + h*(F5 + F4 + F3 + F2 + F1 + F0)
			yni += (-863/60480)
		elif degree == 6:
			Tn1 = T + h 
			Yn1 = yAdam[i-1] + h*expression.subs([(t, T), (y, yAdam[i-1])]) 

			F6 = (19087/60480)*expression.subs([(t, Tn1), (y, Yn1)])
			F5 = (2713/2520)*expression.subs([(t, T), (y, yAdam[i-1])])
			F4 = (-15487/20160)*expression.subs([(t, T-h), (y, yAdam[i-2])])
			F3 = (586/945)*expression.subs([(t, T-(2*h)), (y, yAdam[i-3])])
			F2 = (-6737/20160)*expression.subs([(t, T-(3*h)), (y, yAdam[i-4])])
			F1 = (263/2520)*expression.subs([(t, T-(4*h)), (y, yAdam[i-5])])
			F0 = (-863/60480)*expression.subs([(t, T-(5*h)), (y, yAdam[i-6])])
			
			yni = yAdam[i-1] + h*(F6 + F5 + F4 + F3 + F2 + F1 + F0)
			yni += (275/24192)
		elif degree == 7:
			Tn1 = T + h 
			Yn1 = yAdam[i] + h*expression.subs([(t, T), (y, yAdam[i])]) 

			F7 = (5257/17280)*expression.subs([(t, Tn1), (y, Yn1)])
			F6 = (139849/120960)*expression.subs([(t, T), (y, yAdam[i])])
			F5 = (-4511/4480)*expression.subs([(t, T-h), (y, yAdam[i-1])])
			F4 = (123133/120960)*expression.subs([(t, T-(2*h)), (y, yAdam[i-2])])
			F3 = (-88547/120960)*expression.subs([(t, T-(3*h)), (y, yAdam[i-3])])
			F2 = (1537/4480)*expression.subs([(t, T-(4*h)), (y, yAdam[i-4])])
			F1 = (-11351/120960)*expression.subs([(t, T-(5*h)), (y, yAdam[i-5])])
			F0 = (275/24192)*expression.subs([(t, T-(6*h)), (y, yAdam[i-6])])

			yni = yAdam[i] + h*(F7 + F6 + F5 + F4 + F3 + F2 + F1 + F0)
			yni += (-33953/3628800)
		outputAppend(str(i+1) + " " + str(yni) + "\n")
		yAdam.append(yni)
		T += h
		xAdam.append(T)
	outputAppend("\n")
	return xAdam, yAdam # WORK IN PROGRESS

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

	elif strct[0] == "adam_bashforth":
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

		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0, t0, h, step, expression, degree-1, methodTitle, False)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_euler":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Euler de Ordem " + str(degree)
		fileImageTitle = "bashforthEuler" + str(degree) + "th"
		x0Euler, y0Euler = euler(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0Euler, x0Euler, h, step, expression, degree, methodTitle, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_euler_aprimorado":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Euler aprimorado de Ordem " + str(degree)
		fileImageTitle = "bashforthEulerImproved" + str(degree) + "th"
		x0EulerImproved, y0EulerImproved = eulerImproved(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0EulerImproved, x0EulerImproved, h, step, expression, degree, methodTitle, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_euler_inverso":  
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Euler Inverso de Ordem " + str(degree)
		fileImageTitle = "bashforthEulerInverse" + str(degree) + "th"
		x0EulerInverse, y0EulerInverse = eulerInverse(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0EulerInverse, x0EulerInverse, h, step, expression, degree, methodTitle, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_bashforth_by_runge_kutta": # Faltando: + 1
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Bashforth por Runge-Kutta de Ordem " + str(degree)
		fileImageTitle = "bashforthRungeKutta" + str(degree) + "th"
		x0RungeKutta, y0RungeKutta = rungeKutta(y0, t0, h, step, expression, False)
		plotxAdamBashforth, plotyAdamBashforth = adamBashforth(y0RungeKutta, x0RungeKutta, h, step, expression, degree, methodTitle, True)
		plotGraphic(plotxAdamBashforth,plotyAdamBashforth, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon":
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

		plotxMulton, plotyMulton = adamMulton(y0, t0, h, step, expression, degree-1, methodTitle, False)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon_by_euler":
		y0 = float(strct[1])
		t0 = float(strct[2])
		h =  float(strct[3])
		step = int(strct[4])
		expression = strct[5]
		degree = int(strct[6])

		methodTitle = "Metodo de Adam-Multon por Euler de Ordem " + str(degree)
		fileImageTitle = "bashforthMulton" + str(degree) + "th"
		x0Euler, y0Euler = euler(y0, t0, h, step, expression, False)
		plotxMulton, plotyMulton = adamMulton(y0Euler, x0Euler, h, step, expression, degree, methodTitle, True)
		plotGraphic(plotxMulton,plotyMulton, methodTitle, fileImageTitle)

	elif strct[0] == "adam_multon_by_euler_inverso":
		adamMultonEulerInverse(strct)
	elif strct[0] == "adam_multon_by_euler_aprimorado":
		adamMultonEulerImproved(strct)
	elif strct[0] == "adam_multon_by_runge_kutta":
		adamMultonRungeKutta(strct)
	elif strct[0] == "formula_inversa":
		inverseFormula(strct)
	elif strct[0] == "formula_inversa_by_euler":
		inverseFormulaEuler(strct)
	elif strct[0] == "formula_inversa_by_euler_inverso":
		inverseFormulaEulerInverse(strct)
	elif strct[0] == "formula_inversa_by_euler_aprimorado":
		inverseFormulaEulerImproved(strct)
	elif strct[0] == "formula_inversa_by_runge_kutta":
		inverseFormulaRungeKutta(strct)
	elif strct[0] == "all": 
		euler(strct)
		eulerInverse(strct)
		eulerImproved(strct)
		rungeKutta(strct)
		adamBashforth(strct)
		adamBashforthEuler(strct)
		adamBashforthEulerImproved(strct)
		adamBashforthEulerInverse(strct)
		adamBashforthRungeKutta(strct)
		adamMultonEuler(strct)
		adamMultonEulerInverse(strct)
		adamMultonEulerImproved(strct)
		adamMultonRungeKutta(strct)
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