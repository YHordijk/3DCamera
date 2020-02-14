#partial derivatives
import numpy as np

def function(x,y):
	return 5*x**2 - y**2
	#deriv (x) = 10x
	#deriv (y) = -2y
	#second deriv (xx) = 10
	#second deriv (xy) = 0
	#second deriv (yx) = 0
	#second deriv (yy) = -2

def first_deriv(func, inputs, d=1e-8):
	inputs = np.asarray(inputs).astype(float)
	f0 = func(*inputs)
	J = np.zeros(len(inputs))
	for i in range(len(inputs)):
		inputs_copy = inputs.copy()
		inputs_copy[i] += d
		J[i] = func(*inputs_copy)

	return (J-f0)/d

def hessian(func, inputs, d=1e-3):
	inputs = np.asarray(inputs).astype(float)
	f0 = func(*inputs)
	H = np.zeros((len(inputs), len(inputs)))
	# print(H)
	for i in range(len(inputs)):
		inputs_copy = inputs.copy()
		inputs_copy[i] += d
		f1 = func(*inputs_copy)
		for j in range(len(inputs)):
			inputs_copy[j] += d
			print(inputs_copy)
			f2 = func(*inputs_copy)

			H[i,j] = (f2-f1)

			inputs_copy[j] -= d

	return H

# print(first_deriv(function, [5, 1]))
print(hessian(function, [5,1]))