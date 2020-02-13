import os, time
import inspect

os.system('color 07')

ps = cm = pt = False
vb = 0
def ff_print_source(val):
	global ps 
	ps = val

def ff_use_colours(val):
	global cm
	cm = val

def ff_print_time(val):
	global pt
	pt = val

def ff_verbosity(val):
	global vb
	vb = val


####colours:

def c(val):
	if cm:
		if val == "purple": return '\033[95m'
		if val == "blue": return '\033[94m'
		if val == "green": return '\033[92m'
		if val == "beige": return '\033[93m'
		if val == "red": return '\033[91m'
		if val == "end": return '\033[0m'
		if val == "white": return '\033[1m'
		if val == "under": return '\033[4m'
	else: return ''

		
def ts():
	return f'{c("purple")}[{time.strftime("%H:%M:%S", time.localtime())}]{c("end")}'

def message(text='', verbosity=0, colour='end'):
	try:
		#check if multiple texts or verbosities are given
		_ = [e for e in verbosity]

	except:
		#if there is only one
		if vb >= verbosity:
			source = inspect.stack()
			source_file = source[1][1].split('\\')[-1][:-3]
			source_func = source[1][3]
			
			if source_func == '<module>':
				source = source_file  + '.py'
			elif source_func == '<listcomp>':
				source = source_file
			else:
				source = source_file  + '.' + source_func

			if pt:
				if ps:
					print(f'{ts()}{c("beige")}({source}):{c("end")} {c(colour)}{text}{c("end")}')
				else:
					print(f'{ts()} {c(colour)}{text}{c("end")}')
			else:
				if ps:
					print(f'{c("beige")}({source}):{c("end")} {c(colour)}{text}{c("end")}')
				else:
					print(f'{c(colour)}{text}{c("end")}')

	else:
		#if there are multiple verbosities and texts
		assert len(text) == len(verbosity)

		#get correct verbosity:
		if vb >= verbosity[-1]: m = text[-1]
		else:
			try:
				v = verbosity.index(vb)
				m = text[v]
			except:
				d = [v for v in verbosity if vb -v > 0]
				if len(d) == 0:
					return
				m = text[verbosity.index(d[-1])]

		source = inspect.stack()
		source_file = source[1][1].split('\\')[-1][:-3]
		source_func = source[1][3]
		
		if source_func == '<module>':
			source = source_file  + '.py'
		elif source_func == '<listcomp>':
			source = source_file
		else:
			source = source_file  + '.' + source_func

		if pt:
			if ps:
				print(f'{ts()}{c("beige")}({source}):{c("end")} {c(colour)}{m}{c("end")}')
			else:
				print(f'{ts()} {c(colour)}{m}{c("end")}')
		else:
			if ps:
				print(f'{c("beige")}({source}):{c("end")} {c(colour)}{m}{c("end")}')
			else:
				print(f'{c(colour)}{m}{c("end")}')