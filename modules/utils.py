import os, time

os.system('color 07')


def suppress_source(val):
	global ss 
	ss = val

####colours:
purple = '\033[95m'
blue = '\033[94m'
green = '\033[92m'
beige = '\033[93m'
red = '\033[91m'
end = '\033[0m'
white = '\033[1m'
under = '\033[4m'

def ts():
	return f'{purple}[{time.strftime("%H:%M:%S", time.localtime())}]{end}'

def message(source, text, colour='end'):
	if not ss:
		print(f'{ts()}{beige}({source}):{end} {globals()[colour]}{text}{end}')
	else:
		print(f'{ts()} {globals()[colour]}{text}{end}')

