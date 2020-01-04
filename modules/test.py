def hash(n):
	return (n%3>0, (n+1)%3>0, (n+2)%2>0)

print(hash(-1))