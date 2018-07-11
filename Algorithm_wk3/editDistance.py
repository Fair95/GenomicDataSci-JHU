def editDistRecursive(x, y):
	# This implementation is very slow
	if len(x) == 0:
		return len(y)
	elif len(y) == 0:
		return len(x)
	else:
		distHor = editDistRecursive(x[:-1],y) + 1
		distVer = editDistRecursive(x, y[:-1]) + 1
		if x[-1] == y[-1]:
			distDiag = editDistRecursive(x[:-1],y[:-1])
		else:
			distDiag = editDistRecursive(x[:-1],y[:-1]) + 1

		return min(distHor,distVer,distDiag)

def editDistance(x, y):
	D = []
	# initialize the matrix D with size (len(x)+1) * (len(y)+1)
	for i in range(len(x)+1):
		D.append([0]* (len(y)+1))

	# the distance between empty string and a string L is always len(L)
	for i in range(len(x)+1):
		D[i][0] = i
	for i in range(len(y)+1):
		D[0][i] = i

	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			distHor = D[i][j-1] + 1
			distVer = D[i-1][j] + 1
			if x[i-1] == y[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + 1

			D[i][j] = min(distHor, distVer, distDiag)
	# # Edit distance is the value in the bottom right corner of the matrix
	return D[-1][-1]

def editDistance_modified(x, y):
	D = []
	# initialize the matrix D with size (len(x)+1) * (len(y)+1)
	for i in range(len(x)+1):
		D.append([0]* (len(y)+1))

	# the distance between empty string and a string L is always len(L)
	for i in range(len(x)+1):
		D[i][0] = i
	for i in range(len(y)+1):
		D[0][i] = 0

	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			distHor = D[i][j-1] + 1
			distVer = D[i-1][j] + 1
			if x[i-1] == y[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + 1

			D[i][j] = min(distHor, distVer, distDiag)
	# # Edit distance is the value in the bottom right corner of the matrix
	return min(D[i])
