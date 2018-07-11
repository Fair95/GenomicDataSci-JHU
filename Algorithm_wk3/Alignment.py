alphabet = ['A','C','G','T']
## A<->G and C<->T have score 2, skipping has score 8, others have score 4
score = [[0, 4, 2, 4, 8], #A
		 [4, 0, 4, 2, 8], #C
		 [2, 4, 0, 4, 8], #G
		 [4, 2, 4, 0, 8]] #T
		 #A #C #G #T #skip
def globalAlignment(x, y):
	D = []
	# initialize the matrix D with size (len(x)+1) * (len(y)+1)
	for i in range(len(x)+1):
		D.append([0]* (len(y)+1))
	# top left cornor will always be 0, alphabet.index(x[i-1]) gives the char we at
	# initialize the first row and first colomn all be 8
	for i in range(1, len(x)+1):
		D[i][0] = D[i-1][0] + score[alphabet.index(x[i-1])][-1]
	for i in range(1, len(y)+1):
		D[0][i] = D[0][i-1] + score[-1][alphabet.index(y[i-1])]

	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			# skip a distance in either x or y will have penalty 8
			distHor = D[i][j-1] + score[-1][alphabet.index(y[j-1])]
			distVer = D[i-1][j] + score[alphabet.index(x[i-1])][-1]
			if x[i-1] == y[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

			D[i][j] = min(distHor, distVer, distDiag)
	return D[-1][-1]