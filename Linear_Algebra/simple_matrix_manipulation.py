def DotProduct(x,y):
	
# Simple Case
	if (type(x) is int) and (type(y) is int):
		return x*y

# Data Checks		
	if len(x)!=len(y):
		return None, "Error: vectors of differing lengths"
	for i in x:
		if (((type(i) is int) == False) and ((type(i) is float) == False)):
			return None, "Error: one or both vectors contain non-numeric values"
	for i in y:
		if (((type(i) is int) == False) and ((type(i) is float) == False)):
			return None, "Error: one or both vectors contain non-numeric values"

# Evaluate Dot Product			
	dot_product = 0	
	i = 0	
	while(i<len(x)):
		dot_product = dot_product + x[i]*y[i]
		i = i + 1
	return dot_product

def zeroMatrix(p,q):
	A = [0.]*p
	for i in range(p):
		A[i] = [0.]*q		
	return A

def transpose(A):
	p = len(A)
	q = len(A[0])	
	B = zeroMatrix(q,p)
	
	for i in range(q):
		for j in range(p):
			B[i][j] = A[j][i]			
	return B
	
def matrixMult(A,B):

	# get dimensions
	p = len(A)
	m_A = len(A[0])
	m_B = len(B)
	q = len(B[0])
	
	# check that multiplication is defined
	if (m_A!=m_B):
		return None, "Error: matrix multiplication undefined for given matrices"	

	# instantiate matrix
	C = zeroMatrix(p,q)
	
	# transpose B
	Bt = transpose(B)
	
	for i in range(p):
		for j in range(q):
			C[i][j] = DotProduct(A[i],Bt[j])
			
	return C
	
A_1 = [[1,3,2],[4,5,6],[8,7,9]]
B_1 = [[0.2,-0.866667,0.533333],[0.8,-0.466667,1.33333],[-0.8,1.13333,-0.466667]]

A_2 = [[1,2],[3,0],[5,7]]
B_2 = [[7,2],[6,8]]

print matrixMult(A_1,B_1)
print matrixMult(A_2,B_2)
print matrixMult(B_2,A_2)