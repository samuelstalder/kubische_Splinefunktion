import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def spline(x,y,xx):
    n = len(x)-1
    #calculate a
    a = y
    #initialize parameters with zeros
    B = np.zeros(n)
    C = np.zeros(n+1) #C(n+1) necessary for algorithm
    D = np.zeros(n)
    h = np.zeros(n)
    #calculate h
    for i in range(n):
        h[i] = x[i+1] - x[i]
    
    #calculate c
    #initialize A & b with given shape of zeros
    A = np.zeros((n-1,n+1))
    b = np.zeros((n-1,1))
    for c in range(n-1):
        #insert into A to create diagonal-like matrix
        #inserts at row c with row-offset c
        A[c, c:3+c] = np.array([h[c], 2*(h[c]+h[c+1]), h[c+1]])
        #creates b
        b[c] = np.array([3*(y[c+2]-y[c+1])/h[c+1] - 3*(y[c+1]-y[c])/h[c]])
    A = A[:, 1:-1] #drop first and last column (only used for easier construction)
    print("A", A)
    print("b", b)
    #solve Ac = b
    c1 = np.linalg.solve(A,b)
    for i in range(len(c1)): #add solution to C array
        C[i+1] = c1[i]
    print("c", C)

    #calculate b & d
    for i in range(n):
        B[i] = (y[i+1] - y[i])/h[i] - h[i]*(C[i+1] + 2*C[i])/3
        D[i] = (C[i+1] - C[i])/(3*h[i])

    print("n", n)
    print("a", a)
    print("B", B)
    print("D", D)
    print("h", h)

    #plotting
    for i in range(n):
        xx = np.arange(x[i],x[i+1],0.01)
        yy = a[i] + B[i]*(xx-x[i]) + C[i]*(xx-x[i])**2 + D[i]*(xx-x[i])**3
        plt.plot(xx,yy, label="S"+str(i))
        
    #ev
    yy = np.zeros(xx.shape)
    for i in range(0,n):
        ind = np.where((xx >= x[i]) & (xx <= x[i+1]))
        t = xx[ind] - x[i]
        yy[ind] = a[i] + t*(B[i] + t*(C[i] + t*D[i]))
    
    plt.scatter(x,y)
    plt.legend()
    plt.xlabel("year")
    plt.ylabel("USA population (mio)")
    plt.show()
    
    return yy;

#Test
#x = np.array([0,1,2,3], dtype=np.float64)
#y = np.array([2,1,2,2], dtype=np.float64)

#2
x = np.array([4,6,8,10], dtype=np.float64)
y = np.array([6,3,9,0], dtype=np.float64)
spline(x,y,2)

#3a
x = np.array([1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010], dtype=np.float64)
y = np.array([75.995, 91.972, 105.711, 123.203, 131.669, 150.697, 179.323, 203.212, 226.505, 249.633, 281.422, 308.745], dtype=np.float64)
spline(x,y,2)

#3b
cs = CubicSpline(x,y)
xs = np.arange(1900,2010, 0.1)
plt.scatter(x,y)
plt.plot(xs, cs(xs))
plt.show()
#der plot sieht identisch aus

#3c
p = np.polyfit(x,y,11)
x_axis = np.arange(1900, 2010, 0.1)
plt.plot(x_axis, np.polyval(p, x_axis))
plt.scatter(x,y)
plt.show()


modified_x = x - 1900
modified_p = np.polyfit(modified_x,y,11)
modified_x_axis = x_axis - 1900
plt.scatter(x,y)
plt.plot(x_axis, np.polyval(modified_p, modified_x_axis))
plt.show()

#der plot sieht ähnlich aus aber "runder/weicher" oder weniger genau


