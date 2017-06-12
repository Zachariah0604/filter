import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

x = np.random.rand(50,30)
print(type(x))
#basic
f1 = plt.figure(1)
plt.subplot(211)
x1=x[:,1]
x2=x[:,0]
plt.scatter(x1,x2,s=0.01)
plt.show()