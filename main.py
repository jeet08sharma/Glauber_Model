import numpy as np
import matplotlib.pyplot as plt
import random

#parametrs
imp_parameter = 6
sigma = 42
radius_of_nucleus_A = 6.38
radius_of_nucleus_B = 6.38
mass_number_of_nucleus_A = 197
mass_number_of_nucleus_B = 197
skin_depth_of_nucleus_A = 0.535
skin_depth_of_nucleus_B = 0.535


#random number generator
def generate_random_numbers(lower_limit, upper_limit, num_points):
    random_numbers = []

    for _ in range(num_points):
        random_number = random.uniform(lower_limit, upper_limit)
        random_numbers.append(random_number)

    return random_numbers

#position of nucleons
def positions_of_nucleons(radius, number_of_nucleons):
    x_values = []
    y_values = []
    z_values = []
    theta = np.random.uniform(0, 360, number_of_nucleons)
    phi = np.random.uniform(0, 180, number_of_nucleons)
    x = radius * np.sin(phi) * np.cos(theta)
    x_values.append(x)
    y = radius * np.sin(phi) * np.sin(theta)
    y_values.append(y)
    z = radius * np.cos(phi)
    z_values.append(z)
    positions = np.column_stack((x, y, z))
    return positions,x_values,y_values,z_values

#radius of each nucleons assuming that, nucleons are moving in z direction
def radii(nucleus):
    r = []
    for i in range(len(nucleus)):
        sum = 0
        for j in range(len(nucleus[i])):
            sum += nucleus[i][j] ** 2
            if j == 1:
                break
        r.append(sum)
    return r

#density of the nucleus
def density(r, R, a):
    rho = []
    for i in range(len(r)):
        den = 1 / (1 + (np.exp((r[i] - R) / a)))
        rho.append(den)
    return rho

#distance between the nucloens
def dist(x1values,y1values,x2values,y2values):
    dist_=[]
    for i in range(len(x1values)):
        for j in range(len(x1values)):
            xdiff=x1values[i]-x2values[j]
            ydiff=y1values[i]-y2values[j]
            r=np.sqrt(pow(xdiff,2)+pow(ydiff,2))
            dist_.append(r)
    return dist_

nucleus_A,xA,yA,zA = positions_of_nucleons(radius_of_nucleus_A,mass_number_of_nucleus_A)
nucleus_B,xB,yB,zB = positions_of_nucleons(radius_of_nucleus_B,mass_number_of_nucleus_B)
r_A = radii(nucleus_A)
r_B = radii(nucleus_B)
x_A = np.array(xA).flatten()
x_B = np.array(xB).flatten()
Y_A = np.array(yA).flatten()
y_B = np.array(yB).flatten()
z_A = np.array(zA).flatten()
z_B = np.array(zB).flatten()


# change of frame
y_A = [x+imp_parameter for x in Y_A]

# distance between the nucleons
def dist(x1values,y1values,x2values,y2values):
    dist_=[]
    for i in range(len(x1values)):
        for j in range(len(x1values)):
            xdiff=x1values[i]-x2values[j]
            ydiff=y1values[i]-y2values[j]
            r=np.sqrt(pow(xdiff,2)+pow(ydiff,2))
            dist_.append(r)
    return dist_
distance = dist(x_A,y_A,x_B,y_B)

#number of collision
def coll(distance,sigma):

    collision = 0
    for i in range(len(distance)):
        if distance[i] < np.sqrt(sigma/(np.pi)):
            collision+=1
    t = 1 - np.exp(-(collision/(208*208)))
    return collision,t
col,pro = coll(distance,42)


# Scatter plot of nuclei A and B

point = 100
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax.scatter(nucleus_A[:, 0], nucleus_A[:, 1], nucleus_A[:, 2],s = point, color='blue',label ='Au')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Nucleus A')
plt.legend()

ax = fig.add_subplot(122, projection='3d')
ax.scatter(nucleus_B[:, 0], nucleus_B[:, 1], nucleus_B[:, 2],s=point, color='red',label = 'Au')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Nucleus B')
plt.legend()
plt.show()

# Scatter plot of density of nuclei A and B

plt.scatter(r_A, density(r_A, radius_of_nucleus_A, skin_depth_of_nucleus_A), color='blue',label = 'Au')
plt.scatter(r_B, density(r_B, radius_of_nucleus_B, skin_depth_of_nucleus_B), color='red',label = 'Au')


plt.xlabel('radius in fm')
plt.ylabel('density rho(r)/rho_0')
plt.axvline(x=radius_of_nucleus_A,color = 'blue',label='R = 6.38 in fm')
plt.axvline(x=radius_of_nucleus_B,color='red',label ='R = 6.38 in fm')

plt.legend()
plt.show()

#colliding nucleus
point_size = 300
color = np.random.rand(mass_number_of_nucleus_A)
colors = np.random.rand(mass_number_of_nucleus_B)
plt.scatter(y_A,x_A, s = point_size, c=color,cmap = 'Blues',label ='Au')
plt.scatter(y_B,x_B, s = point_size,c =colors,cmap  = 'Reds',label = 'Au')
plt.title("Cu nuclues at impact parameter of 6fm")

plt.xlabel('y')
plt.ylabel('x')
plt.legend()

plt.show()

b = generate_random_numbers(0,20,70)
proba,cool,sigp=[],[],[]
for i in range(len(b)):
    yy_A = [x + b[i] for x in Y_A]
    D=dist(x_A,yy_A,x_B,y_B)
    c,p=coll(D,sigma)
    cool.append(c)
    proba.append(p)
    sigp.append(p*b[i])

plt.scatter(b,sigp)
plt.xlabel('impact parameter')
plt.ylabel('cross section')
plt.title('Au+Au')
plt.show()

plt.scatter(b,cool)
plt.show()


