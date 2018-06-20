import matplotlib.pyplot as plt
import csv
import numpy as np

planned = [[None], [None], [None], [None]]
measured = [[None], [None], [None], [None]]
t_measured = np.linspace(0, 8, 8.0/0.01)
t_planned = []

with open('times.csv','r') as csvfile:
	data = csv.reader(csvfile, delimiter=',')
	print data
	for row in data:
		for t in range(len(row)):
			t_planned.append(float(row[t]))

with open('states.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    s_idx = 0
    for row in data: # states
    	for s in range(len(row)):
        	planned[s_idx].append(float(row[s]))
        s_idx = s_idx+1

with open('states_measured.csv','r') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    s_idx = 0
    for row in data: # states
    	for s in range(len(row)):
        	measured[s_idx].append(float(row[s]))
        s_idx = s_idx+1

print t_measured
print len(t_measured)
print len(measured[0])

plt.plot(t_planned, planned[0][0:-1], label='planned')
plt.plot(t_measured, measured[0][0:len(t_measured)], label='measured')
plt.xlabel('t')
plt.ylabel('x')
plt.title('cart x-pos over time')
plt.legend()
plt.show()

plt.plot(t_planned, planned[1][0:-1], label='planned')
plt.plot(t_measured, measured[1][0:len(t_measured)], label='measured')
plt.xlabel('t')
plt.ylabel('theta')
plt.title('cart theta-pos over time')
plt.legend()
plt.show()