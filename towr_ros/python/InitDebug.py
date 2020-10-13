# author: Chao NI
# for debug use:
# plot the end effector position (z direction) and the decision flags

import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ee', default="all", help="which foot to plot: LF,RF,LH,RH")
parser.add_argument('--name', default="Initdebug", help="csv filename")

# reading data
opt = parser.parse_args()
filename = opt.name + ".csv"
data = pd.read_csv(filename, header=None)
data = data.values;

node_file_name = "NodesValue.csv"
data_node = pd.read_csv(node_file_name, header=None)
data_node = data_node.values

manual_node_file_name = "ManualNodesValue.csv"
manual_data_node = pd.read_csv(manual_node_file_name, header=None)
manual_data_node = manual_data_node.values

manual_ee_name = "manualEEPos.csv"
manual_ee_data = pd.read_csv(manual_ee_name, header=None)
manual_ee_data = manual_ee_data.values
# print(len(data)) 

# number of instances
n_instance = (len(data)) / 28

instance_axis = np.linspace(0,2.4,n_instance)
EE_z_0 = [] # end effector position in z direction
EE_z_1 = []
EE_z_2 = []
EE_z_3 = []

EE_x_0 = []
EE_x_1 = []
EE_x_2 = []
EE_x_3 = []

EE_d_0 = []  # decision
EE_d_1 = []
EE_d_2 = []
EE_d_3 = []

EE_c_0 = [] # contact flag
EE_c_1 = []
EE_c_2 = []
EE_c_3 = []

Node_z_0 = [] #11
Node_z_1 = [] #10
Node_z_2 = [] #8
Node_z_3 = [] #8

manual_Node_z_0 = [] #11
manual_Node_z_1 = [] #10
manual_Node_z_2 = [] #8
manual_Node_z_3 = [] #8

manual_ee_0 = []
manual_ee_1 = []
manual_ee_2 = []
manual_ee_3 = []


for i in range(n_instance):

	EE_x_0.append(data[0+28*i])
	EE_x_1.append(data[7+28*i])
	EE_x_2.append(data[14+28*i])
	EE_x_3.append(data[21+28*i])

	EE_z_0.append(data[2+28*i])
	EE_z_1.append(data[9+28*i])
	EE_z_2.append(data[16+28*i])
	EE_z_3.append(data[23+28*i])

	EE_d_0.append(data[3+28*i])
	EE_d_1.append(data[10+28*i])
	EE_d_2.append(data[17+28*i])
	EE_d_3.append(data[24+28*i])

	EE_c_0.append(data[6+28*i])
	EE_c_1.append(data[13+28*i])
	EE_c_2.append(data[20+28*i])
	EE_c_3.append(data[27+28*i])

	manual_ee_0.append(manual_ee_data[2+12*i])
	manual_ee_1.append(manual_ee_data[5+12*i])
	manual_ee_2.append(manual_ee_data[8+12*i])
	manual_ee_3.append(manual_ee_data[11+12*i])


len_0 = int(data_node[-4])
print(len_0)
len_1 = int(data_node[-3])
len_2 = int(data_node[-2])
len_3 = int(data_node[-1])

for i in range(len_0):
	Node_z_0.append(data_node[2+3*i])
	manual_Node_z_0.append(manual_data_node[2+3*i])
for i in range(len_1):
	Node_z_1.append(data_node[len_0*3+2+3*i])
	manual_Node_z_1.append(manual_data_node[len_0*3+2+3*i])
for i in range(len_2):
	Node_z_2.append(data_node[(len_0+len_1)*3+2+3*i])
	manual_Node_z_2.append(manual_data_node[(len_0+len_1)*3+2+3*i])
for i in range(len_3):
	Node_z_3.append(data_node[(len_0+len_1+len_2)*3+2+3*i])
	manual_Node_z_3.append(manual_data_node[(len_0+len_1+len_2)*3+2+3*i])


# Plot the data

plt.hlines(0.7,0,2.4,linestyles='dashed', label='Step Position')

if opt.ee =="all":

	plt.plot(instance_axis, EE_z_0, label='EE_z_0')
	plt.plot(instance_axis, EE_d_0, label='EE_d_0')
	plt.plot(instance_axis, EE_c_0, label='EE_c_0')


	plt.plot(instance_axis, EE_z_1, label='EE_z_1')
	plt.plot(instance_axis, EE_d_1, label='EE_d_1')
	plt.plot(instance_axis, EE_c_1, label='EE_c_1')



	plt.plot(instance_axis, EE_z_2, label='EE_z_2')
	plt.plot(instance_axis, EE_d_2, label='EE_d_2')
	plt.plot(instance_axis, EE_c_2, label='EE_c_2')


	plt.plot(instance_axis, EE_z_3, label='EE_z_3')
	plt.plot(instance_axis, EE_d_3, label='EE_d_3')
	plt.plot(instance_axis, EE_c_3, label='EE_c_3')

if (opt.ee=="LF"):

	plt.plot(instance_axis, EE_x_0, label='EE_x_0')
	plt.plot(instance_axis, EE_z_0, label='EE_z_0')
	plt.plot(instance_axis, EE_d_0, label='EE_d_0')
	plt.plot(instance_axis, EE_c_0, label='EE_c_0')
	instance_axis_node = np.linspace(0,2.4,len_0)
	# plt.plot(instance_axis_node, Node_z_0, label='Node_z_0' )
	# plt.plot(instance_axis_node, manual_Node_z_0, label='manual_Node_z_0' )
	# plt.plot(instance_axis, manual_ee_0, label='manual_ee_0')

if (opt.ee=="RF"):

	plt.plot(instance_axis, EE_x_1, label='EE_x_1')
	plt.plot(instance_axis, EE_z_1, label='EE_z_1')
	plt.plot(instance_axis, EE_d_1, label='EE_d_1')
	plt.plot(instance_axis, EE_c_1, label='EE_c_1')
	instance_axis_node = np.linspace(0,2.4,len_1)
	# plt.plot(instance_axis_node, Node_z_1, label='Node_z_1' )
	# plt.plot(instance_axis_node, manual_Node_z_1, label='manual_Node_z_1' )
	# plt.plot(instance_axis, manual_ee_1, label='manual_ee_1')

if (opt.ee=="LH"):

	plt.plot(instance_axis, EE_x_2, label='EE_x_2')
	plt.plot(instance_axis, EE_z_2, label='EE_z_2')
	plt.plot(instance_axis, EE_d_2, label='EE_d_2')
	plt.plot(instance_axis, EE_c_2, label='EE_c_2')
	instance_axis_node = np.linspace(0,2.4,len_2)
	# plt.plot(instance_axis_node, Node_z_2, label='Node_z_2' )
	# plt.plot(instance_axis_node, manual_Node_z_2, label='manual_Node_z_2' )
	# plt.plot(instance_axis, manual_ee_2, label='manual_ee_2')

if (opt.ee=="RH"):

	plt.plot(instance_axis, EE_x_3, label='EE_x_3')
	plt.plot(instance_axis, EE_z_3, label='EE_z_3')
	plt.plot(instance_axis, EE_d_3, label='EE_d_3')
	plt.plot(instance_axis, EE_c_3, label='EE_c_3')
	instance_axis_node = np.linspace(0,2.4,len_3)
	# plt.plot(instance_axis_node, Node_z_3, label='Node_z_3' )
	# plt.plot(instance_axis_node, manual_Node_z_3, label='manual_Node_z_3' )
	# plt.plot(instance_axis, manual_ee_3, label='manual_ee_3')


plt.legend()
plt.xlabel("time")
plt.ylabel("related values")
plt.title(opt.name)
plt.show()