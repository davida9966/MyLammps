###########################################################################

# read trj file to determine how many Li enters into and how many drags Si out
# the index of Si is from 1 to n

# Author: Ya0 Mingze   Date: June 8th 2019

###########################################################################


import re
import numpy as np
import os.path

amount_of_Si = 864

input_trj_name = 'Silicene_Lattice_NPT.trj'
output_fname = 'how_many_Li_new_var_new_status.txt'
fin = open(input_trj_name,'r')
fout = open(output_fname,'w')




# read the fore-header
a = fin.readline()
a = re.sub("\D","",a)
a = int(a)
for i in range(0,a):
	t = fin.read(1)
text = fin.readline()
temp, char1, line1 = map(str,text.split())
line1 = int(re.sub("\D","",line1)) # read how many atoms (including Si and Li) in total
for i in range(1,line1+1):
	temp = fin.readline()




# read the 0 timestep to determine the "layer" of Si atoms according to z coordinate
# frame header
a = fin.readline()
a = re.sub("\D","",a)
a = int(a)
step = fin.readline()
char_length = len(step)
step = re.sub("\D","",step)
timestep = int(step)
print("Timestep %d is processed! Reading atom coordinate" %(timestep))
fin.read(a-char_length)
# coordinate
temp, char1, natom = map(str,fin.readline().split())
natom = int(re.sub("\D","",natom))
layer_of_Si = [[]*1 for i in range(0,9)] # initializing an empty matrix
Li_z_coordinate_0 = [0]*(natom+1) # initializing an empty matrix to store z coordinate of Li atoms
for i in range(1,natom+1):
	temp = [float(s) for s in fin.readline().split()]
	if round(temp[0])<=amount_of_Si:
		# this is a Si atom, it will be classified according to its z coordinate
		if temp[3]>22:
			layer_of_Si[0].append(round(temp[0]))
		elif temp[3]<22 and temp[3]>18:
			layer_of_Si[1].append(round(temp[0]))
		elif temp[3]<18 and temp[3]>15:
			layer_of_Si[2].append(round(temp[0]))
		elif temp[3]<15 and temp[3]>12:
			layer_of_Si[3].append(round(temp[0]))
		elif temp[3]<12 and temp[3]>9:
			layer_of_Si[4].append(round(temp[0]))
		elif temp[3]<9  and temp[3]>7:
			layer_of_Si[5].append(round(temp[0]))
		elif temp[3]<7  and temp[3]>4:
			layer_of_Si[6].append(round(temp[0]))
		elif temp[3]<4  and temp[3]>2:
			layer_of_Si[7].append(round(temp[0]))
		elif temp[3]<2  and temp[3]>0:
			layer_of_Si[8].append(round(temp[0]))
	elif round(temp[0])> amount_of_Si:
		# this is a Li atom, its z coordinate will be stroed in the matrix Li_z_coordinate
		Li_z_coordinate_0[round(temp[0])] = temp[3]
# bonds
temp, char1, nbond = map(str,fin.readline().split())
nbond = int(re.sub("\D","",nbond))
Si_Si_bond_intralayer = []  # initializing a empty matrix
bond_info_Si_Li_0 = [] # initializing, write down Si-Li bonds, the first element is Si and the second is Li
bond_info_Si_Si_0 = [] # initializing, write down Si-Si bonds
for i in range(1,nbond+1):
	temp = [float(s) for s in fin.readline().split()]
	if round(temp[0])<=amount_of_Si and round(temp[1])<=amount_of_Si:
		# this is a Si-Si bond, we need to specify whether it is intralayer
		bond_info_Si_Si_0.append( [ round(temp[0]),round(temp[1]) ] )
		if abs(round(temp[0])-round(temp[1]))<96:
			# this is a intralayer Si-Si bond
			Si_Si_bond_intralayer.append([round(temp[0]),round(temp[1])])
	elif round(temp[0])<=amount_of_Si and round(temp[1])>amount_of_Si:
		# this is a Si-Li bond
		bond_info_Si_Li_0.append( [ round(temp[0]),round(temp[1]) ] )
	elif round(temp[0])>amount_of_Si and round(temp[1])<=amount_of_Si:
		# this is a Li-Si bond
		bond_info_Si_Li_0.append( [ round(temp[1]),round(temp[0]) ] )






# initializing comparing matrices
Si_Si_bond_intralayer_0 = [[0,0] for i in range(0,len(Si_Si_bond_intralayer))] # the first element would always be larger than the second
Si_Si_bond_intralayer_1 = []  # these two matrices will be compared with each other to decide which Si-Si bonds break
for i in range(0,len(Si_Si_bond_intralayer)):
	for j in range(0,2):
		Si_Si_bond_intralayer_0[i][j] = Si_Si_bond_intralayer[i][j]
	if Si_Si_bond_intralayer_0[i][0]<Si_Si_bond_intralayer_0[i][1]:
		temp1 = Si_Si_bond_intralayer_0[i][0]
		temp2 = Si_Si_bond_intralayer_0[i][1]
		Si_Si_bond_intralayer_0[i][0] = temp2
		Si_Si_bond_intralayer_0[i][1] = temp1

bond_info_Si_Li_1 = [] # initializing, write down Si-Li bonds, the first element is Si and the second is Li
bond_info_Si_Si_1 = [] # initializing, write down Si-Si bonds

Li_z_coordinate_1 = [0]*(natom+1) # initializing an empty matrix to store z coordinate of Li atoms

Li_enter = 0 # initializing how many Li enters into
Li_drag = 0 # initializing how many Li drags Si out

Si_layer_larger_than_me_is_contact = 0

Li_status = [[0,0,0] for i in range(0,natom+1)]  # index, is_enter, is_drag; if is_enter == 1, then is_drag must be 0
for i in range(0,len(Li_status)):
	Li_status[i][0] = i



# read other timesteps to determine how many Li
# go into text
while timestep<=400000:
	# frame header
	a = fin.readline()
	if a == '':
		print('Exec end!')
		break
	a = re.sub("\D","",a)
	a = int(a)
	step = fin.readline()
	char_length = len(step)
	step = re.sub("\D","",step)
	timestep = int(step)
	print("Timestep %d is processed! Now layer larger than %3d is intact" %(timestep,Si_layer_larger_than_me_is_contact))
	fin.read(a-char_length)
	# coordinate
	temp, char1, natom = map(str,fin.readline().split())
	natom = int(re.sub("\D","",natom))
	Si_z_coordinate = [0]*len(layer_of_Si)
	Si_z_min_max = [[0,100,0] for i in range(0,len(layer_of_Si))] # max,min, max-min
	Si_layer_for_var = [[]*1 for i in range(0,len(layer_of_Si))]
	for i in range(1,natom+1):
		temp = [float(s) for s in fin.readline().split()]
		if round(temp[0])<=amount_of_Si:
			# this is a Si atom
			for i in range(0,len(layer_of_Si)):
				if round(temp[0]) in layer_of_Si[i]:
					Si_z_coordinate[i] = Si_z_coordinate[i]+temp[3]
					Si_layer_for_var[i].append(temp[3])
					if temp[3]>Si_z_min_max[i][0]:
						Si_z_min_max[i][0] = temp[3]
					if temp[3]<Si_z_min_max[i][1]:
						Si_z_min_max[i][1] = temp[3]
		else:
			# this is a Li atom
			Li_z_coordinate_1[round(temp[0])] = temp[3]
	for i in range(0,len(Si_z_coordinate)):
		Si_z_coordinate[i] = Si_z_coordinate[i] / len(layer_of_Si[i])
	for i in range(0,len(Si_z_min_max)):
		Si_z_min_max[i][2] = Si_z_min_max[i][0] - Si_z_min_max[i][1]

	Si_layer_larger_than_me_is_contact = 0 # 有些层随着时间变大会坏掉，这种情况下就不要考虑Li破坏这些已经被破坏的层的Si-Si键
	for i in range(0,len(Si_z_min_max)):
		layer_cal_var_np = np.array(Si_layer_for_var[i])
		var_now = np.var(layer_cal_var_np)
		if var_now<=9:
			Si_layer_larger_than_me_is_contact = i
			break

	# bonds
	temp, char1, nbond = map(str,fin.readline().split())
	nbond = int(re.sub("\D","",nbond))
	for i in range(1,nbond+1):
		temp = [float(s) for s in fin.readline().split()]
		if round(temp[0])<=amount_of_Si and round(temp[1])<=amount_of_Si:
			# this is a Si-Si bond, but we don't know whether it is intralayer
			bond_info_Si_Si_1.append( [ round(temp[0]),round(temp[1]) ] )
			if abs(round(temp[0])-round(temp[1]))<96:
				# this is an intralayer Si-Si bond
				if round(temp[0])<round(temp[1]):
					Si_Si_bond_intralayer_1.append( [ round(temp[1]),round(temp[0]) ] )
				else:
					Si_Si_bond_intralayer_1.append( [ round(temp[0]),round(temp[1]) ] )
		elif round(temp[0])<=amount_of_Si and round(temp[1])>amount_of_Si:
			# this is a Si-Li bond
			bond_info_Si_Li_1.append( [ round(temp[0]),round(temp[1]) ] )
		elif round(temp[0])>amount_of_Si and round(temp[1])<=amount_of_Si:
			# this is a Li-Si bond
			bond_info_Si_Li_1.append( [ round(temp[1]),round(temp[0]) ] )

	# processing bond info
	# decide whether Li atoms go deeper
	if timestep%1000 == 0:
		class_of_Li_z_coordinate_0 = [-1000]*(natom+1)
		class_of_Li_z_coordinate_1 = [-1000]*(natom+1)
		for i in range(amount_of_Si+1,natom):
			# classify each Li atom
			for j in range(0,len(Si_z_coordinate)):
				if j == 0:
					if Li_z_coordinate_0[i] > (Si_z_coordinate[0]+Si_z_coordinate[1])/2 + 1 and Li_z_coordinate_0[i] <Si_z_coordinate[j] + 2.5:
						class_of_Li_z_coordinate_0[i] = j
					if Li_z_coordinate_1[i] > (Si_z_coordinate[0]+Si_z_coordinate[1])/2 + 1 and Li_z_coordinate_1[i] <Si_z_coordinate[j] + 2.5:
						class_of_Li_z_coordinate_1[i] = j
				else:
					if Li_z_coordinate_0[i]<(Si_z_coordinate[j-1]+Si_z_coordinate[j])/2 + 1 and Li_z_coordinate_0[i]>(Si_z_coordinate[j-1]+Si_z_coordinate[j])/2 - 1:
						class_of_Li_z_coordinate_0[i] = j
					if Li_z_coordinate_1[i]<(Si_z_coordinate[j-1]+Si_z_coordinate[j])/2 + 1 and Li_z_coordinate_1[i]>(Si_z_coordinate[j-1]+Si_z_coordinate[j])/2 - 1:
						class_of_Li_z_coordinate_1[i] = j
		for i in range(amount_of_Si+1,natom):
			if class_of_Li_z_coordinate_0[i] != -1000 and class_of_Li_z_coordinate_1[i] != -1000 and (class_of_Li_z_coordinate_1[i]-class_of_Li_z_coordinate_0[i])>=1:
				Li_enter = Li_enter + 1
				Li_status[i][1] = 1
				Li_status[i][2] = 0 # By definition, is_enter = 1 means is_drag = 0


	# decide whether Li atoms drag Si out (that is, whether Li atoms participate in the process of breaking Si-Si chemical bonds)
	for i in range(0,len(Si_Si_bond_intralayer_0)):
		if Si_Si_bond_intralayer_0[i] in Si_Si_bond_intralayer_1:
			# this Si-Si bond is not broken, so nothing would happen
			pass
		else:
			# 有些层破坏的太严重了，就不要再考虑它了！
			Si_a = Si_Si_bond_intralayer_0[i][0]
			Si_b = Si_Si_bond_intralayer_0[i][1]
			whether_a_in_contact = 0
			whether_b_in_contact = 0
			for i in range(Si_layer_larger_than_me_is_contact,len(layer_of_Si)):
				if Si_a in layer_of_Si[i]:
					whether_a_in_contact = 1
					break
			for i in range(Si_layer_larger_than_me_is_contact,len(layer_of_Si)):
				if Si_b in layer_of_Si[i]:
					whether_b_in_contact = 1
					break
			if whether_a_in_contact == 1 and whether_b_in_contact == 1:
				Li_drag_a = []
				Li_drag_b = []
				for i in range(0,len(bond_info_Si_Li_0)):
					if Si_a in bond_info_Si_Li_0[i]:
						if bond_info_Si_Li_0[i][0]>amount_of_Si:
							Li_drag_a.append(bond_info_Si_Li_0[i][0])
						else:
							Li_drag_a.append(bond_info_Si_Li_0[i][1])
					if Si_b in bond_info_Si_Li_0[i]:
						if bond_info_Si_Li_0[i][0]>amount_of_Si:
							Li_drag_b.append(bond_info_Si_Li_0[i][0])
						else:
							Li_drag_b.append(bond_info_Si_Li_0[i][1])
				Li_drag_a_b = list(set(Li_drag_a).union(set(Li_drag_b)))
				Li_drag = Li_drag + len(Li_drag_a_b)
				for i in range(0,len(Li_drag_a_b)):
					Li_status[Li_drag_a_b[i]][2] = 1 # is_drag = 1



	# iteration progressing
	Si_Si_bond_intralayer_0 = [[0,0] for i in range(0,len(Si_Si_bond_intralayer_1))]
	for i in range(0,len(Si_Si_bond_intralayer_1)):
		for j in range(0,2):
			Si_Si_bond_intralayer_0[i][j] = Si_Si_bond_intralayer_1[i][j]
	Si_Si_bond_intralayer_1 = []

	bond_info_Si_Li_0 = [[0,0] for i in range(0,len(bond_info_Si_Li_1))]
	for i in range(0,len(bond_info_Si_Li_1)):
		for j in range(0,2):
			bond_info_Si_Li_0[i][j] = bond_info_Si_Li_1[i][j]
	bond_info_Si_Li_1 = []

	bond_info_Si_Si_0 = [[0,0] for i in range(0,len(bond_info_Si_Si_1))]
	for i in range(0,len(bond_info_Si_Si_1)):
		for j in range(0,2):
			bond_info_Si_Si_0[i][j] = bond_info_Si_Si_1[i][j]
	bond_info_Si_Si_1 = []

	if timestep%1000 == 0:
		for i in range(0,len(Li_z_coordinate_1)):
			Li_z_coordinate_0[i] = Li_z_coordinate_1[i]
		Li_z_coordinate_1 = [0]*(natom+1)

		for i in range(amount_of_Si+1,len(Li_status)):
			if Li_status[i][1] == 1:
				Li_status[i][2] = 0

		new_Li_enter = 0
		new_Li_drag = 0

		for i in range(amount_of_Si+1,len(Li_status)):
			if Li_status[i][1] == 1:
				new_Li_enter = new_Li_enter + 1
			elif Li_status[i][2] == 1:
				new_Li_drag = new_Li_drag + 1



	# output
	if timestep%1000 == 0:
		fout.write('%10d %6d %6d\n' %(timestep,new_Li_drag,new_Li_enter))








fin.close()
fout.close()
print('Exec ended!')
