from brial import *
import copy
import pdb
import xdrlib ,sys
import random
import time
import gc



lanesize=8

varnames = sum ([[ v+ str (i) for i in range (200)] for v in 'ck'], [])
Keccak = BooleanPolynomialRing ( len ( varnames ), varnames , order = 'deglex')
Cvar, Kvar = list (Keccak.gens())[:200], list (Keccak.gens())[200:]



def invchi(state):

	tempstate=[[] for i in range(25)]
	for i in range(5):
		for j in range(5):
			for k in range(lanesize):
				tempstate[5*i+j].append(state[5*i+j][k]+ (state[5*i+(j+1)%5][k]+1) * (state[5*i+(j+2)%5][k]+ (state[5*i+(j+3)%5][k]+1)*state[5*i+(j+4)%5][k] ))

	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]

def invpi(state):
	tempstate=[[] for i in range(25)]
	for i in range(25):
		y=floor(i/5)
		x=i%5
		x1=y
		y1=(2*x+3*y)%5
		temp=5*y1+x1
		for j in range(lanesize):
			tempstate[i].append(state[temp][j])
	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]

def invrio(state):
	rot=[0,1,62,28,27,36,44,6,55,20,3,10,43,25,39,41,45,15,21,8,18,2,61,56,14]
	tempstate=[[] for i in range(25)]
	for i in range(25):
		for j in range(lanesize):
			tempstate[i].append(state[i][(j+rot[i])%lanesize])

	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]


def invtheta(state):

	xor=[0xe8,0xca,0x3e,0x5b,0x44]
	columnstate=[[] for i in range(5)]
	sumstate = [[Keccak(0) for j in range(lanesize)] for i in range(5)]
	for i in range(5):
		for j in range(5):
			for k in range(lanesize):
				columnstate[i].append(state[5*j+i][k]+state[5*((j+1)%5)+i][k]+state[5*((j+2)%5)+i][k]+state[5*((j+3)%5)+i][k]+state[5*((j+4)%5)+i][k])
	
	for i in range(5):
		for j in range(lanesize):
			#sumstate[i].append(Keccak(0))
			for m in range(5):
				for n in range(lanesize):
					if((xor[m]>>n)&0x1):
						sumstate[i][j] += columnstate[(m+i)%5][(j+n)%lanesize]
	for i in range(5):
		for j in range(5):
			for k in range(lanesize):
				state[5*i+j][k] += sumstate[j][k]


			
def theta(state):
	tempstate=[[] for i in range(25)]
	for i in range(25):
		for j in range(lanesize):
			tempstate[i].append(state[i][j])
			for k in range(5):
				tempstate[i][j]+=state[(i%5+5-1)%5+5*k][j]+state[(i%5+1+5)%5+5*k][(j-1+lanesize)%lanesize]

	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]



def rio(state):
	rot=[0,1,62,28,27,36,44,6,55,20,3,10,43,25,39,41,45,15,21,8,18,2,61,56,14]
	tempstate=[[] for i in range(25)]
	for i in range(25):
		for j in range(lanesize):
			tempstate[i].append(state[i][(j-rot[i]+lanesize)%lanesize])

	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]

def pi(state):
	tempstate=[[] for i in range(25)]
	for i in range(25):
		y=floor(i/5)
		x=i%5
		x1=y
		y1=(2*x+3*y)%5
		temp=5*y1+x1
		for j in range(lanesize):
			tempstate[temp].append(state[i][j])
	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]

def chi(state):
	tempstate=[[] for i in range(25)]
	for i in range(5):
		for j in range(5):
			for k in range(lanesize):
				tempstate[5*i+j].append(state[5*i+j][k]+(state[5*i+(j+1)%5][k]+1)*state[5*i+(j+2)%5][k])

	for i in range(25):
		for j in range(lanesize):
			state[i][j]=tempstate[i][j]




begin_time = time.time()

state=[[] for i in range(25)]
			
for i in range(200):
	state[i/lanesize].append(Kvar[i]+Cvar[i])


invchi(state)
invpi(state)
invrio(state)
invtheta(state)

#invchi()
#state[0][0] = state[0][0]+ (state[1][0]+1) * (state[2][0]+(state[3][0]+1) * state[4][0])
x00 = (state[3][0]+1) * state[4][0]

#coefficients of (state[3][0]+1) * state[4][0]
x = set([])
mon = x00.monomials()
for i in range(0,len(mon)):
	for ci in Cvar:
		if(mon[i]/Keccak(ci)!=0):
			mon[i] = mon[i]/Keccak(ci)
	x.add(mon[i])
del mon
gc.collect()

#coefficients of state[2][0]+(state[3][0]+1) * state[4][0]
mon = state[2][0].monomials()
for i in range(0,len(mon)):
	for ci in Cvar:
		if(mon[i]/Keccak(ci)!=0):
			mon[i] = mon[i]/Keccak(ci)
	x.add(mon[i])
del mon
gc.collect()



#state[1][0]+1
xx = set([])
mon = (state[1][0]+1).monomials()
for i in range(0,len(mon)):
	for ci in Cvar:
		if(mon[i]/Keccak(ci)!=0):
			mon[i] = mon[i]/Keccak(ci)
	xx.add(mon[i])
del mon
gc.collect()



#coefficients of (state[1][0]+1) * (state[2][0]+(state[3][0]+1) * state[4][0])
sumset = set([])
for i in x:
	for j in xx:
		sumset.add(i*j)
del x
del xx
gc.collect()


mon = state[0][0].monomials()

for i in range(0,len(mon)):
	for ci in Cvar:
		if(mon[i]/Keccak(ci)!=0):
			mon[i] = mon[i]/Keccak(ci)
	sumset.add(mon[i])
del mon
gc.collect()

print '\n','\n',len(sumset),'\n','\n'





