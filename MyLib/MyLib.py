# -*- coding: utf-8 -*-

import numpy as np
import os
import time
import math
import sys

#线性内插法
def Linear_Interpolation(x1,y1,x2,y2,xx):
	return (y2-y1)/(x2-x1)*xx+(x2*y1-x1*y2)/(x2-x1)

#指数分区，用于在对数分布中适当取点，radix常取1.1，1.2，1.3
def Create_Partition(m_max,m_min,radix=1.1):
	if m_min==0:
		base=1
	elif m_min>0:
		base=m_min
	else:
		print "MyLib_Create_Partition error, m_min>0"
		sys.exit()
	index_max=int(math.log(m_max/base,radix))+1
	if m_min==0:
		index_min=0
	else:
		index_min=int(math.log(m_min/base,radix))

	partition=[]
	for i in range(index_min,index_max+1):
		partition.append(base*radix**i)
	return partition

def Create_Standard_Partition(x_min,x_max):
	if x_min==0:
		x_min=1
		partition=[0]
	elif x_min>0:
		partition=[]

	log_min=int(math.log(x_min,10))
	log_max=int(math.log(x_max,10))+1
	
	for i in range(log_min,log_max):
		for j in range(9):
			partition.append((j+1)*10**i)
	partition.append(10**log_max)
	return partition

def Fib_Search(x,arr): 
    n=len(arr)
    fibMMm2 = 0 # (m-2)'th Fibonacci No. 
    fibMMm1 = 1 # (m-1)'th Fibonacci No. 
    fibM = fibMMm2 + fibMMm1 # m'th Fibonacci 
  
    # Fibonacci Number greater than or equal to n  
    while (fibM < n): 
        fibMMm2 = fibMMm1 
        fibMMm1 = fibM 
        fibM = fibMMm2 + fibMMm1 
  
    # Marks the eliminated range from front 
    offset = -1; 
  
    # while there are elements to be inspected. 
    # Note that we compare arr[fibMm2] with x. 
    # When fibM becomes 1, fibMm2 becomes 0  
    while (fibM > 1): 
          
        # Check if fibMm2 is a valid location 
        i = min(offset+fibMMm2, n-1) 
  
        # If x is greater than the value at  
        # index fibMm2, cut the subarray array  
        # from offset to i  
        if (arr[i] < x): 
            fibM = fibMMm1 
            fibMMm1 = fibMMm2 
            fibMMm2 = fibM - fibMMm1 
            offset = i 
  
        # If x is greater than the value at  
        # index fibMm2, cut the subarray  
        # after i+1 
        elif (arr[i] > x): 
            fibM = fibMMm2 
            fibMMm1 = fibMMm1 - fibMMm2 
            fibMMm2 = fibM - fibMMm1 
  
        # element found. return index  
        else : 
            return i 
  
    # comparing the last element with x */ 
    if(fibMMm1 and arr[offset+1] == x): 
        return offset+1; 
  
    # element not found. return -1  
    return -1












