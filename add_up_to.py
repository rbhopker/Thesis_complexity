#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:48:08 2022

@author: ricardobortothopker
"""

# Python3 program to find out all
# combinations of positive
# numbers that add upto given number

# arr - array to store the combination
# index - next location in array
# num - given number
# reducedNum - reduced number
def findCombinationsUtil(arr, index, num,
							reducedNum,tot_arr):

	# Base condition
    if (reducedNum < 0):
        return

	# If combination is
	# found, print it
    if (reducedNum == 0):
        # for i in range(index):
        #     print(arr[i], end = " ")
        # print("")
        result = list(filter(lambda val: val !=  1, arr[:index].copy()))
        # result = list(result)
        tot_arr.append(result)
        return

	# Find the previous number stored in arr[].
	# It helps in maintaining increasing order
    prev = 1 if(index == 0) else arr[index - 1]

	# note loop starts from previous
	# number i.e. at array location
	# index - 1
    for k in range(prev, num + 1):
		
		# next element of array is k
        arr[index] = k

		# call recursively with
		# reduced number
        findCombinationsUtil(arr, index + 1, num,
								reducedNum - k,tot_arr)

# Function to find out all
# combinations of positive numbers
# that add upto given number.
# It uses findCombinationsUtil()
def findCombinations(n):
	
	# array to store the combinations
	# It can contain max n elements
    arr = [0] * n
    out_arr =[]
	# find all combinations
    findCombinationsUtil(arr, 0, n, n,out_arr)
    return out_arr[1:-1]
def abstractions(n):
    lvl = findCombinations(n)
    out = []
    for i in lvl:
        temp = []
        count = 0
        for j in i:
            abst1 = []
            for k in range(j):
                abst1.append(count)
                count+=1
            temp.append(abst1)
        out.append(temp)
    return out
            

# Driver code
n = 13;
my_arr = findCombinations(n)
abstr = abstractions(n)
# This code is contributed by mits
