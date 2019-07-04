# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import mylist;

nested_list = mylist.Results

def identity(x):
    return x

def my_max(sequence, key_func=lambda x: x[4]):
    """
    Return the maximum element of a sequence.
    key_func is an optional one-argument ordering function.
    """
    if not sequence:
        raise ValueError('empty sequence')

    if not key_func:
        key_func = identity

    maximum = sequence[0]

    for item in sequence:
        # Ask the key func which property to compare
        if key_func(item) != 1:
            if key_func(item) > key_func(maximum):
                maximum = item
        '''if key_func(item) == 1:
            print("edge(", end ="")
            print(item[0], end=", ")
            print(item[1], end=", ")            
            print(item[2], end=", ")            
            print(item[3], end="")
            print(").")'''



    return maximum

'''print(*my_max(nested_list), sep = ", ")'''
mymax = my_max(nested_list)
print("edge(", end ="")
print(mymax[0], end=", ")
print(mymax[1], end=", ")            
print(mymax[2], end=", ")            
print(mymax[3], end="")
print(").")
'''my_max(nested_list)'''