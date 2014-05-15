#test.py
#Copyright (C) 2014  Mathias Christensen Cronjaeger KitKat2.0@gmail.com

#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

class testClass(object):
	
	def __init__(self,x):
		self.foo = x
		print "Initialized"
  

def fib(n):
    if n == 0 or n==1:
        return 1
    else:
        return fib(n-1) + fib(n-2)

def listCalc(n):
    l = []
    listBuild(l,n,n,0)
    return l

def listBuild(l,n,N,length):
#    print n,N,l
    if length == N-1:
        return l.append(N)
    else:
        l.append(n)
        listBuild(l,n-1,N,length+1)