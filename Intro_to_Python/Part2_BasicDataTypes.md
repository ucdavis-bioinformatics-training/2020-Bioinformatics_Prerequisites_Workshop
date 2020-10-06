# Part 2: Basics of Python and Basic Data Types

## Outline 
- [Basic data types](#datatypes)
- [Basic arithmetic](#arithmetic)
- [Comparisons](#comparisons)
- [Basic data structures](#datastructures)
- [Boolean Operations](#operations)
- [Range](#range)
- [String Formatting](#stringformat)
- [Group Exercise](#exercise)

A few notes before we get started:
- Learning all the nuances of python takes a long time! Our goal here is to introduce you to as many concepts as possible
but if you are serious about mastering python you will need to apply yourself beyond this introduction. 
- We bring up a lot of concepts to expose you to them but we encourage you to have a "scientific" mentality and highly
encourage you to continue testing the waters beyond these materials. 


Setup
```
cd /share/workshop/prereq_workshop/$USER/python
python3.8
```

---

# <a name="datatypes"></a> Basic Data Types: Integers, Floating-point numbers, booleans, strings.

### Integers
whole numbers, negative or positive

```
a = 1
type(a)
a
```
<div class="output">>>> a = 1
>>> type(a)
<class 'int'>
>>> a
1

</div>


### Floats
similar to decimal field

```
b = 1.2
type(b)
b
```

<div class="output">>>> b = 1.2
>>> type(b)
<class 'float'>
>>> b
1.2
</div>


### Booleans

"In computer science, the Boolean data type is a data type that has one of two possible values (usually denoted true 
and false) which is intended to represent the two truth values of logic and Boolean algebra. It is named after George 
Boole, who first defined an algebraic system of logic in the mid 19th century." [-wikipedia](https://en.wikipedia.org/wiki/Boolean_data_type)

```
c = True
d = False
type(c)
c
int(d)
# TODO float(d), long(d),
```

<div class="output">>>> c = True
>>> d = False
>>> type(c)
<class 'bool'>
>>> c
True
>>> int(d)
0
</div>

### Strings
```
my_string = "Hello world!"
type(my_string)
my_string
print(my_string)
my_string + my_string
my_string[1] # this will make more sense when we learn lists later
```


<div class="output">
>>> my_string = "Hello world!"
>>> type(my_string)
<class 'str'>
>>> my_string
'Hello world!'
>>> print(my_string)
Hello world!
>>> my_string + my_string
'Hello world!Hello world!'
>>> my_string[1] # this will make more sense when we learn lists later
'e'
</div>

---


<h3><font color="red">Challenge Questions:</font></h3>
1. What is an empty string and string in terms of a boolean value? (Hint: `bool(some_string)`)


---

# <a name="arithmetic"></a> Arithmetic: Adding, subtracting, multiplication, assignment arithmetic (assignment operators).
<img src="figures/ed9e3c89.png" alt="if flow" width="600px"/>

```
# This is a comment
type(a)
type(b)

# Addition
a + b
type(a+b)

# Subtraction
b - a
type(b-a)

# Division
a/b
type(a/b)

# Exponents
4**b
#or
pow(4,b)
type(a**b)

# Remainder 
4 % 3

# Absolute value
abs(22-32)

# Round, Floor, Ceiling
round(3.2)
int(3.2)
import math
math.ceil(3.2)
math.floor(3.7)
int(3.7)

```

<div class="output">
>>> # This is a comment
>>> type(a)
<class 'int'>
>>> type(b)
<class 'float'>
>>> 
>>> # Addition
>>> a + b
2.2
>>> type(a+b)
<class 'float'>
>>> 
>>> # Subtraction
>>> b - a
0.19999999999999996
>>> type(b-a)
<class 'float'>
>>> 
>>> # Division
>>> a/b
0.8333333333333334
>>> type(a/b)
<class 'float'>
>>> 
>>> # Exponents
>>> 4**b
5.278031643091577
>>> type(a**b)
<class 'float'>
</div>

Notice the error in the subtraction that produces 0.19999999999996 that should be 0.2
- [Python subtraction error](https://stackoverflow.com/questions/14120340/python-error-in-basic-subtraction)


## Assignment operators
<img src="figures/e652b5fc.png" alt="if flow" width="600px"/>

```
a += 1
a

a -= 1
a
```

<div class="output">
>>> a += 1
>>> a
2
>>> 
>>> a -= 1
>>> a
1
</div>

---


<h3><font color="red">Challenge Questions:</font></h3>
1. We have a value of 200 gene counts for a certain sample. One gene has a count of 7. 
    - What is the relative proportion of reads for this gene?
    - We find some novel information saying that all values have had `sqrt()` or `x**0.5` performed. Now what is the relative proportion?

---


# <a name="comparisons"></a> Comparisons: <, >, <= , >= , ==, !=
<img src="figures/de235da1.png" alt="if flow" width="600px"/>

```
1<1
1<2
2>1
1<=1
2>=1
1==1
0==1
```

<div class="output">
>>> 1<1
False
>>> 1<2
True
>>> 2>1
True
>>> 1<=1
True
>>> 2>=1
True
>>> 1==1
True
>>> 0==1
False

</div>

# <a name="datastructures"></a> Basic Data Structures: Lists, Sets, Tuples, Dictionaries.

## Lists
+  <img src="figures/44b15102.png" alt="if flow" width="600px"/>
+  <img src="figures/2b36588e.png" alt="if flow" width="600px"/>
+  <img src="figures/93dd9342.png" alt="if flow" width="600px"/>

```
my_list = [1,2,3,4,5,6]
type(my_list)

# getting the first element in a list (0 index)
my_list[0]

# getting the last element in a list
my_list[-1]
# OR
my_list[5]

# getting a range of the list
my_list[-3:]
my_list[1:3]
my_list[:3]

# some other features of integer lists
sum(my_list)
len(my_list)

my_list2 = [7,8,9]

# adding and subtracting lists
my_list+my_list2
my_list-my_list2


# string lists, double indexing
my_string_list = ['the', 'dog', 'says', 'woof']
my_string_list[0][1:]

# every x number of indexes
my_string_list[::2]
my_string_list[::3]

# join, splits, replace, 
" ".join(my_string_list)
list(my_string_list[0])
" ".join(my_string_list).split()
" ".join(my_string_list).split('s')
" ".join(my_string_list).replace(" ", "-")
```

<div class="output">
>>> my_list = [1,2,3,4,5,6]
>>> type(my_list)
<class 'list'>
>>> 
>>> # getting the first element in a list (0 index)
>>> my_list[0]
1
>>> 
>>> # getting the last element in a list
>>> my_list[-1]
6
>>> # OR
>>> my_list[5]
6
>>> 
>>> # getting a range of the list
>>> my_list[-3:]
[4, 5, 6]
>>> my_list[1:3]
[2, 3]
>>> my_list[:3]
[1, 2, 3]
>>> 
>>> # some other features of integer lists
>>> sum(my_list)
21
>>> len(my_list)
6
>>> 
>>> my_list2 = [7,8,9]
>>> 
>>> # adding and subtracting lists
>>> my_list+my_list2
[1, 2, 3, 4, 5, 6, 7, 8, 9]
>>> my_list-my_list2
Traceback (most recent call last):
  File 
TypeError: unsupported operand type(s) for -: 'list' and 'list'
>>> 
>>> 
>>> # string lists, double indexing
>>> my_string_list = ['the', 'dog', 'says', 'woof']
>>> my_string_list[0][1:]
'he'
>>> 
>>> # every x number of indexes
>>> my_string_list[::2]
['the', 'says']
>>> my_string_list[::3]
['the', 'woof']
>>> 
>>> # join, splits, replace, 
>>> " ".join(my_string_list)
'the dog says woof'
>>> list(my_string_list[0])
['t', 'h', 'e']
>>> " ".join(my_string_list).split()
['the', 'dog', 'says', 'woof']
>>> " ".join(my_string_list).split('s')
['the dog ', 'ay', ' woof']
>>> " ".join(my_string_list).replace(" ", "-")
'the-dog-says-woof'

</div>



### Some more list features: count, pop, append, reassignment
```
my_string_list = ['the', 'dog', 'says', 'woof']
my_string_list.append('the')
my_string_list.count('the')
my_string_list.pop()
my_string_list.count('the')
my_string_list.reverse()
my_string_list
```

<div class="output"> 
>>> my_string_list.append('the')
>>> my_string_list.count('the')
2
>>> my_string_list.pop()
'the'
>>> my_string_list.count('the')
1

</div>


---

<h3><font color="red">Challenge Questions:</font></h3>
1. We have the following dna sequence. We are curious about the adapters of the sequence (first and last 15 bps in the sequence.)
    ```
    some_dna = "ATCAATGCGCGCATACGATCAATGCGCGCATACGATCAATGCGCGCATACGGGTCCATACGCAATCAATGCGCGCATA"
    ```
   Can you tell me if they match?
   What is the GC content of each adapter?

---


## Tuples


```
my_tuple = (1,2,3,4,5,6)
type(my_tuple)

# getting the first element in a tuple (0 index)
my_tuple[0]

# getting the last element in a tuple
my_tuple[-1]
# OR
my_tuple[5]

# getting a range of the tuple
my_tuple[-3:]
my_tuple[1:3]
my_tuple[:3]

my_tuple2 = (7,8,9)

# adding and subtracting tuples
my_tuple+my_tuple2
my_tuple-my_tuple2
```

<div class="output">
>>> my_tuple = (1,2,3,4,5,6)
>>> type(my_tuple)
<class 'tuple'>
>>> 
>>> # getting the first element in a tuple (0 index)
>>> my_tuple[0]
1
>>> 
>>> # getting the last element in a tuple
>>> my_tuple[-1]
6
>>> # OR
>>> my_tuple[5]
6
>>> 
>>> # getting a range of the tuple
>>> my_tuple[-3:]
(4, 5, 6)
>>> my_tuple[1:3]
(2, 3)
>>> my_tuple[:3]
(1, 2, 3)
>>> 
>>> my_tuple2 = (7,8,9)
>>> 
>>> # adding and subtracting tuples
>>> my_tuple+my_tuple2
(1, 2, 3, 4, 5, 6, 7, 8, 9)
>>> my_tuple-my_tuple2
Traceback (most recent call last):
  File 
TypeError: unsupported operand type(s) for -: 'tuple' and 'tuple'
</div>

## Tuples vs Lists
<img src="figures/b783230c.png" alt="if flow" width="600px"/>


```
my_list[0] = 135
my_list
my_tuple[0] = 135
my_tuple

```
<div class="output">
>>> my_list = [1,2,3,4,5,6]
>>> my_list[0] = 135
>>> my_list
[135, 2, 3, 4, 5, 6]
>>> my_tuple[0] = 135
Traceback (most recent call last):
  File
TypeError: 'tuple' object does not support item assignment
>>> my_tuple
(1, 2, 3, 4, 5, 6)

</div>

There are quite a few differences between the two and some other features we talked about earlier for lists may work with tuples. 
I don't use tuples a ton but feel free to experiment and see what you can do with one but not the other. Check out the docs or see 
what features a class has using `list.__dict__.keys()` and `tuple.__dict__.keys()`. (These commands will make more sense later when 
we talk about attributes and dictionaries)

## Sets
Lets consider the events in my_set and my_set2 with respect to events A and B commonly used in denoting events in set theory. 
With respect to the pictures in the graph below where GREEN is the results objects in the set, consider the following events:

-  The relative complement is considered set difference in statistics. That is events in one group (B or `my_set2`) not in the other group (A or` my_set`).
    + <img src="https://latex.codecogs.com/gif.latex? \text{B}-\text{A}   " />    OR   <img src="https://latex.codecogs.com/gif.latex? \text{   my\_set2}-\text{my\_set   }" />

-  The relative complement is considered set difference in statistics. That is events in one group (A or `my_set`) not in the other group (B or `my_set2`).
    + <img src="https://latex.codecogs.com/gif.latex? \text{A}-\text{B}   " />    OR   <img src="https://latex.codecogs.com/gif.latex? \text{   my\_set}-\text{my\_set2   }" />

- In stats a union of two events is consider events in situation A and B. Similarly we consider `my_set` and `my_set2`.
    + <img src="https://latex.codecogs.com/gif.latex? \text{A}\cup\text{B}   " />    OR   <img src="https://latex.codecogs.com/gif.latex? \text{   my\_set}\cup\text{my\_set2}" />

- In stats a union of two events is consider events in situation A or B. Similarly we consider `my_set` or `my_set2`.
    + <img src="https://latex.codecogs.com/gif.latex? \text{A}\cup\text{B}   " />    OR   <img src="https://latex.codecogs.com/gif.latex? \text{   my\_set}\cup\text{my\_set2}" />


<img src="figures/sets.png" alt="if flow" width="600px"/>
<img src="figures/b028c7d8.png" alt="if flow" width="600px"/>



```
my_set = {1,2,3,4,5,5}
my_set
#OR 
my_set = set([1,2,3,4,5,5]) 
my_set
type(my_set)


# Can we index sets like we do lists? No
my_set[0]


my_set2 = {5,7,8,9}

# intersection and union of sets etc.

# union
my_set|my_set2
# intersection
my_set&my_set2
# subtracting sets
my_set-my_set2
my_set2-my_set

```


<div class="output">

>>> my_set = {1,2,3,4,5,5}
>>> my_set
{1, 2, 3, 4, 5}
>>> #OR 
>>> my_set = set([1,2,3,4,5,5]) 
>>> my_set
{1, 2, 3, 4, 5}
>>> type(my_set)
<class 'set'>
>>> 
>>> 
>>> # Can we index sets like we do lists? No
>>> my_set[0]
Traceback (most recent call last):
  File 
TypeError: 'set' object is not subscriptable
>>> 
>>> 
>>> my_set2 = {5,7,8,9}
>>> 
>>> # intersection and union of sets etc.
>>> 
>>> # union
>>> my_set|my_set2
{1, 2, 3, 4, 5, 7, 8, 9}
>>> # intersection
>>> my_set&my_set2
{5}
>>> # subtracting sets
>>> my_set-my_set2
{1, 2, 3, 4}
>>> my_set2-my_set
{8, 9, 7}
</div>

---

<h3><font color="red">Challenge Questions:</font></h3>
1. What is the union, intersection, and difference between the following two sets of amino acids in protein A and B?
   ```
   A = {'G', 'F', 'P', 'S', 'T', 'Y', 'C', 'M', 'M', 'Q'}
   B = { 'Y', 'C', 'Q', 'G', 'F', 'P', 'T', 'W', 'E', 'N'}
   ```

---


## Dictionaries
- A list of codes, terms, keys, etc., and their meanings, used by a computer program or system.
- Dictionary values are pointed too by the keys. Values can be anything from int, float, and bool to lists, tuples, and dictionaries.
    + <img src="figures/abda21dc.png" alt="if flow" width="600px"/>


```
my_dict = {'a':1, 'b':2, 'c':3}
type(my_dict)

my_dict.keys()
my_dict.values()
my_dict['a']

# adding a new key value set to the dictionary
my_dict['d'] = 4
my_dict
```

<div class="output">
>>> my_dict = {'a':1, 'b':2, 'c':3}
>>> type(my_dict)
<class 'dict'>
>>> 
>>> my_dict.keys()
dict_keys(['a', 'b', 'c'])
>>> my_dict.values()
dict_values([1, 2, 3])
>>> my_dict['a']
1
>>> my_dict['d']
{'a':1, 'b':2, 'c':3, 'd':4}
</div>


# <a name="operations"></a> Boolean Operations: and, or, not 
<img src="figures/6d14333b.png" alt="if flow" width="600px"/>
<img src="figures/33ec5c8c.png" alt="if flow" width="600px"/>

```
1 and 1
0 and 1
0 or 1 
not 1
not 0
not True
not False
True and False
True or False
None and None
None or None


bool(None) # we will see why this is important later when we talk about if statements in part3
bool(0)
bool(1)
```

<div class="output">
>>> 1 and 1
1
>>> 0 and 1
0
>>> 0 or 1 
1
>>> not 1
False
>>> not 0
True
>>> not True
False
>>> not False
True
>>> True and False
False
>>> True or False
True
>>> None and None
>>> None or None
>>> 
>>> 
>>> bool(None) # we will see why this is important later when we talk about if statements in part3
False
>>> bool(0)
False
>>> bool(1)
True

</div>



# <a name="range"></a> Range
The `range()` function returns a sequence of numbers, starting from 0 by default, and increments by 1 (by default), and stops before a specified number.

<img src="figures/d9454cac.png" alt="if flow" width="600px"/>

Much more info is available on this [here](https://docs.python.org/3/library/stdtypes.html#range).

```
type(range(0,10))

# first arg is start, 2nd is finish, 3rd is gap size (default is 1)
list(range(0,10))
list(range(0,10,5))
```

<div class="output">
<class 'range'>
>>> # first arg is start, 2nd is finish, 3rd is gap size (default is 1)
>>> list(range(0,10))
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
>>> list(range(0,10,5))
[0, 5]
</div>



# <a name="stringformatting"></a> String Formatting

This is just a brief intro. There are lots of features to this topic but this seems to solve 95% of tasks/cases for me.
If you interested in more info regarding the other featrures checkout [Pyformat](https://pyformat.info/).

<img src="figures/1a9dfa16.png" alt="if flow" width="600px"/>

```
extra = "Goodbye."
print("It was nice to meet you. %s" %extra)
print("It was nice to meet you. %s %s" %(extra,extra))
print(f"It was nice to meet you. {extra}")
# or
print("It was nice to meet you. {}".format(extra))
```

<div class="output">
>>> extra = "Goodbye."
>>> print("It was nice to meet you. %s" %extra)
It was nice to meet you. Goodbye.
>>> print("It was nice to meet you. %s %s" %(extra,extra))
It was nice to meet you. Goodbye. Goodbye.
>>> print(f"It was nice to meet you. {extra}")
It was nice to meet you. Goodbye.
>>> print("It was nice to meet you. {}".format(extra))
It was nice to meet you. Goodbye.
</div>


---

# <a name="exercise" style="color:red"></a> Group Exercises (~30 mins)
1. Create a list of all even values from 0 to 100.
- What is the length of the list?
- What is the average of the list?  
- What is the sum of the 13 to the 17th elements in the list?
- What is the 16th element from the end?
- Replace the 23rd element with the value 200. Now what is the average?


2. Learning to use python it is good to have an experimental mentality to learn more about edge cases
- What is the empty list in terms of a boolean value?
- What is an empty dictionary in terms of a boolean value?
- What about a float value of 0?


3. Create a string in python with some information about yourself (Name, background, etc.) 
- Take the string you created and make a list out of it where each word is an element.
- Now recreate the string with two spaces where there was originally one space. There are two ways, one using the list 
and one using the original string. Can you think of both?




---


