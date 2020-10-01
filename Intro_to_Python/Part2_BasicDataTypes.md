# Part 2: Basics of Python and Basic Data Types

## Outline 
- [Basic data types](#datatypes)
- [Basic arithmetic](#arithmetic)
- [Comparisons](#comparisons)
- [Boolean Operations](#operations)
- [Basic data structures](#datastructures)


---

# <a name="datatypes"></a> Basic Data Types: Integers, Floating-point numbers, booleans, strings.

## Integers

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


## Floats

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


## Booleans

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





# <a name="arithmetic"></a> Arithmetic: Adding, subtracting, multiplication, etc.

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
type(a**b)

# TODO (//, abs, pow(x,y))

```

<div class="output">>>> # This is a comment
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


## A few other things

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

# <a name="comparisons"></a> Comparisons: <, >, <= , >= , ==,

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

# <a name="dattastructures"></a> Basic Data Structures: Lists, Sets, Tuples, Dictionaries.

## Lists
```
my_list = [1,2,3,4,5,6]
type(my_list)

# getting the first element in a list (0 index)
my_list[0]

# getting the last element in a list
my_list[-1]
# OR
my_list[5]

my_list2 = [7,8,9]

# adding and subtracting lists
my_list+my_list2
my_list-my_list2
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
>>> my_list2 = [7,8,9]
>>> 
>>> # adding and subtracting lists
>>> my_list+my_list2
[1, 2, 3, 4, 5, 6, 7, 8, 9]
>>> my_list-my_list2
Traceback (most recent call last):
  File
TypeError: unsupported operand type(s) for -: 'list' and 'list'
</div>


## Tuples


```

```

## Sets



## Dictionaries




# <a name="operations"></a> Boolean Operations: and, or, not 






-  [Basic data types](https://docs.python.org/3/library/datatypes.html)
    +  Integers, Floating-point numbers, boolean, 
        - abs(), 
    +  Strings
        -  Single and double quotes, escape characters, etc
        -  print()
        -  inserting variables 
        -  
    +  Comments (#)




*  Part 2: Basics of Python and Basic Data Types
    -  [Basic data types](https://docs.python.org/3/library/datatypes.html)
        +  Integers, Floating-point numbers, boolean, 
            - abs(), 
        +  Strings
            -  Single and double quotes, escape characters, etc
            -  print()
            -  inserting variables 
            -  
        +  Comments (#)

*  Part 2: Hands on exercises (1:30 - 2pm)
    -  TODO
*  Part 2: Continued more advanced data types
    -  Lists and Tuples
        -  Zero-based index, variables vs objects
    -  Dictionaries
    -  Maybe Counter() and OrderedDict() from collections (https://docs.python.org/2/library/collections.html)
    -  Hands on Exercises 2
*  Part 3: Flow control
        -  if statements
        -  range
        -  break, continue, pass
        -  For Loops, while etc
        -  Iterators - what is an iterator, how does it work?
        -  Maybe list comprehension, lambda etc?
    -  Basic file read and write
    -  Simple functions
        -  What is "Scope" and how does it work?
        -  defining variables, passing variables, returning variables
    -  Hands on Exercises 3
