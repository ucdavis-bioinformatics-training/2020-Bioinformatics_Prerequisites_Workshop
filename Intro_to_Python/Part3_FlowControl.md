# Part 3: Flow Control

## Outline 
- [If, else if, else](#ifs)
- [For loops](#fors)
- [List comprehension](#listcomprehension)
- [While, Iterators](#while)
- [Break, Continue, Pass, Try, Except](#stoppers)

TODO 
talk about None, string formatting (%s and f)
- [Simple Functions](#functions)
- [Group Exercise](#exercise)


---

# <a name="ifs"></a> If, else if, else.

```
if True:
    print("Hi there!")

if 22:
    print("Hello world")

if True or False:
    print("Binary operators work here.")


if 0:
    print("This won't be printed")
elif 1:
    print("Any value other then 0 is treated as True")
else:
    print("This won't be printed")

# Be sure to include a couple of enters here so it all executes!

```
<div class="output">

>>> if True:
...     print("Hi there!")
... 
Hi there!
>>> if 22:
...     print("Hello world")
... 
Hello world
>>> if True or False:
...     print("Binary operators work here.")
... 
Binary operators work here.
>>> 
>>> if 0:
...     print("This won't be printed")
... elif 1:
...     print("Any value other then 0 is treated as True")
... else:
...     print("This won't be printed")
... 
Any value other then 0 is treated as True
>>> # Be sure to include a couple of enters here so it all executes!
</div>



# <a name="fors"></a> For loops

## Iterating through lists

```
wizard_list = ["Harry", "Ron", "Hermione", "Fred", "George"]
for name in wizard_list:
    print(name + " says Wingardium Leviosa")
```


<div class="output">
>>> wizard_list = ["Harry", "Ron", "Hermione", "Fred", "George"]
>>> for name in wizard_list:
...     print(name + " says Wingardium Leviosa")
... 
Harry says Wingardium Leviosa
Ron says Wingardium Leviosa
Hermione says Wingardium Leviosa
Fred says Wingardium Leviosa
George says Wingardium Leviosa
</div>


## Iterating through dictionaries

```
wizard_dict = {"Harry": "Lumos", "Ron": "Alohomora", "Hermione": "Wingardium Leviosa", "Fred": "Riddikulus", "George": "Sectumsempra"} 
for wizard in wizard_dict.keys():
    print(wizard + " says " + wizard_dict[wizard]) 


for spell in wizard_dict.values():
    print(f"A wizard says {spell}")
```

<div class="output">
>>> wizard_dict = {"Harry": "Lumos", "Ron": "Alohomora", "Hermione": "Wingardium Leviosa", "Fred": "Riddikulus", "George": "Sectumsempra"} 
>>> for wizard in wizard_dict.keys():
...     print(wizard + " says " + wizard_dict[wizard]) 
... 
Harry says Lumos
Ron says Alohomora
Hermione says Wingardium Leviosa
Fred says Riddikulus
George says Sectumsempra
>>> 
>>> for spell in wizard_dict.values():
...     print(f"A wizard says {spell}")
... 
A wizard says Lumos
A wizard says Alohomora
A wizard says Wingardium Leviosa
A wizard says Riddikulus
A wizard says Sectumsempra
</div>



# <a name="listcomprehension"></a> List Comprehension

```
wizard_dict = {"Harry Potter": "Lumos", "Ron Weasley": "Alohomora", "Hermione Granger": "Wingardium Leviosa", \
"Fred Weasley": "Riddikulus", "George Weasley": "Sectumsempra"} 
weasley_list = [wizard for wizard in wizard_dict.keys() if "Weasley" in wizard]
weasley_list
h_list = [wizard if "H" in wizard else "Voldemort" for wizard in wizard_dict.keys()]
h_list
```

<div class="output">
>>> wizard_dict = {"Harry Potter": "Lumos", "Ron Weasley": "Alohomora", "Hermione Granger": "Wingardium Leviosa", \
... "Fred Weasley": "Riddikulus", "George Weasley": "Sectumsempra"} 
>>> weasley_list = [wizard for wizard in wizard_dict.keys() if "Weasley" in wizard]
>>> weasley_list
['Ron Weasley', 'Fred Weasley', 'George Weasley']
>>> h_list = [wizard if "H" in wizard else "Voldemort" for wizard in wizard_dict.keys()]
>>> h_list
['Harry Potter', 'Voldemort', 'Hermione Granger', 'Voldemort', 'Voldemort']
</div>


# <a name="while"></a> While and Iterators

## 

```
h_iter = iter(h_list)
cont = True
while cont:
    if "Voldemort" == next(h_iter):
        cont = False
    else:
        print(next(h_iter))

h_iter = iter(h_list)
cont = True
while cont:
    value = next(h_iter)
    if "Voldemort" == value:
        cont = False
    else:
        print(value)
```

<div class="output">
>>> 
>>> h_iter = iter(h_list)
>>> cont = True
>>> while cont:
...     if "Voldemort" == next(h_iter):
...         cont = False
...     else:
...         print(next(h_iter))
... 
Voldemort
Voldemort
>>> h_iter = iter(h_list)
>>> cont = True
>>> while cont:
...     value = next(h_iter)
...     if "Voldemort" == value:
...         cont = False
...     else:
...         print(value)
... 
Harry Potter
</div>

# <a name="stoppers"></a> Break, continue, pass, try, except

## Continue vs pass vs break

```
a = [0, 1, 2]
for element in a:
    if not element:
        pass
    print(element)


for element in a:
    if not element:
        continue
    print(element)


for element in a:
    if not element:
        break
    print(element)
```


<div class="output">
>>> a = [0, 1, 2]
>>> for element in a:
...     if not element:
...         pass
...     print(element)
... 
0
1
2
>>> for element in a:
...     if not element:
...         continue
...     print(element)
... 
1
2
>>> for element in a:
...     if not element:
...         break
...     print(element)
... 
</div>



# <a name="functions"></a> Simple Functions: parameters, scope, returning, passing

## Lets define a simple function

```
def simple_function(
```




# <a name="exercise"></a> Group Exercises (~30 mins)
1. 
4. 






<div class="output">
</div>

<div class="output">
</div>

<div class="output">

</div>

*  Part 3: Flow control
        -  range
        -  For Loops, while etc
        -  Iterators - what is an iterator, how does it work?
        -  Maybe list comprehension, lambda etc?
    -  Basic file read and write
    -  Simple functions
        -  What is "Scope" and how does it work?
        -  defining variables, passing variables, returning variables
    -  Hands on Exercises 3