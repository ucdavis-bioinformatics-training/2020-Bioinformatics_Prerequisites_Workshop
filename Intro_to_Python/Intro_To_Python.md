
## Intro to Python
<img src="./figures/python-logo2x.png" alt="Python" width="30%" align="right"/> 



### Schedule
*  Part 1: Intro To Python and first program (1 - 1:30pm)
*  Part 2: Basics of Python and data types (1:30pm-3pm)

3pm Break (15 minutes)

*  Part 3: Flow control (3:15 - 3:45pm) 
    *  Breakout session 3, Practice Flow Control (3:45 - 4pm)
*  Part 3: Continued, files, functions, documentation 
    *  Hands on Exercises 3
*  Part 4: Special topics (4 - 5pm)
    -  Biopython
    -  CSV package maybe https://docs.python.org/3/library/csv.html
    -  Plotting Seaborn / GGPlot





## Part 1: Intro to Python

### Why learn Python?


<img src="figures/programming_languages_recommended.png" alt="Python" width="70%" align="center"/>

https://businessoverbroadway.com/2019/01/13/programming-languages-most-used-and-recommended-by-data-scientists/


-  Python is extremely popular and widely used, especially for data science.
    +  Popular and getting more so in Bioinformatics, especially for building tools.
    +  For analysis, R (which you will learn later in the week) is arguably more useful currently due to the huge number of packages available from [Bioconductor](http://bioconductor.org/) and [CRAN](https://cran.r-project.org/). 
    +  The best option is to learn [Python, R, and bash](http://omgenomics.com/programming-languages/). A little of each will go a long way.
-  Freely available to [download](https://www.python.org/downloads/) for Windows, Linux, Mac OS X, etc.
-  Python is extremely versatile
    +  Used for a wide range of purposes from automating simple tasks to massive software projects with wide adoption [by many large companies.](https://realpython.com/world-class-companies-using-python/)
-  Installed on almost every Linux server.
-  Vast number of resources online: If you can Google for it you can learn how to do it.


------


### Goals for this section

-  Lower the barrier to entry by resolving basic "getting started" hurdles.
-  Provide tips and pointers on things to be aware of.
-  Provide a foundation in basic Python and some hands-on experience.
-  Give you some basic tools and recipes that you can build on in the future.
-  Get you hooked on programming and solving problems with Python.

-----

### Question
1)  Click "YES" if you have ever programmed in Python before.

2)  Click "YES" if you have ever programmed in any language before.

----- 

##  Background

###  What is a programming language and why do we need it?

Speaking to a computer in its native language is tedious and complicated. A programming language is a way for humans to describe a set of operations to a computer in a more abstract and understandable way. A helper program then translates our description of the operations into a set of instructions (machine code) for the computer to carry out.<br> 

Some day we may develop a programming language that allows us to communicate our instructions to the computer in our native language (Alexa, turn on the TV). Except for simple cases, this option doesn't exist yet, largely because human languages are complicated and instructions can be difficult to understand (even for other humans).


In order for the helper program to work properly, we need to use a concise language:
-  Well defined vocabulary for describing the basic set of supported operations.
-  Well defined set of Data Types that have a defined set of valid operations (add, subtract, etc).
-  Well defined syntax that leaves no ambiguity about when the computer should carry out each instruction.

Specifically in Python:

<img src="figures/python_interpreter.png" alt="PythonInterpreter" width="90%" align="center"/> <br><br>



-----

### A brief history of Python
-  Initially developed during the late 1980's by [Guido van Rossum](https://en.wikipedia.org/wiki/Guido_van_Rossum), BDFL until 2018.
-  First development version released in 1991. Version 1 released in 1994.
-  Python 2.0.0 released June, 2001
    -  Python 2.x end-of-life Jan 1, 2020.
    -  This version was so popular and widely used that many Bioinformatics programs were written using it. Some of these tools have been converted to support v3.x, others are in the process of being upgraded or have been abandoned and will stay on v2.x. The last Python 2.x release is still available for [download](https://www.python.org/downloads/release/python-2718/).
-  Python 3.x (December 2008) was a significant re-design and broke compatibility with some parts of v2.x.
-  The current version is 3.8.


------

### Interesting features of Python
-  High level: It hides a lot of the complicated details.
-  Interpreted: programs are compiled to byte code and run on a virtual machine instead of being compiled to native machine language
    -  This provides the option of running Python interactively, or writing Python scripts.
-  Garbage Collected: memory is allocated and freed for you automatically
-  Spaces matter in Python and are part of the language syntax. Be careful with copy/paste!
-  In Python, "Readability counts". 
    *  There is a style guide called [Python Enhancement Proposal 8 (PEP8)](https://www.python.org/dev/peps/pep-0008/) that documents and encourages consistent, readable code. If you plan to share your code with others, it is good practice to follow the style guide (or adopt the style used by the rest of the team). 
    *  These best practices are also known as writing "pythonic" or "idiomatic" python, [this guide](https://docs.python-guide.org/writing/style/) has more details. Try ```import this``` in your Python interpreter if you are a fan of programmer philosophy.

----

### Base Python and the extensive package ecosystem
-  Python has been extremely successful partly because it is modular and highly extensible. The core of Python is relatively small and compact, but this is supplemented with a large ["standard library"](https://docs.python.org/3/library/) that adds a large amount of additional functionality.
    +  Thousands of additional packages are available from the [PyPI](https://pypi.org/) repository.
    +  PythonPath variable
    +  Where do libraries live?
    +  Virtual Environments 
    +  Conflicts and package versions
        -  Virtual environments
        -  Conda


---------

### Writing Python Code

You can work directly on the interactive Python terminal, but eventually you might want a record of your work. One option is to edit your code in a text editor, then run it in the terminal. This is also a great way to test your code piece by piece as you work through developing or editing your program.


Lots of tools have been developed to help with writing computer code. Two popular options are [Visual Studio Code](https://code.visualstudio.com/), and [Atom Editor](https://atom.io/). Both work well and are available for free. Download and install either one of these if you want, or use a text editor you already have.

When working directly on the server, [Nano](https://www.nano-editor.org/), [Vim](https://www.vim.org/), and [Emacs](https://www.gnu.org/software/emacs/) are all popular editors.


-------


### Your first Python program

[Hello, World!](https://en.wikipedia.org/wiki/%22Hello,_World!%22_program) is traditionally the first program to write in any new programming language.


Connect to the server<br>
```bash
ssh <user>@tadpole.genomecenter.ucdavis.edu
```

Then follow this example:

```bash
export PATH=/share/biocore/shunter/opt/bin/:$PATH
alias ll='ls -alFh --color'

mkdir -p  /share/workshop/prereq_workshop/$USER/python

cd /share/workshop/prereq_workshop/$USER/python

ipython3
```

https://ipython.readthedocs.io/en/stable/index.html
<div class="output">
Python 3.8.2 (default, May  5 2020, 12:05:43) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.14.0 -- An enhanced Interactive Python. Type '?' for help.<br>

In [1]:   
</div>

Note that [IPython](https://ipython.org/) is an enhanced version of the Python shell. Among [many other things](https://ipython.readthedocs.io/en/stable/interactive/python-ipython-diff.html), it has better support for pasting Python code into the interactive terminal without messing up indentation. Some systems don't have ipython installed, in which case you can usually run "python" or "python3".


Use the interactive Python interpreter to run a short program:


```python
>>> print("Hello World!")
Hello World!
```


<h3><font color="red">Once you have successfully run `python3` and printed "Hello World!", mark "Yes" in zoom. Post questions or problems to the Slack channel.</font></h3>

Now we need to exit the Python interpreter. This can be done by typing ```exit()``` or **CTRL+D**. 


Next, lets try writing a "Hello World" Python script.

```bash
nano helloworld.py
```

Enter the same print statement we used above, press Ctrl+O to write the changes, then Ctrl+X to exit the file.

Check that the file you created is present and that your changes were saved successfully, then run the script:

```bash
ll

cat helloworld.py

python3 helloworld.py
```

<div class="output">
shunter@tadpole:python$ ll
total 5.0K
drwxrwxr-x  2 shunter shunter 2.0K Oct  4 20:11 ./
drwxrwx--- 26 shunter shunter 2.0K Oct  4 20:14 ../
-rw-rw-r--  1 shunter shunter   23 Oct  4 20:09 helloworld.py<br>

shunter@tadpole:python$ cat helloworld.py 
print("Hello, World!")

shunter@tadpole:python$ python3 helloworld.py
Hello, World!

</div>

<h3><font color="red">Once you have successfully created your Python script and run it, mark "Yes" in zoom. Post questions or problems to the Slack channel.</font></h3>



----

This is the end of the Intro section. Next we will talk about data types and do more hands-on exercies.

