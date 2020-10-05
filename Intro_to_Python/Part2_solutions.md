1. Create a list of all even values from 0 to 100.
- What is the length of the list?
- What is the average of the list?  
- What is the sum of the 13 to the 17th elements in the list?
- What is the 16th element from the end?
- Replace the 23rd element with the value 200. Now what is the average?

    ```
    my_list = range(0,102,2)
    type(my_list)
    my_list = list(my_list)
    type(my_list)
    len(my_list)
    sum(my_list)/len(my_list)
    sum(my_list[13:18])
    my_list[-17] # note when working backwards the index starts at -1 rather the 0 like when forward
    my_list[22] = 200
    sum(my_list)/len(my_list)
    ```

2. Learning to use python it is good to have an experimental mentality to learn more about edge cases
- What is the empty list in terms of a boolean value?
- What is an empty dictionary in terms of a boolean value?
- What about a float value of 0?

    ```
    bool([])
    bool({})
    bool(0.0)
    ```


3. Create a string in python with some information about yourself (Name, background, etc.) 
- Take the string you created and make a list out of it where each word is an element.
- Now recreate the string with two spaces where there was originally one space. There are two ways, one using the list 
and one using the original string. Can you think of both?

    ```
    bio = "Hi my name is Keith Mitchell and I am a grad student in genetics/genomics and biostats."
    bio
    bio_list = bio.split()
    bio_list
    new_bio_string = '  '.join(bio_list)
    new_bio_string
    new_bio_string2 = bio.replace(' ', '  ')
    new_bio_string2
    ```


