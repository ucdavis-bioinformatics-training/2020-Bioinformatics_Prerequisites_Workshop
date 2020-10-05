# Jupyter, Pandas, and Seaborn
1. [Open this link](https://mybinder.org/v2/gh/ipython/ipython-in-depth/master?filepath=binder/Index.ipynb)
- binder is a web hosted temporary environment where you can run jupyter notebooks. The memory is limited to 2048MB and doesnt remember your libraries. 
- If you find this helpful I suggest you look into running it locally on your computer with the offical [Jupyter Notebook software](https://jupyter.org/)
2. Your screen should look something like this:
    ![](figures/71ea6011.png)
    
3. Click the + button to add a new cell and copy the following into the cell:
    ```
    import os
    os.system('wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-Bioinformatics_Prerequisites_Workshop/master/Intro_to_Python/JupyterPandasSeabornIntro.ipynb')
    ```
   
4. Then go File > Open and select the new .ipynb in the directory (JupyterPandasSeabornIntro.ipynb). 
The rest of the tutorial will take place in this notebook.