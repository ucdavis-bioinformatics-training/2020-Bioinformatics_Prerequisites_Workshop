# Controlling Environments with Conda

## What is conda and why do we care about it?
- Conda is a cross platform (Windows, Linux, MacOSx) environment management software. 
- Works for any language, version, and library as specified by the user. 
- Allows for easy control and maintenance of multiple environments with different specifications.
- Allows for software packages and their dependencies to be easily loaded "automatically" via channels such as Bioconda, which specifically hosts a variety of bioinformatics software.  

    ![](./conda.png)

     - Installability of 98 randomly selected published software tools across 22 life-science journals over a span of 15 years. Error bars, where present, indicate SEM. (A) Pie chart showing the percentage of tools with various levels of installability. (B) A pie chart showing the proportion of evaluated tools that required no deviation from the documented installation procedure. (C) Tools that require no manual intervention (pass automatic installation test) exhibit decreased installation time. (D) Tools installed exhibit increased citation per year compared with tools that were not installed (Kruskal-Wallis, p-value = 0.035). (E) Tools that are easy to install include a decreased portion of undocumented commands (Not Installed versus Easy Install: Mann-Whitney U test, p-value = 0.01, Easy Install versus Complex Install: Mann-Whitney U test, p-value = 8.3 × 10 −8 ). (F) Tools available in well-maintained package managers such as Bioconda were always installable, whereas tools not shipped via package managers were prone to problems in 32% of the studied cases. SEM, standard error of the mean. https://doi.org/10.1371/journal.pbio.3000333.g002

## Getting and using conda

### Download miniconda

Anaconda is the name of the software where we get conda from. However, it is quite large, so we are going to instead download a smaller version called miniconda. Go to the [miniconda website](https://docs.conda.io/en/latest/miniconda.html) and get the link for the Linux 64 bit python 3.8 installer. Then download the installer:

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Run the installer:

	sh Miniconda3-latest-Linux-x86_64.sh

Press Enter to continue the installer and then scroll through the EULA using the spacebar. Then type "yes" to accept the terms. Now, we want to change the default install location, so enter "/share/workshop/prereq_workshop/$USER/software/miniconda3":

<div class="output">Do you accept the license terms? [yes|no]
[no] >>> yes       

Miniconda3 will now be installed into this location:
/home/joshi/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/joshi/miniconda3] >>> /share/workshop/prereq_workshop/$USER/software/miniconda3
PREFIX=/share/workshop/prereq_workshop/joshi/software/miniconda3
Unpacking payload ...
Collecting package metadata (current_repodata.json): done                                                          
Solving environment: done
</div>

You will see a bunch of messages about packages being installed, and then it will ask you if you want to initialize Miniconda3. Type "yes". The initialization just creates a file called "\~/.bashrc" which we will need to source to actually initialize:

	source ~/.bashrc

You will see that your prompt changes to have "(base)" at the beginning.

### Installing software

Conda is useful for software installation because you can create an "environment" for each piece of software and that environment will have exactly the right packages it needs with the exact versions it needs. Different software will have different requirements in terms of which packages and versions it needs, so it is useful to compartmentalize those. Now let's create an environment and install some software. Here we want to install HTStream, a tool for High Throughput Sequencing Read Processing.

[HTSream Homepage](https://github.com/s4hts/HTStream)

[HTStream Bioconda Documentation](https://bioconda.github.io/recipes/htstream/README.html)

First, let's search for htstream:

	conda search htstream

<div class="output">(base) joshi@tadpole:/share/workshop/prereq_workshop/joshi/software$ conda search htstream
Loading channels: done
# Name                       Version           Build  Channel             
htstream                       1.0.0      hf4b34a0_0  bioconda            
htstream                       1.1.0      h5ca1c16_0  bioconda            
htstream                       1.2.0      h5ca1c16_0  bioconda            
htstream                       1.3.1      h5ca1c16_0  bioconda            
htstream                       1.3.2      h5ca1c16_0  bioconda            
</div>

You see that there are multiple versions of htstream you could install. And that the software is available on the "bioconda" channel. If you run "conda info" you can see information about the install including the channels that are searched.

<div class="output">(base) joshi@tadpole:/share/workshop/prereq_workshop/joshi/software$ conda info

     active environment : base
    active env location : /share/workshop/prereq_workshop/joshi/software/miniconda3
            shell level : 1
       user config file : /home/joshi/.condarc
 populated config files : /home/joshi/.condarc
          conda version : 4.8.3
    conda-build version : not installed
         python version : 3.8.3.final.0
       virtual packages : __glibc=2.23
       base environment : /share/workshop/prereq_workshop/joshi/software/miniconda3  (writable)
           channel URLs : https://conda.anaconda.org/conda-forge/linux-64
                          https://conda.anaconda.org/conda-forge/noarch
                          https://conda.anaconda.org/bioconda/linux-64
                          https://conda.anaconda.org/bioconda/noarch
                          https://repo.anaconda.com/pkgs/main/linux-64
                          https://repo.anaconda.com/pkgs/main/noarch
                          https://repo.anaconda.com/pkgs/r/linux-64
                          https://repo.anaconda.com/pkgs/r/noarch
                          https://conda.anaconda.org/r/linux-64
                          https://conda.anaconda.org/r/noarch
                          https://conda.anaconda.org/etetoolkit/linux-64
                          https://conda.anaconda.org/etetoolkit/noarch
          package cache : /share/workshop/prereq_workshop/joshi/software/miniconda3/pkgs
                          /home/joshi/.conda/pkgs
       envs directories : /share/workshop/prereq_workshop/joshi/software/miniconda3/envs
                          /home/joshi/.conda/envs
               platform : linux-64
             user-agent : conda/4.8.3 requests/2.23.0 CPython/3.8.3 Linux/4.15.0-107-generic ubuntu/16.04.6 glibc/2.23
                UID:GID : 20331:2019
             netrc file : None
           offline mode : False
</div>

So since the bioconda channel is already available, we do not have to specify it when we create the environment and install the software:

	conda create -n htstream-1.3.2 htstream

Enter "y" when prompted. This command is creating an environment called "htstream-1.3.2" and installing the "htstream" package into it. By default, it installs the latest version. This might take a little while to run.

<div class="output">(base) joshi@rafter-5:/share/workshop/prereq_workshop/joshi/software$ conda create -n htstream-1.3.2 htstream
Collecting package metadata (current_repodata.json): done
Solving environment: done


==> WARNING: A newer version of conda exists. <==
  current version: 4.8.3
  latest version: 4.8.5

Please update conda by running

    $ conda update -n base -c defaults conda



## Package Plan ##

  environment location: /share/workshop/prereq_workshop/joshi/software/miniconda3/envs/htstream-1.3.2

  added / updated specs:
    - htstream


The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    _libgcc_mutex-0.1          |      conda_forge           3 KB  conda-forge
    _openmp_mutex-4.5          |            1_gnu          22 KB  conda-forge
    boost-1.70.0               |   py38h9de70de_1         355 KB  conda-forge
    boost-cpp-1.70.0           |       h7b93d67_3        21.1 MB  conda-forge
    bzip2-1.0.8                |       h516909a_3         398 KB  conda-forge
    ca-certificates-2020.6.20  |       hecda079_0         145 KB  conda-forge
    certifi-2020.6.20          |   py38h32f6830_0         151 KB  conda-forge
    htstream-1.3.2             |       h5ca1c16_0         1.7 MB  bioconda
    ld_impl_linux-64-2.35      |       h769bd43_9         617 KB  conda-forge
    libblas-3.8.0              |      17_openblas          11 KB  conda-forge
    libcblas-3.8.0             |      17_openblas          11 KB  conda-forge
    libffi-3.2.1               |    he1b5a44_1007          47 KB  conda-forge
    libgcc-ng-9.3.0            |      h5dbcf3e_17         7.8 MB  conda-forge
    libgfortran-ng-7.5.0       |      hae1eefd_17          22 KB  conda-forge
    libgfortran4-7.5.0         |      hae1eefd_17         1.3 MB  conda-forge
    libgomp-9.3.0              |      h5dbcf3e_17         378 KB  conda-forge
    liblapack-3.8.0            |      17_openblas          11 KB  conda-forge
    libopenblas-0.3.10         |pthreads_hb3c22a3_4         7.8 MB  conda-forge
    libstdcxx-ng-9.3.0         |      h2ae2ef3_17         4.0 MB  conda-forge
    lz4-c-1.9.2                |       he1b5a44_3         203 KB  conda-forge
    ncurses-6.2                |       he1b5a44_1         993 KB  conda-forge
    numpy-1.19.1               |   py38hbc27379_2         5.3 MB  conda-forge
    openssl-1.1.1h             |       h516909a_0         2.1 MB  conda-forge
    pip-20.2.3                 |             py_0         1.1 MB  conda-forge
    python-3.8.5               |h1103e12_9_cpython        21.9 MB  conda-forge
    python_abi-3.8             |           1_cp38           4 KB  conda-forge
    readline-8.0               |       he28a2e2_2         281 KB  conda-forge
    setuptools-49.6.0          |   py38h32f6830_1         940 KB  conda-forge
    sqlite-3.33.0              |       h4cf870e_0         1.4 MB  conda-forge
    tk-8.6.10                  |       hed695b0_0         3.2 MB  conda-forge
    wheel-0.35.1               |     pyh9f0ad1d_0          29 KB  conda-forge
    xz-5.2.5                   |       h516909a_1         343 KB  conda-forge
    zlib-1.2.11                |    h516909a_1009         106 KB  conda-forge
    ------------------------------------------------------------
                                           Total:        83.6 MB

The following NEW packages will be INSTALLED:

  _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge
  _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-1_gnu
  boost              conda-forge/linux-64::boost-1.70.0-py38h9de70de_1
  boost-cpp          conda-forge/linux-64::boost-cpp-1.70.0-h7b93d67_3
  bzip2              conda-forge/linux-64::bzip2-1.0.8-h516909a_3
  ca-certificates    conda-forge/linux-64::ca-certificates-2020.6.20-hecda079_0
  certifi            conda-forge/linux-64::certifi-2020.6.20-py38h32f6830_0
  htstream           bioconda/linux-64::htstream-1.3.2-h5ca1c16_0
  icu                conda-forge/linux-64::icu-67.1-he1b5a44_0
  ld_impl_linux-64   conda-forge/linux-64::ld_impl_linux-64-2.35-h769bd43_9
  libblas            conda-forge/linux-64::libblas-3.8.0-17_openblas
  libcblas           conda-forge/linux-64::libcblas-3.8.0-17_openblas
  libffi             conda-forge/linux-64::libffi-3.2.1-he1b5a44_1007
  libgcc-ng          conda-forge/linux-64::libgcc-ng-9.3.0-h5dbcf3e_17
  libgfortran-ng     conda-forge/linux-64::libgfortran-ng-7.5.0-hae1eefd_17
  libgfortran4       conda-forge/linux-64::libgfortran4-7.5.0-hae1eefd_17
  libgomp            conda-forge/linux-64::libgomp-9.3.0-h5dbcf3e_17
  liblapack          conda-forge/linux-64::liblapack-3.8.0-17_openblas
  libopenblas        conda-forge/linux-64::libopenblas-0.3.10-pthreads_hb3c22a3_4
  libstdcxx-ng       conda-forge/linux-64::libstdcxx-ng-9.3.0-h2ae2ef3_17
  lz4-c              conda-forge/linux-64::lz4-c-1.9.2-he1b5a44_3
  ncurses            conda-forge/linux-64::ncurses-6.2-he1b5a44_1
  numpy              conda-forge/linux-64::numpy-1.19.1-py38hbc27379_2
  openssl            conda-forge/linux-64::openssl-1.1.1h-h516909a_0
  pip                conda-forge/noarch::pip-20.2.3-py_0
  python             conda-forge/linux-64::python-3.8.5-h1103e12_9_cpython
  python_abi         conda-forge/linux-64::python_abi-3.8-1_cp38
  readline           conda-forge/linux-64::readline-8.0-he28a2e2_2
  setuptools         conda-forge/linux-64::setuptools-49.6.0-py38h32f6830_1
  sqlite             conda-forge/linux-64::sqlite-3.33.0-h4cf870e_0
  tk                 conda-forge/linux-64::tk-8.6.10-hed695b0_0
  wheel              conda-forge/noarch::wheel-0.35.1-pyh9f0ad1d_0
  xz                 conda-forge/linux-64::xz-5.2.5-h516909a_1
  zlib               conda-forge/linux-64::zlib-1.2.11-h516909a_1009
  zstd               conda-forge/linux-64::zstd-1.4.5-h6597ccf_2


Proceed ([y]/n)? y


Downloading and Extracting Packages
python_abi-3.8       | 4 KB      | ##################################### | 100% 
libopenblas-0.3.10   | 7.8 MB    | ##################################### | 100% 
ca-certificates-2020 | 145 KB    | ##################################### | 100% 
openssl-1.1.1h       | 2.1 MB    | ##################################### | 100% 
boost-cpp-1.70.0     | 21.1 MB   | ##################################### | 100% 
wheel-0.35.1         | 29 KB     | ##################################### | 100% 
readline-8.0         | 281 KB    | ##################################### | 100% 
libffi-3.2.1         | 47 KB     | ##################################### | 100% 
htstream-1.3.2       | 1.7 MB    | ##################################### | 100% 
libstdcxx-ng-9.3.0   | 4.0 MB    | ##################################### | 100% 
libgcc-ng-9.3.0      | 7.8 MB    | ##################################### | 100% 
ld_impl_linux-64-2.3 | 617 KB    | ##################################### | 100% 
libblas-3.8.0        | 11 KB     | ##################################### | 100% 
zlib-1.2.11          | 106 KB    | ##################################### | 100% 
xz-5.2.5             | 343 KB    | ##################################### | 100% 
python-3.8.5         | 21.9 MB   | ##################################### | 100% 
libgfortran4-7.5.0   | 1.3 MB    | ##################################### | 100% 
bzip2-1.0.8          | 398 KB    | ##################################### | 100% 
certifi-2020.6.20    | 151 KB    | ##################################### | 100% 
liblapack-3.8.0      | 11 KB     | ##################################### | 100% 
ncurses-6.2          | 993 KB    | ##################################### | 100% 
pip-20.2.3           | 1.1 MB    | ##################################### | 100% 
lz4-c-1.9.2          | 203 KB    | ##################################### | 100% 
sqlite-3.33.0        | 1.4 MB    | ##################################### | 100% 
tk-8.6.10            | 3.2 MB    | ##################################### | 100% 
setuptools-49.6.0    | 940 KB    | ##################################### | 100% 
libgfortran-ng-7.5.0 | 22 KB     | ##################################### | 100% 
numpy-1.19.1         | 5.3 MB    | ##################################### | 100% 
libcblas-3.8.0       | 11 KB     | ##################################### | 100% 
_libgcc_mutex-0.1    | 3 KB      | ##################################### | 100% 
boost-1.70.0         | 355 KB    | ##################################### | 100% 
libgomp-9.3.0        | 378 KB    | ##################################### | 100% 
_openmp_mutex-4.5    | 22 KB     | ##################################### | 100% 
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
#
# To activate this environment, use
#
#     $ conda activate htstream-1.3.2
#
# To deactivate an active environment, use
#
#     $ conda deactivate
</div>

Once that finishes, you need to activate the environment:

	conda activate htstream-1.3.2

Now, you should have access to the tools in the htstream package:

	hts_Stats -h
	hts_Overlapper -h

Lets see what our environment looks like after running this.

	conda list

<div class="output">(htstream-1.3.2) joshi@rafter-5:/share/workshop/prereq_workshop/joshi/software$ conda list
# packages in environment at /share/workshop/prereq_workshop/joshi/software/miniconda3/envs/htstream-1.3.2:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       1_gnu    conda-forge
boost                     1.70.0           py38h9de70de_1    conda-forge
boost-cpp                 1.70.0               h7b93d67_3    conda-forge
bzip2                     1.0.8                h516909a_3    conda-forge
ca-certificates           2020.6.20            hecda079_0    conda-forge
certifi                   2020.6.20        py38h32f6830_0    conda-forge
htstream                  1.3.2                h5ca1c16_0    bioconda
icu                       67.1                 he1b5a44_0    conda-forge
ld_impl_linux-64          2.35                 h769bd43_9    conda-forge
libblas                   3.8.0               17_openblas    conda-forge
libcblas                  3.8.0               17_openblas    conda-forge
libffi                    3.2.1             he1b5a44_1007    conda-forge
libgcc-ng                 9.3.0               h5dbcf3e_17    conda-forge
libgfortran-ng            7.5.0               hae1eefd_17    conda-forge
libgfortran4              7.5.0               hae1eefd_17    conda-forge
libgomp                   9.3.0               h5dbcf3e_17    conda-forge
liblapack                 3.8.0               17_openblas    conda-forge
libopenblas               0.3.10          pthreads_hb3c22a3_4    conda-forge
libstdcxx-ng              9.3.0               h2ae2ef3_17    conda-forge
lz4-c                     1.9.2                he1b5a44_3    conda-forge
ncurses                   6.2                  he1b5a44_1    conda-forge
numpy                     1.19.1           py38hbc27379_2    conda-forge
openssl                   1.1.1h               h516909a_0    conda-forge
pip                       20.2.3                     py_0    conda-forge
python                    3.8.5           h1103e12_9_cpython    conda-forge
python_abi                3.8                      1_cp38    conda-forge
readline                  8.0                  he28a2e2_2    conda-forge
setuptools                49.6.0           py38h32f6830_1    conda-forge
sqlite                    3.33.0               h4cf870e_0    conda-forge
tk                        8.6.10               hed695b0_0    conda-forge
wheel                     0.35.1             pyh9f0ad1d_0    conda-forge
xz                        5.2.5                h516909a_1    conda-forge
zlib                      1.2.11            h516909a_1009    conda-forge
zstd                      1.4.5                h6597ccf_2    conda-forge
</div>


To deactivate your environment:

	conda deactivate

And then to deactivate conda, run it again:

	conda deactivate

Compared to typical installations such as the one shown above, conda installs only requires one command. Conda is quicker, more user friendly, and more commonly results in success. 
    
---
# A few things for your future Conda usage:

<object data="https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>

### When running conda on your own computer you will need to add channels, which is where conda will look when performing package installs

	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
