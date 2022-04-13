# Step-by-step Omnipose installation instructions (command line based, assumes cellpose is not already installed)

### Windows

1. Install miniconda (https://docs.conda.io/en/latest/miniconda.html)

Warning: If added to PATH, can run miniconda directly from Command Window, but could cause interferences. If not added to PATH, need to use the separate Anaconda Prompt to run.

2. Install git (https://git-scm.com/download/win)

3. In Command Window, navigate to desired folder and install Omnipose (ie C:\Users\Name\omnipose\).
```
git clone https://github.com/kevinjohncutler/cellpose.git
```

4. In Command Window/Anaconda Prompt, navigate to folder with environment.yml file (probably `cd cellpose`), activate environment, and install Omnipose.
```
conda env create -f environment.yml
conda activate cellpose
pip install git+https://github.com/kevinjohncutler/cellpose.git
```

May need: 
```
pip install sklearn
pip install -e C:\Users\Name\omnipose\cellpose
```


### Linux (Debian/Ubuntu)

1. Install Anaconda (https://www.anaconda.com/products/distribution)
```
sudo apt-get update
sudo apt-get install curl
sudo apt-get install wget
sudo apt install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
```

Tested and known to work: Anaconda3-2021.05-Linux-x86_64
```
(sudo) wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh 
bash Anaconda3-2021.05-Linux-x86_64.sh 
rm -rf Anaconda3-2021.05-Linux-x86_64.sh 
```

Latest Anaconda: (as of April 2022; find the full availability here https://repo.anaconda.com/archive/)
```
(sudo) wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh 
rm -rf Anaconda3-2021.05-Linux-x86_64.sh 
```

Latest miniconda:
```
(sudo) wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
rm -rf bash Miniconda3-latest-Linux-x86_64.sh
```

"Do you accept the license terms?" >yes

2. Restart shell:
```
source ~/.bashrc
```

Optional: to prevent base environment from loading automatically, run:
```
conda config --set auto_activate_base false
```
and restart shell.

3.  Install git (https://git-scm.com/download/linux): 
```
sudo apt-get install git
```

4. Clone Omnipose/Cellpose from Github
Navigate to folder you want to install Omnipose in (using cd), then run:
```
git clone https://github.com/kevinjohncutler/cellpose.git
```

Navigate to folder containing environment.yml file (ie 'cellpose') and create cellpose environment
```
cd cellpose
conda env create -f environment.yml
```

5. Activate Omnipose/Cellpose environment (repeat step each time you want to use Omnipose)
```
conda activate cellpose
```

6. Install cellpose 
```
pip install git+https://github.com/kevinjohncutler/cellpose.git
```

7. For updating: 
```
git -C ~/cellpose pull origin master
```




### MacOS (not yet tested, assumes homebrew is installed)

1. Install Anaconda (https://www.anaconda.com/products/distribution)
```
brew install curl
brew install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
brew install wget
```

Latest Anaconda: (as of April 2022; find the full availability here https://repo.anaconda.com/archive/)
```
(sudo) wget https://repo.anaconda.com/archive/Anaconda3-2021.11-MacOSX-x86_64.sh
bash Anaconda3-2021.11-MacOSX-x86_64.sh
```

Latest miniconda:
```
(sudo) wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

2. Restart shell:
```
source ~/.bashrc
```

Optional: to prevent base environment from loading automatically, run:
```
conda config --set auto_activate_base false
```
and restart shell.

3.  Install git (https://git-scm.com/download/mac): 
```
brew install git
```

4. Clone Omnipose/Cellpose from Github
Navigate to folder you want to install Omnipose in (using cd), then run:
```
git clone https://github.com/kevinjohncutler/cellpose
```

Navigate to folder containing environment.yml file (ie 'cellpose') and create cellpose environment
```
cd cellpose
conda env create -f environment.yml
```

5. Activate Omnipose/Cellpose environment (repeat step each time you want to use Omnipose)
```
conda activate cellpose
```

6. Install cellpose 
```
pip install git+https://github.com/kevinjohncutler/cellpose.git
```

7. For updating: 
```
git -C ~/cellpose pull origin master
```










