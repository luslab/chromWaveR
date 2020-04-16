# chromWaveR: R package to predict nucleosome occupancy from DNA sequence using pre-trained chromWave models.

## Abstract 


## Requirements 

* Install following python packages:
   
    ```console 
    # CPU version 
    pip install tensorflow==2 keras numpy
    # GPU version 
    pip install tensorflow-gpu==2 keras numpy
    ``` 

* ALTERNATIVE: Setup Tensorflow/Keras CONDA env: 

    1) Download and install [CONDA](https://www.anaconda.com/distribution/).

    2) Install python packages
    
    ```console 
    conda create -n keras python
    # CPU version 
    pip install tensorflow==2 keras numpy
    # GPU version 
    pip install tensorflow-gpu==2 keras numpy
    ``` 

* Install following R packages from CRAN: 
    ```R 
    install.packages('tidyverse', 'reticulate', 'kerasR')
    ``` 

* Install chromWaveR:
    ```R 
    install.packages('devtools')
    devtools::install_github('luslab/chromWaveR')
    ``` 

## Tutorial  

For an example of how to use the package see [here](https://github.com/luslab/chromWaveR/blob/master/inst/example/chromWaveR_example.md). 



## Reference 
**ChromWave: Deciphering the DNA-encoded competition between DNA-binding factors with deep neural networks**

Preprint: Add link upon submission.
