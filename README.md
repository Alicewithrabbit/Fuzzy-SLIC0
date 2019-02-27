# C implementation of Fuzzy SLIC with Matlab Interface
## Auto compactness coefficient selection version
Copyright (c) 2019, Chong WU 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
  
* Neither the name of City University of Hong Kong nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Usage: 

You need Matlab and a C compiler to compile this code before using it.

Input variables are:

1. 8 bit images (color or grayscale)
2. Number of intended superpixels (optional, default is 200)

Ouputs are:

1. Labels (in raster scan order)
2. Number of labels (final output superpixels) in the image

Syntax: 

[labels, numlabels] = fuzzyslic0mex(img,200);

NOTES:

Number of returned superpixels may be different from the input number of superpixels.

## Demo illustration

![Comparison result](https://github.com/Alicewithrabbit/Fuzzy-SLIC0/blob/master/demo.jpg)

<center><font color=grey>**Figure 1. Comparison of Fuzzy SLIC0 and Fuzzy SLIC. Left one is the result obtained by Fuzzy SLIC0 and right one is the result obtained by Fuzzy SLIC**</font></center>
