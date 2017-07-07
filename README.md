This is a MATLAB program which convolves radio-telescope data onto the grid, one of the most time-consuming processing steps to create a sky image. 

## W-projection explanation

The first step of the process is pair of antennas sampling the electromagnetic spectrum at a high rate. The integrated product of such a sample is called a visibility, and has an associated (u,v,w) coordinate, depending on the position of the antennas, the position of the observed source, the frequency of te observed signal and the time. The (u,v,w) coordinates for an antenna pair change over time due to rotation of the earth.

Since the visibilities are sampled in the Fourier domain, the next step is to perform an inverse Fourier transform operation to create a sky image. However, with a lot of data this process consumes an exorbitant amount of time. Therefore these visibilities are first placed on a UV-grid. A two-dimensional FFT can then be used, converting the UV-grid to a sky image in a much more effective fashion. 

A visibility is not simply placed on the UV-grid, since it contributes to neighboring grid points as well. This contribution is obtained by convolving the visibility with a convolution matrix, i.e. by multiplying the visibility with points in a so-called support matrix and adding this result to the grid. In this application the convolution matrix is considered as a precomputed matrix of real weights.

This would be all that is necessary when the earth and the observed part of the sky are considered flat. Using a third coordinate w we can account for this extra dimension by using a technique called wide-field imaging. We implement this extra dimension by making different convolution matrices for different values of w,  typically using several tens of w-planes.

This code basically does the gridding step: it puts every visibility or (u,v,w) point on a grid point and makes it contribute to neighboring grid points by using a convolution matrix dependent on w. It is also able to do the inverse process: converting a grid back to visibilities.

Code explanation


Choosing the oversampling mode leads to the support kernel having a much 
