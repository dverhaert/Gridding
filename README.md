This is a MATLAB program which convolves radio-telescope data onto the grid, one of the most time-consuming processing steps to create a sky image. 

### W-projection explanation

The first step of the process of creating a sky image out of radio-telescope data is pair of antennas sampling the electromagnetic spectrum at a high rate. The integrated product of such a sample is called a visibility, and has an associated (u,v,w) coordinate, depending on the position of the antennas, the position of the observed source, the frequency of te observed signal and the time. The (u,v,w) coordinates for an antenna pair change over time due to rotation of the earth.

Since the visibilities are sampled in the Fourier domain, the next step is take these points and perform an inverse Fourier transform operation operation on all of them to create a sky image. However, with a lot of data this process consumes an exorbitant amount of time, and is not feasible for practical used. Therefore these visibilities are first placed on a UV-grid, so that a two-dimensional FFT can be used. Doing this converts the UV-grid to a sky image in a much more effective fashion. 

While a visibility can be simply rounded to the closest grid point, this is not fully correct: a visibility contributes to neighboring grid points as well. Therefore this contribution is obtained by convolving this visibility with a convolution matrix. Since a convolution represents a multiplication in the Fourier domain, the contribution can easily be obtained by multiplying the visibility with points in the convolution matrix and adding this result to the grid. In this application the convolution matrix is considered as a precomputed matrix of real weights, and is referred to as the support matrix.

This would be all that is necessary when the earth and the observed part of the sky are considered flat. However, since this is not correct, this dimension must be accounted for. This is done by using a third coordinate w and a technique called wide-field imaging. We implement this extra dimension by making different convolution matrices for different values of w,  typically using several tens of w-planes. When all the visibilities have been convolved to the grid, a final IFFT is performed to obtain the desired sky image.

This code basically does the gridding step: it puts every visibility or (u,v,w) point on a grid point and makes it contribute to neighboring grid points by using a convolution matrix dependent on w. It is also able to do the inverse process: converting a grid back to visibilities.

### Code explanation

A implementation of gridding radio-telescope data using the W-projection algorithm has been created in MATLAB. Although some code projects already exist which can do about the same thing, they are rather complicated and inaccessible. Two of such implementations where this program is based on are [1] and [2], both created in C++. 

Choosing the oversampling mode leads to the support kernel having a much 


Generally, oversampling is slightly faster than interpolating. Testing some data we noticed that generally the oversampling mode was about 1,5-2x faster. This can be explained by the extra computations done in the most computationally intensive innter loop of the code.


[1] J. W. Romein. An Efficient Work-Distribution Strategy for Gridding Radio-Telescope Data on GPUs. *ACM, ICS'12*, 2012.
[2] T.J. Cornwell. The impact of convolutional resampling on ASKAP and SKA Phase 1 computing costs. *CONRAD-SW-0008*, 2007.
