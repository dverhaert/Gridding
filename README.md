This is a MATLAB program which convolves radio-telescope data onto the grid, one of the most time-consuming processing steps to create a sky image. 

## W-projection explanation

The first step of the process of creating a sky image out of radio-telescope data is pair of antennas sampling the electromagnetic spectrum at a high rate. The integrated product of such a sample is called a visibility, and has an associated (*u,v,w*) coordinate, depending on the position of the antennas, the position of the observed source, the frequency of te observed signal and the time. The (*u,v,w*) coordinates for an antenna pair change over time due to rotation of the earth.

Since the visibilities are sampled in the Fourier domain, the next step is take these points and perform an inverse Fourier transform operation operation on all of them to create a sky image. However, with a lot of data this process consumes an exorbitant amount of time, and is not feasible for practical use. Therefore these visibilities are first placed on a UV-grid, so that a two-dimensional FFT can be used. Doing this converts the UV-grid to a sky image in a much more effective fashion. 

While a visibility can be simply rounded to the closest grid point, this is not fully correct: a visibility contributes to neighboring grid points as well. Therefore the contribution of a visibility to the grid is obtained by convolving it with a convolution matrix. Since a convolution represents a multiplication in the Fourier domain, the contribution can easily be obtained by multiplying the visibility with points in the convolution matrix and adding this result to the grid. In this application the convolution matrix is considered as a precomputed matrix of real weights, and is referred to as the support matrix.

This would be all that is necessary when the earth and the observed part of the sky are considered flat. However, since this is not correct, this dimension must be accounted for. This is done by using a third coordinate w and a technique called wide-field imaging. We implement this extra dimension by making different convolution matrices for different values of w,  typically using several tens of w-planes. When all the visibilities have been convolved to the grid, a final IFFT is performed to obtain the desired sky image.

This code basically does the gridding step: it puts every visibility or (*u,v,w*) point on a grid point and makes it contribute to neighboring grid points by using a convolution matrix dependent on w. It is also able to do the inverse process: converting a grid back to visibilities.

## Code explanation

A implementation of gridding radio-telescope data using the W-projection algorithm has been created in MATLAB. Although some code projects already exist which can do about the same thing, they are rather complicated and inaccessible. Two of such implementations where this program is based on are [1] and [2], both created in C++. 

### Initialization
At the top of the code the user is able to choose a *mode*. This *mode* variable can be set to either 1, 2 or 3, representing the simple, w-projection and interpolation mode respectively. The functioning of these modes is explained later.

After that, many variables can be customized depending on the user's requirements. The first value to be set is *oversampling*. Oversampling is the process of making the support kernel more dense in comparison to the grid it's in. This way, instead of (u,v) coordinates simply being rounded to the nearest integer, they get scaled according to their fraction depending on how many times we're oversampling. Oversampling with a factor 2 makes (u,v) coordinates round to halves, with a factor 3 to thirds etc. This way, the accuracy is improved since we look at fractions more carefully but the computation time will be longer.

Next, some variables can be set to tweak the program as desired. These determine the size of the support kernel, the size of the total grid, the number of frequency channels, and the number of points taken from the dataset. Real (*u,v,w*) coordinates from a six-hour LOFAR observation with 44 antennas (946 baselines, 10s integration time and one subband of 16 frequency channels) were used as taken from [1]. These are further subdivided in timesteps and blocks, but can be tweaked as necessary. Also a scaling factor is applied so taht the points are adjusted to fit in a ([-1024, 1024, -1024, 1024]) grid.

Frequencies are set, differing slightly for every channel. The frequencies in this case used are around 60 MHz, but of course any values can be used as long as they are in the radio frequency. After this, the values of the visibilities are set. In this example they are all set to one, but of course any (complex) value can be taken, or a dataset of actual visibilities can be used.

### Mode 1: Simple mode
#### Generating the kernel
The simple mode uses a rather standard kernel which has two dimensions, and so matches with the grid. It is generated by the following formula:

*support(i,j) = exp(-((i-ccenter)^2 + (j-ccenter)^2));* (1)

Where ccenter is the center of the support kernel. With this formula, the center of the support kernel will be equal to one, and the the farther away one is from the center the lower the values are. For example, with a 9x9 support kernel, the values at the corners and so the farthest away from the center are 1,2664e-14. 

Formula (1) holds for the general case where we're not *oversampling*. In the case that we're oversampling, we add a scaling factor to it so that the values at the edges remain the same.

#### Gridding
Since no w-coordinates are used in the simple mode, just the first and second column of the data are taken. These *u* and *v* coordinates are scaled according to the defined frequency, cellsize and scaling factor. 

Next the program enters its main execution loop, going over all data points. The current scaled _u_ and _v_ points are set to the variables _uscaled_ and _vscaled_ respectively. The fractional part of the current _u_ and _v_ coordinate is determined by first rounding and then computing the difference.

The index we will update in the grid is equal to the integer value of a point plus half of the grid size (so that we look at the point from the middle of the grid). Now, we convolve the point with the support kernel by multiplying the visibility with some values of the support matrix. Without oversampling we do this with support values taken from around the center, but with oversampling this changes depending on the fraction. The results are added to the grid.

Finally, after this process is repeated for every grid point, the resulting grid is flipped so that _v_ is not inverted and the result is plotted.

### Mode 2: W-projection mode
#### Generating the kernel
The support kernel for the w-projection mode has three dimensions, it can be seen as a two-dimensional kernel for different values of _w_. The formula to generate the support changes to the following one:

*support(i,j,k) = exp(-((i-ccenter)^2 + (j-ccenter)^2)/(k\*fScale));* (2)

*fScale = sqrt(abs(w)\*wcellsize\*frequencies(1))/cellsize;* (3)

This extra scaling factor added makes it so that the kernel will decay normally with low values of _w_, but will will decay at a much slower rate with higher values of _w_. This means that, for example in the highest _w_-plane in the 9x9 example, the value at the edge will be 0.8862 (instead of 0.0004 in the lowest w-plane!). 

#### Gridding
Generally, the gridding procedure is the same as the one for the simple mode. This time all data columns are taken and the *u*, *v* and *w* coordinates are scaled according to the defined frequency and cellsize and their fractional part is determined. In the gridding procedure again the grid is determined by convolving the two-dimensional (*u,v*) support kernel with every visibility, but this time instead of using one fixed kernel it depends on the value of *w*. Again, the results are added to the grid and the result is plotted.

### Mode 3: Interpolation mode
#### Generating the kernel
The third implemented mode is the interpolation mode. The key difference of this mode is explained in the gridding part of this section. The kernel is generated in the same way as in mode 2 but without oversampling, so it's a three-dimensional matrix of which the size will be much smaller. 

#### Gridding

The gridding procedure for the interpolation mode starts in the same way as the w-projection: all data columns are taken and the *u*, *v* and *w* coordinates are scaled according to the defined frequency and cellsize and their fractional part is determined. 

The key difference is that instead of expanding the support kernel and moving over this kernel depending on the fractional part of *u* and *v* as we did in the previous w-projection mode, one accounts for a point in proportion to it's fraction. For example, if a *u*-coordinate is equal to 119.9, we perform the the gridding operation with the smaller (non-oversampled) kernel for *u*-coordinates of both 119 and 120. However, the operation for the point 119 will only be counted 0.1 times, while the operation for the point 120 will be counted 0.9 times. This way, a point is weighted much more when the fraction is near it. Recall that with the use of oversampling in the previous modes the support kernel is made bigger and we go over points depending on the fraction, while here in the interpolation mode we keep the original support kernel size but perform the computation for multiple points: this is the key difference between these modes. Since we can have fractions in the *u*, *v* and *w* directions, this implies a total of 2^3 = 8 computations for every grid point (instead of one with an oversampled support kernel).

Basically, in return for having to make much more computations per grid point we have a much smaller kernel. When oversampling with a small factor this is quite a bit faster (although less precise) than interpolating. However with large oversampling factors interpolating becomes more efficient.  

### Degridding
Degridding is basically the gridding process, but ran in reverse. For every point, we look at the points around that point according to the size of the support kernel and add the sum of these values. Finally, we divide by the sum of the support kernel. Note that this leads to very different visibilities in comparison to those input into the gridding function:  To obtain correct visibilities, a more advanced algorithm has to be used.


[1] J. W. Romein. An Efficient Work-Distribution Strategy for Gridding Radio-Telescope Data on GPUs. *ACM, ICS'12*, 2012.

[2] T.J. Cornwell. The impact of convolutional resampling on ASKAP and SKA Phase 1 computing costs. *CONRAD-SW-0008*, 2007.
