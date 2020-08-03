# EstimatingHomography

Input Files: 1. Image set img1.jpg img2.jpg
Output: 1: image: extracted corners image
           image: matched corners image
           image: warped and aligned image
           
# Feature extraction and matching:
For feature extraction, use the Harris corner detector that was already developed in a previous assignment. Similarly, employ normalized cross correlation (NCC) for matching corner points between pair of images. Figure 2 shows extracted Harris corner points, and matched correspondences by NCC.  

# Homography estimation using RANSAC framework
Here, a homography is used for stitching images. The homography can be computed using matched correspondences from the previous step. However, initially matched correspondences contain lots of outliers, and the incorrect homography could be obtained. Thus, implement a RANSAC framework for estimating the correct homography by rejecting outliers. 

# Image warping
With the computed homography from the previous step, draw the aligned images by inverse warping with a bilinear interpolation.
