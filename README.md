# VanishingPointDetection

2-Line Exhaustive Searching for Real-Time Vanishing Point Estimation in Manhattan World,Xiaohu Lu, JianYao, Haoang Li, Yahui Liu and Xiaofeng Zhang, WACV2017.

http://xiaohulugo.github.io/papers/Vanishing_Point_Detection_WACV2017.pdf

Prerequisites:
---
OpenCV > 2.4.x

Usage:
---
1. build the project with Cmake
2. set the internal parameters of the input image (it's ok to use approximate values, for example: pp(cols/2,rows/2), f=max(cols,rows), but the result will be a little worse)

Performance:
---
40ms on a computer with Intel Core i5-3550p CPU without any optimization and parallel computation in the Release mode.

Please cite this paper if you use this data or code:

    @InProceedings{Lu_2017_WACV,
    author = {Lu, Xiaohu and Yao, Jian and Li, Haoang and Liu, Yahui and Zhang, Xiaofeng},
    title = {2-Line Exhaustive Searching for Real-Time Vanishing Point Estimation in Manhattan World},
    booktitle = {IEEE Winter Conference on Applications of Computer Vision (WACV)},
    month = {March},
    year = {2017}
    }
    
Feel free to correct my code, if you spotted the mistakes. You are also welcomed to Email me: fangzelu@gmail.com
