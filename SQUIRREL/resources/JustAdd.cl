//#pragma OPENCL EXTENSION cl_khr_fp64: enable

_kernel void justAdd(
    _global float* pixels1,
    _global float* pixels2,
    _global float* pixelsSum,
    const int width,
    const int height
){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int p = y*width+x;

    float v1 = pixels1[p];
    float v2 = pixels2[p];

    float vSum = v1 + v2;

    pixelsSum[p] = vSum;

}