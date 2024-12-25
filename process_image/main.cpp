#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include "stb_image.h"
#include "stb_image_write.h"
#include <cmath>
#define MAX_FILE_NAME_SIZE 256


void ImageProcessingGray(unsigned char* pOut,
                         unsigned char* pIn,
                         size_t nWidth,
                         size_t nHeight);

void Convolution2D(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, const float* pKernel, int iKernelSize);

void ResizeBilinear(const unsigned char* pbIn, int lWidthIn, int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut);

void ResizeNearestNeighbor(const unsigned char* pbIn, int lWidthIn, int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut);

void MedianFilter2D(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, int iKernelSize);

void AddNoise(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight);

void ThresholdSegmentation(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, unsigned char threshold);

void ConnectedComponentLabeling(const unsigned char* pBin, unsigned char* pLabel, int iWidth, int iHeight, int connectivity);
                           
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Error: the program must have at least one input parameter - the path to the image file!"
                  << std::endl;
        return 1;
    }

    // get input image file path:
    // null character taken into account:
    if (std::strlen(argv[1]) >= MAX_FILE_NAME_SIZE) {
        std::cerr << "Error: the input file path must be less than "
                  << MAX_FILE_NAME_SIZE << " characters!" << std::endl;
        return 2;
    }
    char input_file_name[MAX_FILE_NAME_SIZE];
    strcpy(input_file_name,
           argv[1]);

    // create output image file path with adding suffix to input file path:
    const char append_suffix[] = "_out.png";
    char output_file_name[MAX_FILE_NAME_SIZE+sizeof(append_suffix)];
    strcpy(output_file_name,
           input_file_name);
    strcat(output_file_name,
           append_suffix);

    // read input image file:
    int input_width, input_height, input_bytes_per_pixel;
    unsigned char *input_image_data;
    input_image_data = stbi_load(input_file_name,
                                 &input_width,
                                 &input_height,
                                 &input_bytes_per_pixel,
                                 0);
    if (input_image_data == NULL) {
        std::cerr << "Error: can not read image from file by path: "
                  << input_file_name << std::endl;
        return 3;
    }
    std::cout << "Readed image with width = "
              << input_width << " pixels, height = "
              << input_height << " pixels, and bit depth = "
              << input_bytes_per_pixel * 8
              << " bits per pixel" << std::endl;

    // initialize data for output:
    int output_width = input_width,
        output_height = input_height,
        output_bytes_per_pixel = input_bytes_per_pixel;
    size_t output_bytes_num = output_width * output_height * output_bytes_per_pixel;
    unsigned char *output_image_data = new unsigned char[output_bytes_num];

    // process input image, result of processing in output image data

    // ImageProcessingGray(output_image_data,
    //                     input_image_data,
    //                     output_width,
    //                     output_height);



            
    // const int kernel_size = 3;
    // float edgeKernel[] = {
    // -1, -1, -1,
    // 2,  2, 2,
    // -1, -1, -1
    // };

    
    // const int kernel_size = 5;
    // float edgeKernel[] = {
    //     1, 4, 7, 4, 1,
    //     4, 16, 26, 16, 4,
    //     7, 26, 41, 26, 7,
    //     4, 16, 26, 16, 4,
    //     1, 4, 7, 4, 1
    // };
    // for (int i = 0; i < kernel_size * kernel_size; ++i) {
    //     edgeKernel[i] /= 273;
    // }

    // const int kernel_size = 3;
    // float edgeKernel[] = {
    //     1, 1, 1,
    //     1, 1, 1,
    //     1, 1, 1
    // };

    // for (int i = 0; i < kernel_size * kernel_size; ++i) {
    //     edgeKernel[i] /= 9;
    // }

    

    // Convolution2D(input_image_data, output_image_data, input_width, input_height, edgeKernel, kernel_size);



    // INTERPOLATION
    // int output_width = input_width * 5,
    //     output_height = input_height * 5,
    //     output_bytes_per_pixel = input_bytes_per_pixel;
    // size_t output_bytes_num = output_width * output_height * output_bytes_per_pixel;
    // unsigned char *output_image_data = new unsigned char[output_bytes_num];
    // ResizeBilinear(input_image_data, input_width, input_height, output_image_data, output_width, output_height);
    // ResizeNearestNeighbor(input_image_data, input_width, input_height, output_image_data, output_width, output_height);

    // MEDIAN FILTERING
    // int output_width = input_width,
    //     output_height = input_height,
    //     output_bytes_per_pixel = input_bytes_per_pixel;
    // size_t output_bytes_num = output_width * output_height * output_bytes_per_pixel;
    // unsigned char *output_image_data = new unsigned char[output_bytes_num];
    // int aperture = 5;
    // AddNoise(input_image_data, output_image_data, input_width, input_height);
    // MedianFilter2D(input_image_data, output_image_data, input_width, input_height, aperture);

    // SEGMENTATION

    ThresholdSegmentation(input_image_data, output_image_data, input_width, input_height, 86);

    unsigned char *labeled_data = new unsigned char[output_bytes_num];
    for (int i = 0; i < output_bytes_num; ++i) labeled_data[i] = 0;
    
    ConnectedComponentLabeling(output_image_data, labeled_data, input_width, input_height, 4);
    output_image_data = labeled_data;

    // write output image to PNG file:
    int write_res2 = stbi_write_png(output_file_name,
                                   output_width,
                                   output_height,
                                   output_bytes_per_pixel,
                                   output_image_data,
                                   0);
    if (!write_res2) {
        std::cerr << "Error: can not write image to file by path: "
                  << output_file_name << std::endl;

        // free data (better to use smart pointers to prevent this code repeat):
        delete[] output_image_data;
        stbi_image_free(input_image_data);

        return 4;
    }
    std::cout << "Resulting output image writed to file with path: "
              << output_file_name << std::endl;

    // free data:
    delete[] output_image_data;
    delete[] labeled_data;
    stbi_image_free(input_image_data);

    return 0;
}

// use this function for write image processing
void ImageProcessingGray(unsigned char* pOut,
                        unsigned char* pIn,
                        size_t nWidth,
                        size_t nHeight)
{
    // ASSERT(pOut != NULL && pIn != NULL && nWidth > 0 && nHeight > 0)

    std::vector<int> histogram(256, 0);

    unsigned char* pTemp = pIn;
    for (size_t y = 0; y < nHeight; ++y) {
        for (size_t x = 0; x < nWidth; ++x) {
            histogram[*pTemp]++;
            pTemp++;
        }
    }

    std::ofstream csvFile("histogram.csv");
    if (csvFile.is_open()) {
        csvFile << "Intensity,Count\n";
        for (size_t i = 0; i < histogram.size(); ++i) {
            csvFile << i << "," << histogram[i] << "\n";
        }
        csvFile.close();
    }

    return;
}

void Convolution2D(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, const float* pKernel, int iKernelSize) {
    int iHalf = iKernelSize / 2;
    memcpy(pRes, pSrc, sizeof(*pRes) * iWidth * iHeight);
    for (int y = iHalf; y < iHeight - iHalf; ++y) {
        for(int x = iHalf; x < iWidth - iHalf; ++x){
            const float* pk = pKernel;
            const unsigned char* ps = &pSrc[(y - iHalf) * iWidth + x - iHalf];
            float iSum = 0;
            for (int v = 0; v < iKernelSize; ++v){
                for (int u = 0; u < iKernelSize; ++u) iSum += ps[u] * pk[u];
                pk += iKernelSize;
                ps += iWidth;
                }
            if (iSum > 255) iSum = 255;
            else if (iSum < 0) iSum = 0;
            pRes[y * iWidth + x] = (unsigned char)iSum;
        }
    }
}

void ResizeBilinear(const unsigned char* pbIn, int lWidthIn, int lHeightIn,
                    unsigned char* pbOut, int lWidthOut, int lHeightOut)
{
    for (int i = 0; i < lHeightOut; ++i)
    {
        double yy = (double)i * (double)lHeightIn / lHeightOut;
        int y = (int)yy; // целая часть yy
        double u = yy - (double)y; // дробная часть yy

        for (int j = 0; j < lWidthOut; ++j)
        {
            double xx = (double)j * (double)lWidthIn / lWidthOut;
            int x = (int)xx; // целая часть xx
            double v = xx - (double)x; // дробная часть xx

            // значения соседних пикселей
            int lP00 = pbIn[y * lWidthIn + x],
                lP01 = pbIn[y * lWidthIn + x + 1],
                lP10 = pbIn[(y + 1) * lWidthIn + x],
                lP11 = pbIn[(y + 1) * lWidthIn + x + 1];

            // билинейная интерполяция
            pbOut[i * lWidthOut + j] = (unsigned char)(
                (1. - u) * (1. - v) * lP00 +
                u * (1. - v) * lP10 +
                v * (1. - u) * lP01 +
                u * v * lP11
            );
        }
    }
}

void ResizeNearestNeighbor(const unsigned char* pbIn, int lWidthIn, int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut)
{
    for (int i = 0; i < lHeightOut; ++i)
    {
        for (int j = 0; j < lWidthOut; ++j)
        {
            double x_in = (double)j * (double)lWidthIn / lWidthOut;
            double y_in = (double)i * (double)lHeightIn / lHeightOut;

            int x = (int)round(x_in);
            int y = (int)round(y_in);

            x = std::min(std::max(x, 0), lWidthIn - 1);
            y = std::min(std::max(y, 0), lHeightIn - 1);

            pbOut[i * lWidthOut + j] = pbIn[y * lWidthIn + x];
        }
    }
}

void AddNoise(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight) {
    for (int i = 0; i < iWidth * iHeight; ++i) {
        double proba = (double)rand() / RAND_MAX;
        if (proba > 0.5) {
            double p = (double)rand() / RAND_MAX;
            if (p > 0.7) pRes[i] = 255;
            else pRes[i] = 0;
        } else {
            pRes[i] = pSrc[i];
        }
    }
}

void MedianFilter2D(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, int iKernelSize) {
    int iHalf = iKernelSize / 2;
    std::vector<unsigned char> window(iKernelSize * iKernelSize);

    memcpy(pRes, pSrc, sizeof(*pRes) * iWidth * iHeight);

    for (int y = iHalf; y < iHeight - iHalf; ++y) {
        for (int x = iHalf; x < iWidth - iHalf; ++x) {
            int index = 0;

            for (int v = -iHalf; v <= iHalf; ++v) {
                for (int u = -iHalf; u <= iHalf; ++u) {
                    window[index++] = pSrc[(y + v) * iWidth + (x + u)];
                }
            }

            std::sort(window.begin(), window.end());
            if (window.size() % 2 == 1) {
                pRes[y * iWidth + x] = window[window.size() / 2];
            } else {
                pRes[y * iWidth + x] = (window[window.size() / 2 - 1] + window[window.size() / 2]) / 2.0;
            }
        }
    }
}

void ThresholdSegmentation(const unsigned char* pSrc, unsigned char* pRes, int iWidth, int iHeight, unsigned char threshold) {
    for (int i = 0; i < iWidth * iHeight; ++i) {
        pRes[i] = (pSrc[i] >= threshold) ? 255 : 0;
    }
}

void ConnectedComponentLabeling(const unsigned char* pBin, unsigned char* pLabel, int iWidth, int iHeight, int connectivity) {
    memset(pLabel, 0, iWidth * iHeight);

    int label = 1;
    std::vector<std::pair<int, int>> directions;

    if (connectivity == 4) {
        directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    } else if (connectivity == 8) {
        directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {1, 1}, {-1, 1}, {1, -1}};
    }

    for (int y = 0; y < iHeight; ++y) {
        for (int x = 0; x < iWidth; ++x) {
            if (pBin[y * iWidth + x] == 255 && pLabel[y * iWidth + x] == 0) {

                std::queue<std::pair<int, int>> q;
                q.push({x, y});
                pLabel[y * iWidth + x] = label;

                while (!q.empty()) {
                    std::pair<int, int> front = q.front();
                    int cx = front.first;
                    int cy = front.second;
                    q.pop();


                    for (const auto& dir : directions) {
                        int dx = dir.first;
                        int dy = dir.second;

                        int nx = cx + dx;
                        int ny = cy + dy;

                        if (nx >= 0 && nx < iWidth && ny >= 0 && ny < iHeight &&
                            pBin[ny * iWidth + nx] == 255 && pLabel[ny * iWidth + nx] == 0) {
                            pLabel[ny * iWidth + nx] = label;
                            q.push({nx, ny});
                        }
                    }
                }

                ++label;
            }
        }
    }

    std::vector<int> region_sizes(label, 0);

    for (int i = 0; i < iWidth * iHeight; ++i) {
        if (pLabel[i] > 0) {
            region_sizes[pLabel[i]]++;
        }
    }

    std::vector<std::pair<int, int>> regions;
    for (int i = 1; i < label; ++i) {
        regions.push_back({region_sizes[i], i});
    }

    std::sort(regions.begin(), regions.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {return a.first > b.first;});

    std::unordered_set<int> top_labels;

    int max_zones = 100;
    for (int i = 0; i < std::min(max_zones, static_cast<int>(regions.size())); ++i) {
        top_labels.insert(regions[i].second); 
    }

    for (int i = 0; i < iWidth * iHeight; ++i) {
        if (top_labels.find(pLabel[i]) == top_labels.end()) pLabel[i] = 0;
    }

    std::unordered_map<int, int> new_labels;
    int new_label = 1;

    for (int label : top_labels) {
        new_labels[label] = new_label++;
    }

    for (int i = 0; i < iWidth * iHeight; ++i) {
        if (pLabel[i] > 0) {
            pLabel[i] = new_labels[pLabel[i]];
        }
    }

    for (int i = 0; i < iWidth * iHeight; ++i) {
        if (pLabel[i] > 0) pLabel[i] = static_cast<unsigned char>(2 * pLabel[i]);
    }


}