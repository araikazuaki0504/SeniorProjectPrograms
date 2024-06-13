#pragma once

#include <opencv2/opencv.hpp>
#include <math.h>

#include "utils.hpp"

using namespace std;

class DataToImage{
    public:
        //Data : データ、Image_width : 画像の幅、Image_height : 画像の高さ
        DataToImage(double* Data, unsigned int Image_width, unsigned int Image_height);
        DataToImage& changeToImage();
        cv::Mat getImage();
        Vector_2D getImage_Size();
        DataToImage& setLineColor(cv::Scalar LineColor);
        DataToImage& setText_height(int Text_height);
        DataToImage& setText_font(int Text_font);
        DataToImage& setText_thickness(int Text_thickness);
        DataToImage& saveImage(const char* outputPath);
        ~DataToImage(){};

    private:
        double* _Data;
        unsigned int _Image_width;
        unsigned int _Image_height;
        cv::Mat _Image_matrix;
        cv::Scalar _LineColor = cv::Scalar(0,0,0);
        int _Text_height = 20;//pixel?
        int _Text_font = cv::FONT_HERSHEY_DUPLEX;
        int _Text_thickness = 1;//pixel?
        const float _LargeCircleLineThickness = 3;
        const float _SmallCircleLineThickness = 3;
        const int _fill = -1;
        const double _Text_scale = cv::getFontScaleFromHeight(_Text_font,_Text_height,_Text_thickness);
        int _Text_baseLine = 0;
        const double _Large_Small_RadiusRatio = 0.14795590447907633;//大きい円と小さい円との半径比(参考：http://www.packomania.com/)
        const double _LargeCircleRadius = (_Image_width > _Image_height ? _Image_height / 2 : _Image_width / 2);  
        const double _SmallCircleRadius = _Large_Small_RadiusRatio * _LargeCircleRadius - _SmallCircleLineThickness;
        const unsigned int _Number_of_DataElement   = 37;
        const unsigned int _Number_of_outerElement  = 18;
        const unsigned int _Number_of_middleElement = 12;
        const unsigned int _Number_of_innerElement   = 6;
        const unsigned int _fiber_index_mapping_table[37] = {
            20,//1
            19,//2
            21,//3
            17,//4
            16,//5
            22,//6
            18,//7
            25,//8
            2,//9
            12,//10
            33,//11
            34,//12
            1,//13
            11,//14
            15,//15
            30,//16
            35,//17
            37,//18
            24,//19
            7,//20
            8,//21
            23,//22
            9,//23
            29,//24
            13,//25
            36,//26
            31,//27
            28,//28
            4,//29
            6,//30
            14,//31
            26,//32
            32,//33
            5,//34
            27,//35          
            10,//36
            3,//37                    
        };
};