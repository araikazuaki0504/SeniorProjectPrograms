#include <math.h>
#include <stdio.h>

#include "DataToImage.hpp"

using namespace std;

DataToImage::DataToImage(double* Data, unsigned int Image_width, unsigned int Image_height):_Data{Data},_Image_width{Image_width},_Image_height{Image_height}
{
     _Image_matrix = cv::Mat(cv::Size(Image_width,Image_height),CV_8UC3,cv::Scalar(255,255,255));
};

DataToImage& DataToImage::changeToImage()
{
    //大きい円の描画
    cv::circle(_Image_matrix,
                cv::Point(_Image_width / 2,  _Image_height / 2), 
                _LargeCircleRadius - _LargeCircleLineThickness - _SmallCircleLineThickness,
                _LineColor,
                _LargeCircleLineThickness,
                cv::LINE_4);

    //小さい円の描画
    for(int i = 0; i < _Number_of_DataElement; i++)
    {
        int Blue = (int)(_Data[_fiber_index_mapping_table[i] - 1] * (double)255);
        int Red = (int)(((double)1 - _Data[_fiber_index_mapping_table[i] - 1]) * (double)255);
        cv::Scalar Fill_BGRCOLOR = cv::Scalar(Blue, 0, Red);
        
        char fiberNumber[3];
        if(_fiber_index_mapping_table[i] < 10)snprintf(fiberNumber,3,"0%d",_fiber_index_mapping_table[i]);
        else if (_fiber_index_mapping_table[i] >= 10)snprintf(fiberNumber,3,"%d",_fiber_index_mapping_table[i]);
        cv::Size Text_size = cv::getTextSize(
            fiberNumber,
            _Text_font,
            _Text_scale,
            _Text_thickness,
            &_Text_baseLine);

        if(i == 0)
        {
            //外枠
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2, _Image_height / 2), 
                _SmallCircleRadius,
                _LineColor,
                _SmallCircleLineThickness,
                cv::LINE_4);

            //塗りつぶし
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2,  _Image_height / 2), 
                _SmallCircleRadius,
                Fill_BGRCOLOR,
                _fill,
                cv::LINE_4);
            
            //番号
            cv::putText(_Image_matrix,
                fiberNumber,
                cv::Point(( _Image_width / 2) - (Text_size.width / 2), (_Image_height / 2) + (Text_size.height / 2)),
                _Text_font,
                _Text_scale,
                _LineColor
            );
        }
        else if (0 < i && i <= 6)
        {
            //外枠
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2 + 2 * _SmallCircleRadius * cos(2 * M_PI * (i - 1) / _Number_of_innerElement),  _Image_height / 2 + 2 * _SmallCircleRadius * sin(2 * M_PI * (i - 1) / _Number_of_innerElement)), 
                _SmallCircleRadius,
                _LineColor,
                _SmallCircleLineThickness,
                cv::LINE_4);
            
            //塗りつぶし
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2 + 2 * _SmallCircleRadius * cos(2 * M_PI * (i - 1) / _Number_of_innerElement),  _Image_height / 2 + 2 * _SmallCircleRadius * sin(2 * M_PI * (i - 1) / _Number_of_innerElement)), 
                _SmallCircleRadius,
                Fill_BGRCOLOR,
                _fill,
                cv::LINE_4);

            //番号
            cv::putText(_Image_matrix,
                fiberNumber,
                cv::Point(( _Image_width / 2 + 2 * _SmallCircleRadius * cos(2 * M_PI * (i - 1) / _Number_of_innerElement)) - (Text_size.width / 2), (_Image_height / 2 + 2 * _SmallCircleRadius * sin(2 * M_PI * (i - 1) / _Number_of_innerElement)) + (Text_size.height / 2)),
                _Text_font,
                _Text_scale,
                _LineColor
            );
        }
        else if (6 < i && i <= 18)
        {
            //外枠
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2 + 4 * _SmallCircleRadius * cos(2 * M_PI * (i - 7) / _Number_of_middleElement + M_PI / _Number_of_middleElement), _Image_height / 2 + 4 * _SmallCircleRadius * sin(2 * M_PI * (i - 7) / _Number_of_middleElement + M_PI / _Number_of_middleElement)), 
                _SmallCircleRadius,
                _LineColor,
                _SmallCircleLineThickness,
                cv::LINE_4);

            //塗りつぶし
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2 + 4 * _SmallCircleRadius * cos(2 * M_PI * (i - 7) / _Number_of_middleElement + M_PI / _Number_of_middleElement), _Image_height / 2 + 4 * _SmallCircleRadius * sin(2 * M_PI * (i - 7) / _Number_of_middleElement + M_PI / _Number_of_middleElement)), 
                _SmallCircleRadius,
                Fill_BGRCOLOR,
                _fill,
                cv::LINE_4);

            //番号
            cv::putText(_Image_matrix,
                fiberNumber,
                cv::Point(( _Image_width / 2 + 4 * _SmallCircleRadius * cos(2 * M_PI * (i - 7) / _Number_of_middleElement + M_PI / _Number_of_middleElement)) - (Text_size.width / 2), (_Image_height / 2  + 4 * _SmallCircleRadius * sin(2 * M_PI * (i - 7) / _Number_of_middleElement + M_PI / _Number_of_middleElement)) + (Text_size.height / 2)),
                _Text_font,
                _Text_scale,
                _LineColor
            );
        }
        else if (18 < i && i <= 36)
        {
            //外枠
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2 + 6 * _SmallCircleRadius * cos(2 * M_PI * (i - 19) / _Number_of_outerElement + M_PI / _Number_of_outerElement),  _Image_height / 2 + 6 * _SmallCircleRadius * sin(2 * M_PI *(i - 19) / _Number_of_outerElement + M_PI / _Number_of_outerElement)), 
                _SmallCircleRadius,
                _LineColor,
                _SmallCircleLineThickness,
                cv::LINE_4);

            //塗りつぶし
            cv::circle(_Image_matrix,
                cv::Point( _Image_width / 2 + 6 * _SmallCircleRadius * cos(2 * M_PI * (i - 19) / _Number_of_outerElement + M_PI / _Number_of_outerElement),  _Image_height / 2 + 6 * _SmallCircleRadius * sin(2 * M_PI *(i - 19) / _Number_of_outerElement + M_PI / _Number_of_outerElement)), 
                _SmallCircleRadius,
                Fill_BGRCOLOR,
                _fill,
                cv::LINE_4);

            //番号
            cv::putText(_Image_matrix,
                fiberNumber,
                cv::Point(( _Image_width / 2 + 6 * _SmallCircleRadius * cos(2 * M_PI * (i - 19) / _Number_of_outerElement + M_PI / _Number_of_outerElement)) - (Text_size.width / 2), (_Image_height / 2 + 6 * _SmallCircleRadius * sin(2 * M_PI *(i - 19) / _Number_of_outerElement + M_PI / _Number_of_outerElement)) + (Text_size.height / 2)),
                _Text_font,
                _Text_scale,
                _LineColor
            );
        }
        
        //char outputPath[31];
        //sprintf(outputPath, "../outputFiles/TestImage%d.jpg",i);
        //cout << outputPath << endl;
        //cv::imwrite(outputPath, _Image_matrix);
    }
    return *this;
}

cv::Mat DataToImage::getImage()
{
    return _Image_matrix;
}

Vector_2D DataToImage::getImage_Size()
{
    return {static_cast<int>(_Image_width),static_cast<int>(_Image_height)};
}

DataToImage& DataToImage::setLineColor(cv::Scalar LineColor)
{
    _LineColor = LineColor;
    return *this;
}

DataToImage& DataToImage::setText_height(int Text_height)
{
    _Text_height = Text_height;
    return *this;
}

DataToImage& DataToImage::setText_font(int Text_font)
{
    _Text_font = Text_font;
    return *this;
}

DataToImage& DataToImage::setText_thickness(int Text_thickness)
{ 
    _Text_thickness = Text_thickness;
    return *this;
}

DataToImage& DataToImage::saveImage(const char* outputfilePath)
{
    cv::imwrite(outputfilePath, _Image_matrix);
    cout << "Success save Image To " << outputfilePath << endl;
    return *this;
}