#define VP_TRACE

#include <string>

//! \example tutorial-viewer.cpp
//! [Include display]
#include <visp3/gui/vpDisplayD3D.h>
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayGTK.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/core/vpPoint.h>
#include <visp3/core/vpPlane.h>
#include <visp3/gui/vpPlot.h>
#include <visp3/core/vpMeterPixelConversion.h>
#include <visp3/core/vpCameraParameters.h>

#include <visp3/core/vpExponentialMap.h>
#include <visp3/core/vpVelocityTwistMatrix.h>



//! [Include display]
//! [Include io]
#include <visp3/io/vpImageIo.h>
//! [Include io]

using namespace std ;


string toString(int i) // convert int to string
{
    std::stringstream value;
    value << i;
    return value.str();
}

int main(int argc, char** argv)
{

    // Déclaration des variables
    int size = 92 * 112;
    int nbImg = 200;
    vpMatrix I(size, nbImg);
    vpImage<unsigned char> img;
    string fileName;

    // Boucle de construction de la matrice I
    for(int fileNum = 1; fileNum <= 40; fileNum++)
    {
        for(int picNum = 1; picNum <= 5; picNum++)
        {
            fileName = "../faces/s" + toString(fileNum) +"/" + toString(picNum) +".pgm";

            vpImageIo::read(img,fileName);

            for(int pointBit = 0; pointBit < size; pointBit++)
            {
                I[pointBit][(fileNum -1) * (picNum -1)] = (double)(img.bitmap[pointBit])/255;
            }
        }
    }

    vpColVector meanFace(size,0);

    for (int n = 0; n < size; ++n)
    {
        for (int m = 0; m < nbImg; ++m)
        {
            meanFace[n] += I[n][m];
        }
    } 

    meanFace/=nbImg;

    vpImage<double> temp1(112,92);
    vpImage<unsigned char> temp2;


    for(int pointBit = 0; pointBit < size; pointBit++)
    {
        temp1.bitmap[pointBit] = meanFace[pointBit];
    }

    vpImageConvert::convert(temp1, temp2);

    vpDisplayX d(temp2);
    vpDisplay::display(temp2);
    vpDisplay::flush(temp2);
    vpDisplay::getClick(temp2);

}
