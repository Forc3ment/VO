#define VP_TRACE

#include <string>
#include <fstream>

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

int imgHeight = 112;
int imgWidth = 92;

string toString(int i) // convert int to string
{
    std::stringstream value;
    value << i;
    return value.str();
}

void saveEigenFaces(const vpMatrix & U, const vpColVector & w, int nbToSave, bool display = false)
{        
    for (int i = 0; i < nbToSave; ++i)
    {
        vpImage<double> temp1(imgHeight,imgWidth);
        vpImage<unsigned char> temp2;

        for(int pointBit = 0; pointBit < U.getRows(); pointBit++)
        {
            temp1.bitmap[pointBit] = U[pointBit][i];
        }

        vpImageConvert::convert(temp1, temp2);
        
        ofstream fichier("../result/eigenFaceValue.csv", ios::out | ios::app);  //déclaration du flux et ouverture du fichier
        if(fichier)
        {
            fichier << w[i] << endl;
        }
        vpImageIo::write(temp2,"../result/eigenFace"+toString(i)+".pgm");
        if(display)
        {
            vpDisplayX d3(temp2);
            vpDisplay::display(temp2);
            vpDisplay::flush(temp2);
            vpDisplay::getClick(temp2);
        }
        if(fichier)
        {
            fichier.close();
        }
    }
}

vpColVector computeMeanFaces(const vpMatrix & I, bool display, int displayNumber=0)
{
    vpColVector meanFace(I.getRows(),0);

    for (int n = 0; n < I.getRows(); ++n)
    {
        for (int m = 0; m < I.getCols(); ++m)
        {
            meanFace[n] += I[n][m];
        }
    } 

    meanFace/=I.getCols();

    vpImage<double> temp1(imgHeight,imgWidth);
    vpImage<unsigned char> temp2;

    for(int pointBit = 0; pointBit < I.getRows(); pointBit++)
    {
        temp1.bitmap[pointBit] = meanFace[pointBit];
    }

    vpImageConvert::convert(temp1, temp2);
    vpImageIo::write(temp2,"../result/meanFace.pgm");

    vpImage<double> temp3(imgHeight,imgWidth);
    vpImage<unsigned char> temp4;


    for(int pointBit = 0; pointBit < I.getRows(); pointBit++)
    {
        temp3.bitmap[pointBit] = I[pointBit][displayNumber];
    }

    vpImageConvert::convert(temp3, temp4);
    vpImageIo::write(temp4,"../result/face"+toString(displayNumber)+".pgm");


    vpImage<double> temp5(imgHeight,imgWidth);
    vpImage<unsigned char> temp6;

    for(int pointBit = 0; pointBit < I.getRows(); pointBit++)
    {
        temp5.bitmap[pointBit] = I[pointBit][displayNumber] - meanFace[pointBit];
    }

    vpImageConvert::convert(temp5, temp6);
    vpImageIo::write(temp6,"../result/faceCentered"+toString(displayNumber)+".pgm");

    if(display)
    {
        vpDisplayX d1(temp2);
        vpDisplay::display(temp2);
        vpDisplay::flush(temp2);
        

        vpDisplayX d2(temp4);
        vpDisplay::display(temp4);
        vpDisplay::flush(temp4);
        

        vpDisplayX d3(temp6);
        vpDisplay::display(temp6);
        vpDisplay::flush(temp6);
        vpDisplay::getClick(temp6);
    }

    return meanFace;
}

vpMatrix computeA(const vpMatrix & I, const vpColVector & meanFace)
{
    vpMatrix A(I.getRows(), I.getCols());

    for (int m = 0; m < I.getCols(); ++m)
    {
        for(int n = 0; n < I.getRows(); n++)
        {
            A[n][m] = I[n][m] - meanFace[n];
        }
    }

    return A;
}

vpColVector computeWk(const vpColVector & J, const vpMatrix & U, const vpColVector & meanFace, int k)
{
    vpColVector Wk(k);

    for (int i = 0; i < k; ++i)
    {
        vpColVector Uk(U.getRows());
        for (int j = 0; j < U.getRows(); ++j)
        {
            Uk[j] = U[j][i];
        }
        
        Wk[i] = Uk * (J-meanFace);
    }                                                                              
    return Wk;
}

vpColVector computeReconstruction(const vpMatrix & U, const vpColVector & meanFace, const vpColVector & Wk, string folderName, string imgNumber, string k, bool display=false)
{
    vpColVector Wkuk(meanFace.size());
    for (int i = 0; i < Wk.size(); ++i)
    {
        vpColVector Uk(U.getRows());
        for (int j = 0; j < U.getRows(); ++j)
        {
            Uk[j] = U[j][i];
        }
        Wkuk += Wk[i]*Uk;     
    }
    vpColVector Jp = meanFace + Wkuk;

    // vpImage<double> temp1(imgHeight,imgWidth);
    // vpImage<unsigned char> temp2;

    // for(int i = 0; i < U.getRows(); i++)
    // {
    //     temp1.bitmap[i] = Jp[i];
    // }

    //vpImageConvert::convert(temp1, temp2);
        
    // ofstream fichier("../result/Wk.txt", ios::out | ios::app);  
    // if(fichier)
    // {
    //     fichier << Wk << "\n" << endl;
    //     fichier.close();
    // }
    //vpImageIo::write(temp2,"../result/reconstructedFace"+folderName+"_"+imgNumber+"_k"+k+".pgm");
    // if(display)
    // {
    //     vpDisplayX d3(temp2);
    //     vpDisplay::display(temp2);
    //     vpDisplay::flush(temp2);
    //     vpDisplay::getClick(temp2);
    // }

    return Jp;
}

double computeReconstructionError(const vpColVector & J, const vpColVector & Jp)
{    
    //double toto = sqrt((J-Jp).sumSquare());
    double toto = 0;
    for (int i = 0; i < Jp.size(); ++i)
    {
        toto += (J[i] - Jp[i]) * (J[i] - Jp[i]);
    }

    return sqrt(toto);
}

void identify(const vector<vpColVector> wkRef, const vpMatrix U, const vpColVector meanFace, int k, string fileName)
{
    vpImage<unsigned char> img;

    vpImageIo::read(img,fileName);    
    vpDisplayX d(img);
    vpDisplay::display(img);
    vpDisplay::flush(img);
    vpDisplay::getClick(img);

    int size = img.getWidth() * img.getHeight();
    vpColVector JToIdentify(size);

    for (int j = 0; j < size; ++j)
    {
        JToIdentify[j] = (double)(img.bitmap[j])/255;
    }

    vpColVector Wk = computeWk(JToIdentify,U,meanFace,k);

    int counter2 = 0;
    double min = 99999;
    int saveFileNum, savePicNum;

    double error;

    for(int fileNum = 0; fileNum < 40; fileNum++)
    {
       for(int picNum = 0; picNum < 5; picNum++)
        {
            error = computeReconstructionError(wkRef[counter2],Wk);
            if(error < min)
            {
                min = error;
                saveFileNum = fileNum;
                savePicNum = picNum;
            }

            counter2++;
        }
    }

    vpImage<unsigned char> img2;
    string fileNameBis = "../faces/s" + toString(saveFileNum+1) +"/" + toString(savePicNum+1) +".pgm";
    vpImageIo::read(img2,fileNameBis);

    vpDisplayX d3(img2);
    vpDisplay::display(img2);
    vpDisplay::flush(img2);
    vpDisplay::getClick(img2);
}

void identifyAll(const vector<vpColVector> wkRef, const vpMatrix U, const vpColVector meanFace, int k)
{
    vpImage<unsigned char> img;
    vpImage<double> error(wkRef[1].size(),wkRef[1].size());
    int counter1 = 0;
    int sum =0;

    for(int fileNum = 0; fileNum < 40; fileNum++)
    {
        for(int picNum = 5; picNum < 10; picNum++)
        {
            string fileName = "../faces/s" + toString(fileNum+1) +"/" + toString(picNum+1) +".pgm";
            vpImageIo::read(img,fileName);

            int size = img.getWidth() * img.getHeight();

            vpColVector JToIdentify(size);

            for (int j = 0; j < size; ++j)
            {
                JToIdentify[j] = (double)(img.bitmap[j])/255;
            }

            vpColVector Wk = computeWk(JToIdentify,U,meanFace,k);

            int counter2 = 0;
            double minInter = 999999;
            double minExter = 999999;

            for(int fileNum2 = 0; fileNum2 < 40; fileNum2++)
            {
                for(int picNum2 = 0; picNum2 < 5; picNum2++)
                {
                    error[counter1][counter2] = computeReconstructionError(wkRef[counter2],Wk);
                    if(fileNum == fileNum2)
                    {
                        minInter = min(minInter, error[counter1][counter2]);
                    }
                    else
                    {
                        minExter = min(minExter, error[counter1][counter2]);
                    }

                    counter2++;
                }
            }
            ofstream fichier("../result/identification.csv", ios::out | ios::app);  //déclaration du flux et ouverture du fichier
            if(fichier)
            {
                int id = (minInter < minExter) ? 1 : 0;
                sum += id;
                fichier << id << endl;
            }
            counter1++;
        }
    }
    vpImage<unsigned char> errorUC;
    vpImageConvert::convert(error, errorUC);
    vpImageIo::write(errorUC,"../result/error.pgm");

    cout << sum << endl;
}


int main(int argc, char** argv)
{
    // Déclaration des variables
    int size = imgWidth * imgHeight;
    int nbImg = 200;
    int count = 0;
    int k = 200;
    vpMatrix I(size, nbImg);
    vpImage<unsigned char> img;
    string fileName;

    // Boucle de construction de la matrice I
    for(int fileNum = 0; fileNum < 40; fileNum++)
    {
        for(int picNum = 0; picNum < 5; picNum++)
        {
            fileName = "../faces/s" + toString(fileNum+1) +"/" + toString(picNum+1) +".pgm";

            vpImageIo::read(img,fileName);

            for(int pointBit = 0; pointBit < size; pointBit++)
            {
                I[pointBit][count] = (double)(img.bitmap[pointBit])/255;
            }
            count++;
        }
    }

    vpColVector meanFace = computeMeanFaces(I,false);

    vpMatrix A = computeA(I, meanFace);

    vpMatrix U = A;
    vpColVector w;
    vpMatrix V;
    
    cout << "------- Start SVD -------" << endl;
    U.svd(w, V);
    cout << "-------  END SVD  -------" << endl;

    //saveEigenFaces(U,w,200);

    cout << "------- Start Construction of the base -------" << endl;

    int fileNum = 0;
    int picNum = 0;
    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;

    
    vector<vpColVector> wkRef;

    for(int fileNum2 = 0; fileNum2 < 40; fileNum2++)
    {
        for(int picNum2 = 0; picNum2 < 5; picNum2++)
        {
            fileName = "../faces/s" + toString(fileNum2+1) +"/" + toString(picNum2+1) +".pgm";
            vpImageIo::read(img,fileName);
            vpColVector Jref(size);

            for (int j = 0; j < size; ++j)
            {
                Jref[j] = (double)(img.bitmap[j])/255;
            }

            wkRef.push_back(computeWk(Jref,U,meanFace,k));
        }
    }

    cout << "-------  END of the construction -------" << endl;

    int a = 0;

    while(a >= 0)
    {
        cout << "1 :  Identify All"<< endl;
        cout << "2 :  Identify One"<< endl;
        cout << "3 :  Change K"<< endl;
        cout << "k = " << k << endl;
        cout << endl;
        cout << "What do you want to do ?"<< endl;
        cout << endl;
        cin >> a;

       switch (a)
        {
            case 1 :
                identifyAll(wkRef, U, meanFace, k);
                break;
            case 2 :
                cout << "Enter a picture location."<< endl;
                cin >> fileName;
                //fileName = "../faces/s11/8.pgm";
                identify(wkRef, U, meanFace, k, fileName);
                break;
            case 3 : 
                cout << "Enter a value for k."<< endl;
                cin >> k;
                break;
            default:
                cout << "Wrong proposition."<< endl;
                break;
        }
        cout << "-----------------------------------------------------"<< endl; 
    }
    
    
    
    
}
