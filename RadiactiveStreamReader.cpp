#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <random>
#include "DetEventInfo.hh"
using namespace std;


#define PI 3.14159265

int main(){


ifstream ifile("/media/simulation/Windows/StreamsGeometrie/1274er_Tub_2Rx1_Au_sample_12.12.22.stream", ios::binary);
ifstream ifile511("/media/simulation/Windows/StreamsGeometrie/511er_Tub_2Rx1_Au_sample_12.12.22.stream", ios::binary);
fstream myfile("/home/simulation/Schreibtisch/StreamReader/Data/Cylinder/Test.8181mm.txt",ios::out);
fstream myfile1("/home/simulation/Schreibtisch/StreamReader/Data/Random_exponential.txt",ios::out);
fstream myfile2("/home/simulation/Schreibtisch/StreamReader/Data/AbstrahlWinkel_alle.txt",ios::out);
fstream myfile3("/home/simulation/Schreibtisch/StreamReader/Data/IRF_Pyramide_4x4cm_1cm_2x2cm_50Mio.txt",ios::out);

// Vector for theoLifetime

const float bucket_size = 0.001;
int number_of_buckets = (int)ceil(15/ bucket_size); // requires <cmath>
std::vector<int> histogram_theolifetime(number_of_buckets);

// Vector for Lifetime

const float bucket_size1 = 0.005;
int number_of_buckets1 = (int)ceil(35/ bucket_size1); // requires <cmath>
std::vector<int> histogram_lifetime(number_of_buckets1);

// Energy1 and Energy2

const float bucket_size2 = 1.;
int number_of_buckets2 = (int)ceil(2000/ bucket_size2); // requires <cmath>
std::vector<int> histogram_Energy1(number_of_buckets2);


const float bucket_size3 = 0.001;
int number_of_buckets3 = (int)ceil(2/ bucket_size3);
std::vector<int> histogram_Energy2(number_of_buckets2);

int eventIndex = 0;

int rejectCnt = 0;
int Hit10 = 0.;
int HitInteraction10 = 0;
int Hit20 = 0.;
int HitInteraction20 = 0;
int Hit30 = 0.;
int HitInteraction30 = 0;
int Hit40 = 0.;
int HitInteraction40 = 0;
int Hit50 = 0.;
int HitInteraction50 = 0;
int Hit60 = 0.;
int HitInteraction60 = 0;
int Hit70 = 0.;
int HitInteraction70 = 0;
int Hit80 = 0.;
int HitInteraction80 = 0;
int Hit90 = 0.;
int HitInteraction90 = 0;

double Before1274kevTime_1 = 0.;
double Before1274kevTime_2 = 0.;
int TimeWindowCounter = 0;

int Hits1_before = 0;
int Hits2_before = 0;

int Window200ns_511_counter = 0.;
int Window200ns_1274_counter = 0.;

int EventsInBothDetektors = 0;
int ValidEvents = 0;
int False1274keV = 0;

int PartFalse1274keV = 0;
int PartFalse511keV = 0;
if(ifile.is_open())
{
 


    while(!ifile.eof() && !ifile511.eof()){
    
    
    EventInfo info;

    
    
    EventInfo info511;

    ifile.read((char*) &info, sizeof(EventInfo));
    ifile511.read((char*) &info511, sizeof(EventInfo));
   
   
  // Define the decay rates and intensities
    double decay_rate1 = 0.158;
    double intensity1 = 0.825;
    double decay_rate2 = 0.380;
    double intensity2 = 0.172;
    double decay_rate3 = 2.75;
    double intensity3 = 0.003;

    // Compute the overall intensity (must sum up to 1)
    double overall_intensity = intensity1 + intensity2 + intensity3;
    double tolerance = 1e-9; // set a small tolerance value
    if (std::abs(overall_intensity - 1.0) > tolerance) {
        std::cerr << "Error: The intensities must sum up to 1!" << std::endl;
        return 1;
    }

    // Define the random number generator
    std::mt19937 gen(std::random_device{}());

    // Generate a random number from the combined exponential distribution
    std::exponential_distribution<double> exp_dist(overall_intensity);
    double num = 0.0;
    while (num == 0.0) {
        double r = exp_dist(gen);
        if (r <= intensity1 / overall_intensity) {
            std::exponential_distribution<double> exp_dist1(decay_rate1);
            num = exp_dist1(gen);
        } else if (r <= (intensity1 + intensity2) / overall_intensity) {
            std::exponential_distribution<double> exp_dist2(decay_rate2);
            num = exp_dist2(gen);
        } else {
            std::exponential_distribution<double> exp_dist3(decay_rate3);
            num = exp_dist3(gen);
        }
    }

    // Print the result
    //std::cout << "Random number from the combined exponential distribution: " << num << std::endl;

    int uCi = 0;
    double activity = 37000*uCi;

    // Calculate the decay constant based on the activity level
    double decay_constant = activity / std::log(2);

    // Create a random number generator
    std::random_device rd1;
    std::mt19937 gen2(rd1());
    std::exponential_distribution<double> distribution(1 / decay_constant);
    double DecayTime = (1/distribution(gen2))*1e9;
    //std::cout << "Random number from the DecayTime: " << DecayTime << std::endl;
    


        const double Energy1 = info.info1().energy();
        const double Energy2 = info.info2().energy();
        const double Time1 = info.info1().MeanArr();
        const double Time2 = info.info2().MeanArr();

        const int NumberOfCounts1 = info.info1().numberOfCounts();
        const int NumberOfCounts2 = info.info2().numberOfCounts();
        

        const double Energy1_511 = info511.info1().energy();
        const double Energy2_511 = info511.info2().energy();
        const double Time1_511 = info511.info1().MeanArr();
        const double Time2_511 = info511.info2().MeanArr();

        const int NumberOfCounts1_511 = info511.info1().numberOfCounts();
        const int NumberOfCounts2_511 = info511.info2().numberOfCounts(); 

        const double GlTime1274 = info.startTime();
        const double GlTime511 = info511.startTime();
            if(myfile.is_open())
            {
            eventIndex ++;
            

               //Angle
                    const float StartAngleX = info.StartWinkelX();
                    const float StartAngleY = info.StartWinkelY();
                    const float StartAngleZ = info.StartWinkelZ();

                    const float StoppAngle1X = info.StoppWinkel1X();
                    const float StoppAngle1Y = info.StoppWinkel1Y();
                    const float StoppAngle1Z = info.StoppWinkel1Z();

                    const float StoppAngle2X = info.StoppWinkel2X();
                    const float StoppAngle2Y = info.StoppWinkel2Y();
                    const float StoppAngle2Z = info.StoppWinkel2Z();

                    const float AlphaStart = atan ((sqrt((StartAngleX*StartAngleX+StartAngleY*StartAngleY)))/StartAngleZ) * 180 / PI;
                    const float AlphaStopp = atan ((sqrt((StoppAngle1X*StoppAngle1X+StoppAngle1Y*StoppAngle1Y)))/StoppAngle1Z) * 180 / PI;
                    

                        
                    const bool HitDetektor10= ((AlphaStart<10.)&&(AlphaStart>=0.));
                    const bool HitDetektor20= ((AlphaStart<20.)&&(AlphaStart>=10.));
                    const bool HitDetektor30= ((AlphaStart<30.)&&(AlphaStart>=20.));
                    const bool HitDetektor40= ((AlphaStart<40.)&&(AlphaStart>=30.));
                    const bool HitDetektor50= ((AlphaStart<50.)&&(AlphaStart>=40.));
                    const bool HitDetektor60= ((AlphaStart<60.)&&(AlphaStart>=50.));
                    const bool HitDetektor70= ((AlphaStart<70.)&&(AlphaStart>=60.));
                    const bool HitDetektor80= ((AlphaStart<80.)&&(AlphaStart>=70.));
                    const bool HitDetektor90= ((AlphaStart<90.)&&(AlphaStart>=80.));

                    if(HitDetektor10){
                        Hit10++;
                        if(info.info1().isValid()){
                        HitInteraction10++;
                    }
                    }
                    if(HitDetektor20){
                        Hit20++;
                        if(info.info1().isValid()){
                        HitInteraction20++;
                    }
                    }
                    if(HitDetektor30){
                        Hit30++;
                        if(info.info1().isValid()){
                        HitInteraction30++;
                    }
                    }
                    if(HitDetektor40){
                        Hit40++;
                        if(info.info1().isValid()){
                        HitInteraction40++;
                    }
                    }
                    if(HitDetektor50){
                        Hit50++;
                        if(info.info1().isValid()){
                        HitInteraction50++;
                    }
                    }
                    if(HitDetektor60){
                        Hit60++;
                        if(info.info1().isValid()){
                        HitInteraction60++;
                    }
                    }
                    if(HitDetektor70){
                        Hit70++;
                        if(info.info1().isValid()){
                        HitInteraction70++;
                    }
                    }
                    if(HitDetektor80){
                        Hit80++;
                        if(info.info1().isValid()){
                        HitInteraction80++;
                    }
                    }
                    if(HitDetektor90){
                        Hit90++;
                        if(info.info1().isValid()){
                        HitInteraction90++;
                    }
                    }

                //myfile2<<StartAngleX<<" "<<StartAngleY<<" "<<StartAngleZ<<endl; 
                

                if(DecayTime< 200.){
                        if( (info511.info1().isValid() ||
                            info511.info2().isValid())){
                        Window200ns_511_counter ++;
                            }
                    }

                if(DecayTime< 200.){
                        if( (info.info1().isValid() ||
                            info.info2().isValid())){
                        Window200ns_1274_counter ++;
                            }
                    }

                if (((info.info1().isValid() &&
                 info.info2().isValid()) && (info511.info1().isValid() ||
                 info511.info2().isValid()))&&(DecayTime> 200))//
                { 
                PartFalse1274keV ++;
                }
                if (((info.info1().isValid() ||
                 info.info2().isValid()) && (info511.info1().isValid() &&
                 info511.info2().isValid()))&&(DecayTime> 200)) 
                { 
                PartFalse511keV ++;
                }
                if (((info.info1().isValid() ||
                 info.info2().isValid()) && (info511.info1().isValid() ||
                 info511.info2().isValid()))) 
                {
                    int Hits1 = info.info1().numberOfCounts();

                    int Hits2 = info.info2().numberOfCounts();

                    const bool isIn511WindowFalsePositive = ((Hits2>=610)&&(Hits2<1060))||((Hits1>=610)&&(Hits1<1060));

                    const bool isIn1274WindowBefore = ((Hits1_before>=1710)&&(Hits1_before<3300))||((Hits2_before>=1710)&&(Hits2_before<3300));

                    if(isIn1274WindowBefore){
                    if((Before1274kevTime_1<200.) || (Before1274kevTime_1<200.)){
                        if (isIn511WindowFalsePositive){
                               TimeWindowCounter ++;
                           }
                        }
                    }
                    



                    Before1274kevTime_1 = DecayTime +(info.info1().MeanArr()-GlTime1274);
                    Before1274kevTime_2 = DecayTime +(info.info2().MeanArr()-GlTime1274);
                    
                    Hits1_before = info.info1().numberOfCounts();

                    Hits2_before = info.info2().numberOfCounts();
                } 

                if (((info.info1().isValid() ||
                 info.info2().isValid()) && (info511.info1().isValid() ||
                 info511.info2().isValid())))//&&(DecayTime> 200) 
                {  
                //myfile1<<num<<endl;
                

                   EventsInBothDetektors ++;
                   // const bool isInEnergy1Window = ((Energy1>0.17)&&(Energy1<0.4))&&((Energy2>0.56)&&(Energy2<1.1));
                    
                    // const bool isInEnergy2Window = ((Energy2>0.17)&&(Energy2<0.4))&&((Energy1>0.56)&&(Energy1<1.1));

                    const bool isInEnergy1Window = ((Energy1>=0.17)&&(Energy1<0.4))&&((Energy2>=0.56)&&(Energy2<1.1));
                    
                    const bool isInEnergy2Window = ((Energy2>=0.17)&&(Energy2<0.4))&&((Energy1>=0.56)&&(Energy1<1.1));

                    const bool isInEnergyWindow = (isInEnergy1Window || isInEnergy2Window);

                    
                    
                    
                    const bool isInEnergy1Window_511 = ((Energy1_511>=0.17)&&(Energy1_511<1.35))&&((Energy2_511>=0.17)&&(Energy2_511<1.35));
                    
                    const bool isInEnergy2Window_511 = ((Energy2_511>=0.17)&&(Energy2_511<1.35))&&((Energy1_511>=0.17)&&(Energy1_511<1.35));

                    const bool isInEnergyWindow_511 = (isInEnergy1Window_511 || isInEnergy2Window_511);

                    
                    int Hits1 = info.info1().numberOfCounts();

                    int Hits2 = info.info2().numberOfCounts();

                    int Hits1_511 = info511.info1().numberOfCounts();

                    int Hits2_511 = info511.info2().numberOfCounts();

                    const bool isInPulsWindow1 =((Hits1>=610)&&(Hits1<1060))||((Hits2>=1710)&&(Hits2<3300));
                    const bool isInPulsWindow2 =((Hits2>=610)&&(Hits2<1060))||((Hits1>=1710)&&(Hits1<3300));
                    const bool isInPulsWindow = (isInPulsWindow1 || isInPulsWindow2);
                    const bool isIn511WindowFalsePositive = ((Hits2>=610)&&(Hits2<1060))||((Hits1>=610)&&(Hits1<1060));

                    const bool isInPulsWindow1_511 =((Hits1_511>=610)&&(Hits1_511<1060))||((Hits2_511>=1710)&&(Hits2_511<3300));
                    const bool isInPulsWindow2_511 =((Hits2_511>=610)&&(Hits2_511<1060))||((Hits1_511>=1710)&&(Hits1_511<3300));
                    const bool isInPulsWindow_511 = (isInPulsWindow1_511 || isInPulsWindow2_511);


                    const bool isInRightWindow_1 = ((Hits1_511>=610)&&(Hits1_511<1060))||((Hits1>=1710)&&(Hits1<3300));
                    const bool isInRightWindow_2 = ((Hits1_511>=610)&&(Hits1_511<1060))||((Hits2>=1710)&&(Hits2<3300));
                    const bool isInRightWindow_3 = ((Hits2_511>=610)&&(Hits2_511<1060))||((Hits1>=1710)&&(Hits1<3300));
                    const bool isInRightWindow_4 = ((Hits2_511>=610)&&(Hits2_511<1060))||((Hits2>=1710)&&(Hits2<3300));
                    const bool isInRightWindow = (isInRightWindow_1 || isInRightWindow_2 || isInRightWindow_3 || isInRightWindow_1);

                    const bool TrueCoincidence = (isInRightWindow_2 || isInRightWindow_3 );

                    
                 
                    //Interaction Positions
                    const double Position1X = info.info1().InterActionX();
                    const double Position1Y = info.info1().InterActionY();
                    const double Position1Z = info.info1().InterActionZ();

                    const double Position2X = info.info2().InterActionX();
                    const double Position2Y = info.info2().InterActionY();
                    const double Position2Z = info.info2().InterActionZ();

                    const double InteractionLength1 = sqrt((Position1X*Position1X) +(Position1Y*Position1Y) +(Position1Z*Position1Z));
                    const double InteractionLength2 = sqrt((Position2X*Position2X) +(Position2Y*Position2Y) +(Position2Z*Position2Z));
                    //const double Lifetime_with_maxPulse = info.info1().MaxPulse() - info.info2().MaxPulse();
                    
                    const bool isInTimeWindow = ((info.info1().MeanArr()>0.)&&(info.info2().MeanArr()>0.)); // ns

                    const bool isInTimeWindow_511 = ((info511.info1().MeanArr()>0.)&&(info511.info2().MeanArr()>0.)); // ns

                    const double  TheoLifetime = info.lifetime();

                    
                    const double Lifetime_with_mean_1 = ((info511.info1().MeanArr()-GlTime511) - (info.info1().MeanArr()-GlTime1274));
                    const double Lifetime_with_mean_2 = ((info511.info2().MeanArr()-GlTime511) - (info.info1().MeanArr()-GlTime1274));
                    const double Lifetime_with_mean_3 = ((info511.info1().MeanArr()-GlTime511) - (info.info2().MeanArr()-GlTime1274));
                    const double Lifetime_with_mean_4 = ((info511.info2().MeanArr()-GlTime511) - (info.info2().MeanArr()-GlTime1274));

                    const bool TrueLifetime_1 = (info511.info1().isValid() && info.info1().isValid());
                    const bool TrueLifetime_2 = (info511.info2().isValid() && info.info1().isValid());
                    const bool TrueLifetime_3 = (info511.info1().isValid() && info.info2().isValid());
                    const bool TrueLifetime_4 = (info511.info2().isValid() && info.info2().isValid());
                    
                    

                       if ( isInPulsWindow){
                           if (isIn511WindowFalsePositive){
                               False1274keV ++;
                           }
                           
                       }
                       if ( isInPulsWindow_511){
                           
                           
                       }
                          if (TrueCoincidence)   
                        {
                            ValidEvents ++;
                        }
                        /*if (isInRightWindow)   
                        {
                            if(TrueLifetime_1){
                               myfile3<< Lifetime_with_mean_1<<endl;
                            }
                            if(TrueLifetime_2){
                               myfile3<< Lifetime_with_mean_2<<endl;
                            }
                            if(TrueLifetime_3){
                               myfile3<< Lifetime_with_mean_3<<endl;
                            }
                            if(TrueLifetime_4){
                               myfile3<< Lifetime_with_mean_4<<endl;
                            } 
                         
                           
                            
                        }*/
                        /*if(isInEnergy2Window)// && isInTimeWindow)
                        {
                        
                             
                             
                          
                            int bucket3 = (int)round(AbsLifetime_with_mean / bucket_size1);
                            //histogram_lifetime[bucket3] += 1;
                            
                             if (StartAngleZ>0){
                                //myfile<<AlphaStart<<" "<<Position2X<<" "<<Position2Y<<" "<<Position2Z<<" "<<Energy2<<" "<<NumberOfCounts2<<" "<<InteractionLength2<<endl;
                            }
                           

                        }*/
                        else{
                            rejectCnt ++;}
                }       
                else{
                rejectCnt ++;}
            }   
            else
            {
            std::cout<<"Error"<<endl;
            } 
   
    
}
   
    


    std::vector<int>::iterator itr;
    /*for (itr=histogram_theolifetime.begin();itr!=histogram_theolifetime.end();itr++)
    {
        myfile<<*itr<<endl;

    }*/
   /* for (itr=histogram_lifetime.begin();itr!=histogram_lifetime.end();itr++)
    {
        myfile3<<*itr<<endl;

    }
    /*for (itr=histogram_Energy1.begin();itr!=histogram_Energy1.end();itr++)
    {
        myfile<<*itr<<endl;

    }
    for (itr=histogram_Energy2.begin();itr!=histogram_Energy2.end();itr++)
    {
        myfile3<<*itr<<endl;

    }*/
}
double Effizenz10 = (double)HitInteraction10/(double)Hit10;
double Effizenz20 = (double)HitInteraction20/(double)Hit20;
double Effizenz30 = (double)HitInteraction30/(double)Hit30;
double Effizenz40 = (double)HitInteraction40/(double)Hit40;
double Effizenz50 = (double)HitInteraction50/(double)Hit50;
double Effizenz60 = (double)HitInteraction60/(double)Hit60;
double Effizenz70 = (double)HitInteraction70/(double)Hit70;
double Effizenz80 = (double)HitInteraction80/(double)Hit80;
double Effizenz90 = (double)HitInteraction90/(double)Hit90;

    std::cout<<"Number of Events :"<<eventIndex<<endl;
    std::cout<<"Number of Events in Spectrum :"<<ValidEvents<<endl;
    std::cout<<"Number of Events in both detectors :"<<EventsInBothDetektors<<endl;
    std::cout<<"Number of 1274kev in 511kev Window:"<<False1274keV<<endl;
    std::cout<<"Number of 1274kev in 511kev Window in the 200ns Window:"<<TimeWindowCounter<<endl;
    std::cout<<"Number of 511kev in 200ns Window:"<<Window200ns_511_counter<<endl;
    std::cout<<"Number of 1274kev in 200ns Window:"<<Window200ns_1274_counter<<endl;
    std::cout<<"PartFalse1274keV in coincidences:"<<PartFalse1274keV<<endl;
    std::cout<<"PartFalse511keV in coincidences:"<<PartFalse511keV<<endl;
    /*std::cout<<"Effizienz10° "<<Effizenz10<<std::endl;
    std::cout<<"Effizienz20° "<<Effizenz20<<std::endl;
    std::cout<<"Effizienz30° "<<Effizenz30<<std::endl;
    std::cout<<"Effizienz40° "<<Effizenz40<<std::endl;
    std::cout<<"Effizienz50° "<<Effizenz50<<std::endl;
    std::cout<<"Effizienz60° "<<Effizenz60<<std::endl;
    std::cout<<"Effizienz70° "<<Effizenz70<<std::endl;
    std::cout<<"Effizienz80° "<<Effizenz80<<std::endl;
    std::cout<<"Effizienz90° "<<Effizenz90<<std::endl;
    std::cout<<"hits 10° "<<Hit10<<std::endl;
    std::cout<<"hits 20° "<<Hit20<<std::endl;
    std::cout<<"hits 30° "<<Hit30<<std::endl;
    std::cout<<"hits 40° "<<Hit40<<std::endl;
    std::cout<<"hits 50° "<<Hit50<<std::endl;
    std::cout<<"hits 60° "<<Hit60<<std::endl;
    std::cout<<"hits 70° "<<Hit70<<std::endl;
    std::cout<<"hits 80° "<<Hit80<<std::endl;
    std::cout<<"hits 90° "<<Hit90<<std::endl;*/

    

    ifile.close();
    ifile511.close();
    myfile.close();
    myfile1.close();
    myfile2.close();
    myfile3.close();

  
    return 0;

}
