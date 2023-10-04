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


double deformed_gauss_pmt1(double transit_time_value, double transit_time_spread_value){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> selector(0.0, 1.0);
        double randomValue = selector(gen);
        double pmt_time_spread = 0.0;
        double sample;
        if (randomValue < 0.995) {
            // Generieren Sie eine Zufallszahl aus der ersten Gauß'schen Verteilung
            std::normal_distribution<double> distribution1(transit_time_value, transit_time_spread_value);
            pmt_time_spread = distribution1(gen);
        } else {
            // Generieren Sie eine Zufallszahl aus der zweiten Gauß'schen Verteilung
            std::normal_distribution<double> distribution2(transit_time_value + 5*transit_time_spread_value, 4*transit_time_spread_value);
            pmt_time_spread = distribution2(gen);
        }
    return pmt_time_spread;
}

double deformed_gauss_pmt2(double transit_time_value, double transit_time_spread_value){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> selector(0.0, 1.0);
        double randomValue = selector(gen);
        double pmt_time_spread= 0.0;
        double sample;
        if (randomValue < 0.995) {
            // Generieren Sie eine Zufallszahl aus der ersten Gauß'schen Verteilung
            std::normal_distribution<double> distribution1(transit_time_value, transit_time_spread_value);
            pmt_time_spread = distribution1(gen);
        } else {
            // Generieren Sie eine Zufallszahl aus der zweiten Gauß'schen Verteilung
            std::normal_distribution<double> distribution2(transit_time_value + 5*transit_time_spread_value, 4*transit_time_spread_value);
            pmt_time_spread = distribution2(gen);
        }
    return pmt_time_spread;
}

double gauss_pmt1(double transit_time_value, double transit_time_spread_value){
    random_device rd_pmt;
    mt19937 gen_pmt (rd_pmt());

    normal_distribution<> gauss {transit_time_value, transit_time_spread_value};

    double pmt_time_spread = gauss(gen_pmt);

    return pmt_time_spread;
}

double gauss_pmt2(double transit_time_value, double transit_time_spread_value){
    random_device rd_pmt;
    mt19937 gen_pmt (rd_pmt());

    normal_distribution<> gauss {transit_time_value, transit_time_spread_value};

    double pmt_time_spread = gauss(gen_pmt);

    return pmt_time_spread;
}

double NeonDecay(double mean){
    random_device rd_neon;
    mt19937 gen_neon(rd_neon());

    double decay_rate_neon = 1.0/ mean;

    exponential_distribution<double> dist(decay_rate_neon);
 
    return dist(gen_neon);



}

int main(){

ifstream ifile("/home/simulation/Schreibtisch/Stream/CF_ArrTime/Na22_1275_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.stream", ios::binary);
ifstream ifile511("/home/simulation/Schreibtisch/Stream/CF_ArrTime/Na22_511_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.stream", ios::binary);

fstream myfile("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/300BS_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile1("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/300IRF_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile2("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/300LT_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);

fstream myfile3("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/350BS_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile4("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/350IRF_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile5("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/350LT_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);

fstream myfile6("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/400BS_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile7("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/400IRF_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile8("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/400LT_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);

fstream myfile9("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/450BS_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile10("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/450IRF_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile11("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/450LT_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);

fstream myfile12("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/500BS_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile13("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/500IRF_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);
fstream myfile14("/home/simulation/Schreibtisch/StreamReader/Data/SpectraDecomposition/500LT_130_0.625_240_0.2_Na22_ohneEindrintiefe_Cone_D4cm_Z2.79cm_D1.9cm.txt",ios::out);




int eventIndex = 0;
int counterPileUpInBackscatering = 0; 
int counterPileUpInBackscatering_1 = 0; 
int counterPileUpInBackscatering_2 = 0; 
int counterPileUpInBackscatering_3 = 0; 
int counterPileUpInBackscatering_4 = 0; 
int counterPileUpInBackscatering_5 = 0; 

int rejectCnt = 0;


double RealLifetime = 0.0;

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


//generate histogramm vector for the IRF, lifetime and backscattering
const int numBins = 10000;
const int binWidth = 5;
const int maxValue = 50000;

// Initialisieren Sie den Histogramm-Vektor mit 0 für alle Bins
std::vector<int> histogram_IRF(numBins, 0);
std::vector<int> histogram_Lifetime(numBins, 0);
std::vector<int> histogram_Backscattering(numBins, 0);

std::vector<int> _350histogram_IRF(numBins, 0);
std::vector<int> _350histogram_Lifetime(numBins, 0);
std::vector<int> _350histogram_Backscattering(numBins, 0);

std::vector<int> _400histogram_IRF(numBins, 0);
std::vector<int> _400histogram_Lifetime(numBins, 0);
std::vector<int> _400histogram_Backscattering(numBins, 0);

std::vector<int> _450histogram_IRF(numBins, 0);
std::vector<int> _450histogram_Lifetime(numBins, 0);
std::vector<int> _450histogram_Backscattering(numBins, 0);

std::vector<int> _500histogram_IRF(numBins, 0);
std::vector<int> _500histogram_Lifetime(numBins, 0);
std::vector<int> _500histogram_Backscattering(numBins, 0);


for(int i=0;i<200;i++){
    int stream2Index = 11*i+1; // Leseindex für stream2
    int stream1Index = 0;
if(ifile.is_open() || ifile511.is_open())
{
 

    while(!ifile.eof() ){//&& !ifile511.eof()
    
    // Wenn der Index stream2Index das Ende des Streams erreicht hat, setze ihn zurück
    if (stream2Index >= 50000000) {
        stream2Index = 0;
    }
    
    
    EventInfo info;

    EventInfo info511;

    ifile.seekg(stream1Index * sizeof(EventInfo));
    ifile.read(reinterpret_cast<char*>(&info), sizeof(EventInfo));
    ++stream1Index;

    ifile.read((char*) &info, sizeof(EventInfo));

    ifile511.seekg(stream2Index * sizeof(EventInfo));
    ifile511.read(reinterpret_cast<char*>(&info511), sizeof(EventInfo));
    ++stream2Index;
    
    //ifile511.read((char*) &info511, sizeof(EventInfo));



   
  // Define the decay rates and intensities
    double decay_rate1 = 1/0.130;
    double intensity1 = 0.625;
    double decay_rate2 = 1/0.240;
    double intensity2 = 0.20;

    double decay_rate3 = 1/0.380;
    double intensity3 = 0.172;
    double decay_rate4 = 1/2.75;
    double intensity4 = 0.003;
    double num = 0.0;

//Do you want to include Lifetime in your analysis

    bool DoYouWantLifetime = true;


    if(DoYouWantLifetime == true){
    // Compute the overall intensity (must sum up to 1)
    double overall_intensity = intensity1 + intensity2 + intensity3+ intensity4;
    double tolerance = 1e-9; // set a small tolerance value
    if (std::abs(overall_intensity - 1.0) > tolerance) {
        std::cerr << "Error: The intensities must sum up to 1!" << std::endl;
        return 1;
    }

    // Define the random number generator
    std::mt19937 gen(std::random_device{}());

    // Generate a random number from the combined exponential distribution
    std::exponential_distribution<double> exp_dist1(decay_rate1);
    std::exponential_distribution<double> exp_dist2(decay_rate2);
    std::exponential_distribution<double> exp_dist3(decay_rate3);
    std::exponential_distribution<double> exp_dist4(decay_rate4);
    
        double r = ((double)rand() /RAND_MAX);
        if (r <= intensity1 / overall_intensity) {            
            num = exp_dist1(gen);
        } else if (r <= (intensity1 + intensity2) / overall_intensity) {            
            num = exp_dist2(gen);
        } else if (r <= (intensity1 + intensity2+ intensity3) / overall_intensity){            
            num = exp_dist3(gen);
        } else {            
            num = exp_dist4(gen);
        }
    }
    
    // Print the result
    //std::cout << "Random number from the combined exponential distribution: " << num << std::endl;

    int uCi = 0;
    double DecayTime = 0;
    if(uCi!=0){
        double activity = 37000*uCi;

        // Calculate the decay constant based on the activity level
        double decay_constant = 1/activity;

        // Create a random number generator
        std::random_device rd1;
        std::mt19937 gen2(rd1());
        std::exponential_distribution<double> distribution(1 / decay_constant);
        double DecayTime = (distribution(gen2))*1e9;
        //std::cout << "Random number from the DecayTime: " << DecayTime << std::endl;
        
        }
    


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
            

                if (((info.info1().isValid() ||
                 info.info2().isValid()) || (info511.info1().isValid() ||
                 info511.info2().isValid())))//&&(DecayTime> 200) 
                {  
              

                   EventsInBothDetektors ++;
                   // const bool isInEnergy1Window = ((Energy1>0.17)&&(Energy1<0.4))&&((Energy2>0.56)&&(Energy2<1.1));
                    
                    // const bool isInEnergy2Window = ((Energy2>0.17)&&(Energy2<0.4))&&((Energy1>0.56)&&(Energy1<1.1));

                    const bool isInEnergy1Window = ((Energy1>=0.17)&&(Energy1<0.4))&&((Energy2>=0.56)&&(Energy2<1.1));
                    
                    const bool isInEnergy2Window = ((Energy2>=0.17)&&(Energy2<0.4))&&((Energy1>=0.56)&&(Energy1<1.1));

                    const bool isInEnergyWindow = (isInEnergy1Window || isInEnergy2Window);

                    
                    
                    
                    const bool isInEnergy1Window_511 = ((Energy1_511>=0.17)&&(Energy1_511<1.35))&&((Energy2_511>=0.17)&&(Energy2_511<1.35));
                    
                    const bool isInEnergy2Window_511 = ((Energy2_511>=0.17)&&(Energy2_511<1.35))&&((Energy1_511>=0.17)&&(Energy1_511<1.35));

                    const bool isInEnergyWindow_511 = (isInEnergy1Window_511 || isInEnergy2Window_511);

                    const double Position1X = info.info1().InterActionX();
                    const double Position1Y = info.info1().InterActionY();
                    const double Position1Z = info.info1().InterActionZ();

                    const double Position2X = info.info2().InterActionX();
                    const double Position2Y = info.info2().InterActionY();
                    const double Position2Z = info.info2().InterActionZ();

                    const double InteractionLength1 = sqrt((Position1X*Position1X) +(Position1Y*Position1Y) +(Position1Z*Position1Z));
                    const double InteractionLength2 = sqrt((Position2X*Position2X) +(Position2Y*Position2Y) +(Position2Z*Position2Z));

                    const double Position1X511 = info511.info1().InterActionX();
                    const double Position1Y511 = info511.info1().InterActionY();
                    const double Position1Z511 = info511.info1().InterActionZ();

                    const double Position2X511 = info511.info2().InterActionX();
                    const double Position2Y511 = info511.info2().InterActionY();
                    const double Position2Z511 = info511.info2().InterActionZ();

                    const double InteractionLength1511 = sqrt((Position1X511*Position1X511) +(Position1Y511*Position1Y511) +(Position1Z511*Position1Z511));
                    const double InteractionLength2511 = sqrt((Position2X511*Position2X511) +(Position2Y511*Position2Y511) +(Position2Z511*Position2Z511));

                    int Hits1 = info.info1().numberOfCounts();

                    int Hits2 = info.info2().numberOfCounts();

                    int Hits1_511 = info511.info1().numberOfCounts();

                    int Hits2_511 = info511.info2().numberOfCounts();

                    int window_511_low = 300;
                    int window_511_low1 = 350;
                    int window_511_low2 = 400;
                    int window_511_low3 = 450;
                    int window_511_low4 = 500;

                    int window_511_high = 950;

                    int window_1275_low = 1500;
                    int window_1275_high = 2230;

                    const bool isInPulsWindow1 =((Hits1>=window_511_low)&&(Hits1<window_511_high))&&((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInPulsWindow1_1 =((Hits1>=window_511_low1)&&(Hits1<window_511_high))&&((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInPulsWindow1_2 =((Hits1>=window_511_low2)&&(Hits1<window_511_high))&&((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInPulsWindow1_3 =((Hits1>=window_511_low3)&&(Hits1<window_511_high))&&((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInPulsWindow1_4 =((Hits1>=window_511_low4)&&(Hits1<window_511_high))&&((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    

                    const bool isInPulsWindow2 =((Hits2>=window_511_low)&&(Hits2<window_511_high))&&((Hits1>=window_1275_low)&&(Hits1<window_1275_high));
                    const bool isInPulsWindow = (isInPulsWindow1 || isInPulsWindow2);
                    const bool isIn511WindowFalsePositive = ((Hits2>=window_511_low)&&(Hits2<window_511_high))||((Hits1>=window_511_low)&&(Hits1<window_511_high));

                    const bool isInPulsWindow1_511 =((Hits1_511>=window_511_low)&&(Hits1_511<window_511_high))&&((Hits2_511>=window_1275_low)&&(Hits2_511<window_1275_high));
                    const bool isInPulsWindow2_511 =((Hits2_511>=window_511_low)&&(Hits2_511<window_511_high))||((Hits1_511>=window_1275_low)&&(Hits1_511<window_1275_high));
                    const bool isInPulsWindow_511 = (isInPulsWindow1_511 || isInPulsWindow2_511);


                    const bool isInRightWindow_1 = ((Hits1_511>=window_511_low)&&(Hits1_511<window_511_high))||((Hits1>=window_1275_low)&&(Hits1<window_1275_high));
                    
                    
                    const bool isInRightWindow_2 = ((Hits1_511>=window_511_low)&&(Hits1_511<window_511_high)) && ((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInRightWindow_2_1 = ((Hits1_511>=window_511_low1)&&(Hits1_511<window_511_high)) && ((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInRightWindow_2_2 = ((Hits1_511>=window_511_low2)&&(Hits1_511<window_511_high)) && ((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInRightWindow_2_3 = ((Hits1_511>=window_511_low3)&&(Hits1_511<window_511_high)) && ((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInRightWindow_2_4 = ((Hits1_511>=window_511_low4)&&(Hits1_511<window_511_high)) && ((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    
                    
                    const bool isInRightWindow_3 = ((Hits2_511>=window_511_low)&&(Hits2_511<window_511_high)) && ((Hits1>=window_1275_low)&&(Hits1<window_1275_high));
                    const bool isInRightWindow_4 = ((Hits2_511>=window_511_low)&&(Hits2_511<window_511_high))||((Hits2>=window_1275_low)&&(Hits2<window_1275_high));

                    const bool isInRightWindow_2_reverse = ((Hits2>=window_511_low)&&(Hits2<window_511_high)) && ((Hits1_511>=window_1275_low)&&(Hits1_511<window_1275_high));
                    const bool isInRightWindow_3_reverse = ((Hits1>=window_511_low)&&(Hits1<window_511_high)) && ((Hits2_511>=window_1275_low)&&(Hits2_511<window_1275_high));

                    const bool isInRightWindow511 = ((Hits1_511>=window_511_low)&&(Hits1_511<window_511_high));

                    const bool isInRightWindow511_2 = ((Hits1>=window_511_low)&&(Hits1<window_511_high));

                    const bool isInRightWindow511_det2 = ((Hits2_511>=window_511_low)&&(Hits2_511<window_511_high));
                    const bool isInRightWindow511_2_det2 = ((Hits2>=window_511_low)&&(Hits2<window_511_high));


                    const bool isInRightWindow1275 = ((Hits2>=window_1275_low)&&(Hits2<window_1275_high));
                    const bool isInRightWindow1275_2 = ((Hits2_511>=window_1275_low)&&(Hits2_511<window_1275_high));

                    const bool isInRightWindow1275_det1 = ((Hits1>=window_1275_low)&&(Hits1<window_1275_high));
                    const bool isInRightWindow1275_2_det1 = ((Hits1_511>=window_1275_low)&&(Hits1_511<window_1275_high));



                    const bool isInReverseWindow = ((Hits1>=window_511_low)&&(Hits1<window_511_high)) && ((Hits2_511>=window_1275_low)&&(Hits2_511<window_1275_high));


                    const bool isInRightWindow = (isInRightWindow_1 || isInRightWindow_2 || isInRightWindow_3 || isInRightWindow_1);

                    const bool TrueCoincidence = (isInRightWindow_2 || isInRightWindow_3 );

                    
                 
                    
                    const bool isInTimeWindow = ((info.info1().MeanArr()>0.)&&(info.info2().MeanArr()>0.)); // ns

                    const bool isInTimeWindow_511 = ((info511.info1().MeanArr()>0.)&&(info511.info2().MeanArr()>0.)); // ns

                 


                    const double IRF_1274_1 = (info.info1().MeanArr()-GlTime1274);
                   
                    const double IRF_1274_2 = (info.info2().MeanArr()-GlTime1274);
                   
                    

                    const double IRF_511_1 = (info511.info1().MeanArr()-GlTime511);
                   
                    
                    const double IRF_511_2 = (info511.info2().MeanArr()-GlTime511);
                    

                   
                    const bool TrueLifetime_1 = (info511.info1().isValid() && info.info1().isValid());
                    const bool TrueLifetime_2 = (info511.info2().isValid() && info.info1().isValid());
                    const bool TrueLifetime_3 = (info511.info1().isValid() && info.info2().isValid());
                    const bool TrueLifetime_4 = (info511.info2().isValid() && info.info2().isValid());
                    
                    


                        double QE = 0.15;

                        
                            ///////////// Backscatering in right PHS 1275keV
                        if((info.info1().isValid() && info.info2().isValid() && isInPulsWindow1)){
                            double neon_lifetime = NeonDecay(0.0037);

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1*QE)));
                            double Pmt_time1 = IRF_1274_1 + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2 + neon_lifetime + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2) + 7.750;

                            double convertValueBackScattering = IRF / 0.005;

                            int value = round(convertValueBackScattering);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                histogram_Backscattering[binIndex]++;
                                }
                            


                            if(info511.info1().isValid() || info511.info2().isValid()){
                                counterPileUpInBackscatering ++;

                            }


                        }
                        
                        

                        if(TrueLifetime_3  && isInRightWindow_2){// isInRightWindow1275 isInReverseWindow || || TrueLifetime_2 
                            //define the internal lifetime of your radioactive source for Na22 t1/2 for the excited Neon is 3.7ps
                            double neon_lifetime1 = NeonDecay(0.0037);
                            

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1_511*QE)));
                            double Pmt_time1 = IRF_511_1  + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2+ neon_lifetime1  + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2)+ 7.750;
                            RealLifetime = IRF + num;

                            double convertValueIRF = IRF / 0.005;

                            int value = round(convertValueIRF);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                histogram_IRF[binIndex]++;
                                }
                            
                            
                            double convertValueLifetime = RealLifetime / 0.005;

                            int value2 = round(convertValueLifetime);

                            
                            if (value2 >= 0 && value2 <= maxValue) {
                                int binIndex2 = value2;
                                histogram_Lifetime[binIndex2]++;
                                }
                            

                            
          
                        }

                        if((info.info1().isValid() && info.info2().isValid() && isInPulsWindow1_1)){
                            double neon_lifetime = NeonDecay(0.0037);

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1*QE)));
                            double Pmt_time1 = IRF_1274_1 + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2 + neon_lifetime + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2) + 7.750;

                            double convertValueBackScattering = IRF / 0.005;

                            int value = round(convertValueBackScattering);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _350histogram_Backscattering[binIndex]++;
                                }
                            


                            if(info511.info1().isValid() || info511.info2().isValid()){
                                counterPileUpInBackscatering ++;

                            }


                        }
                        
                        

                        if(TrueLifetime_3  && isInRightWindow_2_1){// isInRightWindow1275 isInReverseWindow || || TrueLifetime_2 
                            //define the internal lifetime of your radioactive source for Na22 t1/2 for the excited Neon is 3.7ps
                            double neon_lifetime1 = NeonDecay(0.0037);
                            

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1_511*QE)));
                            double Pmt_time1 = IRF_511_1  + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2+ neon_lifetime1  + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2)+ 7.750;
                            RealLifetime = IRF + num;

                            double convertValueIRF = IRF / 0.005;

                            int value = round(convertValueIRF);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _350histogram_IRF[binIndex]++;
                                }
                            
                            
                            double convertValueLifetime = RealLifetime / 0.005;

                            int value2 = round(convertValueLifetime);

                            
                            if (value2 >= 0 && value2 <= maxValue) {
                                int binIndex2 = value2;
                                _350histogram_Lifetime[binIndex2]++;
                                }
                            

                            
          
                        }

                        if((info.info1().isValid() && info.info2().isValid() && isInPulsWindow1_2)){
                            double neon_lifetime = NeonDecay(0.0037);

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1*QE)));
                            double Pmt_time1 = IRF_1274_1 + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2 + neon_lifetime + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2) + 7.750;

                            double convertValueBackScattering = IRF / 0.005;

                            int value = round(convertValueBackScattering);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _400histogram_Backscattering[binIndex]++;
                                }
                            


                            if(info511.info1().isValid() || info511.info2().isValid()){
                                counterPileUpInBackscatering ++;

                            }


                        }
                        
                        

                        if(TrueLifetime_3  && isInRightWindow_2_2){// isInRightWindow1275 isInReverseWindow || || TrueLifetime_2 
                            //define the internal lifetime of your radioactive source for Na22 t1/2 for the excited Neon is 3.7ps
                            double neon_lifetime1 = NeonDecay(0.0037);
                            

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1_511*QE)));
                            double Pmt_time1 = IRF_511_1  + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2+ neon_lifetime1  + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2)+ 7.750;
                            RealLifetime = IRF + num;

                            double convertValueIRF = IRF / 0.005;

                            int value = round(convertValueIRF);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _400histogram_IRF[binIndex]++;
                                }
                            
                            
                            double convertValueLifetime = RealLifetime / 0.005;

                            int value2 = round(convertValueLifetime);

                            
                            if (value2 >= 0 && value2 <= maxValue) {
                                int binIndex2 = value2;
                                _400histogram_Lifetime[binIndex2]++;
                                }
                            

                            
          
                        }

                        if((info.info1().isValid() && info.info2().isValid() && isInPulsWindow1_3)){
                            double neon_lifetime = NeonDecay(0.0037);

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1*QE)));
                            double Pmt_time1 = IRF_1274_1 + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2 + neon_lifetime + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2) + 7.750;

                            double convertValueBackScattering = IRF / 0.005;

                            int value = round(convertValueBackScattering);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _450histogram_Backscattering[binIndex]++;
                                }
                            


                            if(info511.info1().isValid() || info511.info2().isValid()){
                                counterPileUpInBackscatering ++;

                            }


                        }
                        
                        

                        if(TrueLifetime_3  && isInRightWindow_2_3){// isInRightWindow1275 isInReverseWindow || || TrueLifetime_2 
                            //define the internal lifetime of your radioactive source for Na22 t1/2 for the excited Neon is 3.7ps
                            double neon_lifetime1 = NeonDecay(0.0037);
                            

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1_511*QE)));
                            double Pmt_time1 = IRF_511_1  + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2+ neon_lifetime1  + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2)+ 7.750;
                            RealLifetime = IRF + num;

                            double convertValueIRF = IRF / 0.005;

                            int value = round(convertValueIRF);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _450histogram_IRF[binIndex]++;
                                }
                            
                            
                            double convertValueLifetime = RealLifetime / 0.005;

                            int value2 = round(convertValueLifetime);

                            
                            if (value2 >= 0 && value2 <= maxValue) {
                                int binIndex2 = value2;
                                _450histogram_Lifetime[binIndex2]++;
                                }
                            

                            
          
                        }

                        if((info.info1().isValid() && info.info2().isValid() && isInPulsWindow1_4)){
                            double neon_lifetime = NeonDecay(0.0037);

                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1*QE)));
                            double Pmt_time1 = IRF_1274_1 + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2 + neon_lifetime + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2) + 7.750;

                            double convertValueBackScattering = IRF / 0.005;

                            int value = round(convertValueBackScattering);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _500histogram_Backscattering[binIndex]++;
                                }
                            


                            if(info511.info1().isValid() || info511.info2().isValid()){
                                counterPileUpInBackscatering ++;

                            }


                        }
                        
                        

                        if(TrueLifetime_3  && isInRightWindow_2_4){// isInRightWindow1275 isInReverseWindow || || TrueLifetime_2 
                            //define the internal lifetime of your radioactive source for Na22 t1/2 for the excited Neon is 3.7ps
                          
                            double neon_lifetime1 = NeonDecay(0.0037);
                            


                            double Spread_Pmt1 = deformed_gauss_pmt1(28.00, 0.55* 1/ (sqrt(Hits1_511*QE)));
                            double Pmt_time1 = IRF_511_1  + Spread_Pmt1;

                            double Spread_Pmt2 = deformed_gauss_pmt2(28.00, 0.55* 1/ (sqrt(Hits2*QE))) ; 
                            double Pmt_time2 = IRF_1274_2+ neon_lifetime1  + Spread_Pmt2;
                            
                            
                            double IRF = (Pmt_time1 - Pmt_time2)+ 7.750;
                            RealLifetime = IRF + num;

                            double convertValueIRF = IRF / 0.005;

                            int value = round(convertValueIRF);

                            
                            if (value >= 0 && value <= maxValue) {
                                int binIndex = value;
                                _500histogram_IRF[binIndex]++;
                                }
                            
                            
                            double convertValueLifetime = RealLifetime / 0.005;

                            int value2 = round(convertValueLifetime);

                            
                            if (value2 >= 0 && value2 <= maxValue) {
                                int binIndex2 = value2;
                                _500histogram_Lifetime[binIndex2]++;
                                }
                            

                            
          
                        }
                        




                        
                       
                    
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



   
    

        
}
cout<<i<<endl;
ifile.clear();

ifile511.clear();


}

for (int i = 0; i < numBins; i++) {
        
        myfile<<histogram_Backscattering[i]<<endl;
        myfile1<<histogram_IRF[i]<<endl;
        myfile2<<histogram_Lifetime[i]<<endl;

        myfile3<<_350histogram_Backscattering[i]<<endl;
        myfile4<<_350histogram_IRF[i]<<endl;
        myfile5<<_350histogram_Lifetime[i]<<endl;

        myfile6<<_400histogram_Backscattering[i]<<endl;
        myfile7<<_400histogram_IRF[i]<<endl;
        myfile8<<_400histogram_Lifetime[i]<<endl;

        myfile9<<_450histogram_Backscattering[i]<<endl;
        myfile10<<_450histogram_IRF[i]<<endl;
        myfile11<<_450histogram_Lifetime[i]<<endl;

        myfile12<<_500histogram_Backscattering[i]<<endl;
        myfile13<<_500histogram_IRF[i]<<endl;
        myfile14<<_500histogram_Lifetime[i]<<endl;

    }



    std::cout<<"Number of Events :"<<eventIndex<<endl;
    std::cout<<"Number of Events in Spectrum :"<<ValidEvents<<endl;
   

    

    ifile.close();
    ifile511.close();
    myfile.close();
    myfile1.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
    myfile5.close();
    myfile6.close();
    myfile7.close();
    myfile8.close();
    myfile9.close();
    myfile10.close();
    myfile11.close();
    myfile12.close();
    myfile13.close();
    myfile14.close();


    
    

  
    return 0;

}
