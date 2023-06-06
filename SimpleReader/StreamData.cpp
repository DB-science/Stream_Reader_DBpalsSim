#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>


using namespace std;

class DetectorInfo { // A || B
    friend class EventInfo;
public:
    DetectorInfo() {
        clear();      
    }

    DetectorInfo(int id) {
		clear(); 
		
        m_id = id;
    }

    ~DetectorInfo(){}
    
    void clear() {
		m_id = -1;
        m_depEnerg = 0;
        m_meanArrTime = 0.;
        m_globalTime = 0;
        m_nofHits = false;
	}
	
	
	void AddEnergyDep(double fEdep, double globalTime, double meanArrTime, bool bHit, int id) {
		m_depEnerg = fEdep;
		m_meanArrTime = meanArrTime;
		m_globalTime = globalTime;
		m_nofHits = bHit;
		m_id = id;
	}
	
    int id() const {
        return m_id;
    }

    double energy() const {
        return m_depEnerg;
    }
    
    double meanArrivalTime() const {
        return m_meanArrTime;
    }

    bool isValid() const {
        return m_nofHits;
    }

private:
    int m_id; // event-id
    double m_depEnerg; // dep. energy at detector
    double m_meanArrTime; // arrival time at detector

    double m_globalTime; //globalTime counter

    bool m_nofHits; // incident gammay ray?
};


class EventInfo {
public:
    EventInfo(){
        
        
        clear();
    }

    ~EventInfo(){
        
    }
    
    void clear() {
		m_lifetime = 0.;
		m_startTime = 0.;
		m_detector1.clear();
		m_detector2.clear();
	}

	void setInfo(double lifetime, double startTime) {
		m_lifetime = lifetime;
		m_startTime = startTime;
	}
	
    void attach(const DetectorInfo& detector1, const DetectorInfo& detector2) {
        m_detector1 = detector1;
        m_detector2 = detector2;
    }

    DetectorInfo info1() const {
        return m_detector1;
    }

    DetectorInfo info2() const {
        return m_detector2;
    }

    double lifetime() const {
        return m_lifetime;
    }

    double startTime() const {
        return m_startTime;
    }

private:
    DetectorInfo m_detector1;
    DetectorInfo m_detector2;
	
    double m_lifetime;
    double m_startTime;
};

int main(){
ifstream ifile("/home/simulation/Schreibtisch/Stream/Stream.stream", ios::binary);
//fstream myfile("/home/simulation/Schreibtisch/StreamReader/ReferenzGlobalTime/Time1.txt",ios::out);
fstream myfile1("/home/simulation/Schreibtisch/StreamReader/klammer.txt",ios::out);
fstream myfile2("/home/simulation/Schreibtisch/StreamReader/_klammer.txt",ios::out);
//fstream myfile3("/home/simulation/Schreibtisch/StreamReader/Referenz/Energy21.txt",ios::out);

// Vector for theoLifetime

const float bucket_size = 0.005;
int number_of_buckets = (int)ceil(50/ bucket_size); // requires <cmath>
std::vector<int> histogram_theolifetime(number_of_buckets);

// Vector for Lifetime

const float bucket_size1 = 0.005;
int number_of_buckets1 = (int)ceil(50/ bucket_size1); // requires <cmath>
std::vector<int> histogram_lifetime(number_of_buckets1);

// Energy1 and Energy2

const float bucket_size2 = 0.001;
int number_of_buckets2 = (int)ceil(2/ bucket_size2); // requires <cmath>
std::vector<int> histogram_Energy1(number_of_buckets2);


const float bucket_size3 = 0.001;
int number_of_buckets3 = (int)ceil(2/ bucket_size3);
std::vector<int> histogram_Energy2(number_of_buckets2);

if(ifile.is_open())
{
    
    // streamed file operations
    ifile.seekg(0, std::ios_base::end);
    auto length = ifile.tellg();
    ifile.seekg(0, std::ios_base::beg);

    std::vector<EventInfo> input;
    
    input.resize(static_cast<size_t>(length)); // for schleife 
    ifile.read(reinterpret_cast<char*>(input.data()), length);
    
    auto success = !ifile.fail() && length == ifile.gcount();
    
    int eventIndex = 0;
    int rejectCnt = 0;
    
               
    for (EventInfo info : input) {
        // reader reinbauen
        const double Energy1 = info.info1().energy();
        const double Energy2 = info.info2().energy();
            
            if(myfile1.is_open())
            {
                
             /* int bucket1 = (int)floor(Energy1 / bucket_size2);
                    histogram_Energy1[bucket1] += 1;
              if (Energy2<1000){
                        int bucket2 = (int)floor(Energy2 / bucket_size3);
                        histogram_Energy2[bucket2] += 1;
                        }*/

               
                        
  
                if (info.info1().isValid()
                && info.info2().isValid()) 
                {
                    const bool isInEnergy1Window = ((Energy1>0.18)&&(Energy1<0.39))||((Energy1>0.66)&&(Energy1<1.1));
                    
                    const bool isInEnergy2Window = ((Energy2>0.18)&&(Energy2<0.39))||((Energy2>0.66)&&(Energy2<1.1));

                    const bool isInEnergyWindow = (isInEnergy1Window && isInEnergy2Window);
                    const double Lifetime = info.info1().meanArrivalTime() - info.info2().meanArrivalTime();

                    const bool isInTimeWindow = abs(Lifetime) < 200.; // ns

                    const double  TheoLifetime = info.lifetime();
                    const double AbsLifetime = abs(Lifetime);

                 /*if (TheoLifetime<50){
                         int bucket = (int)floor(TheoLifetime / bucket_size);
                        histogram_theolifetime[bucket] += 1;
                         }*/

                        if (isInEnergyWindow
                            && isInTimeWindow) 
                        {
                         // add to file 
                            if (AbsLifetime<50){
                            int bucket3 = (int)floor(AbsLifetime / bucket_size1);
                            histogram_lifetime[bucket3] += 1;
                            }

                            myfile2<<AbsLifetime<<endl;
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
    eventIndex ++;
    }
    std::vector<int>::iterator itr;
    /*for (itr=histogram_theolifetime.begin();itr!=histogram_theolifetime.end();itr++)
    {
        myfile<<*itr<<endl;

    }*/
    for (itr=histogram_lifetime.begin();itr!=histogram_lifetime.end();itr++)
    {
        myfile1<<*itr<<endl;

    }
   /* for (itr=histogram_Energy1.begin();itr!=histogram_Energy1.end();itr++)
    {
        myfile2<<*itr<<endl;

    }
    for (itr=histogram_Energy2.begin();itr!=histogram_Energy2.end();itr++)
    {
        myfile3<<*itr<<endl;

    }*/

    const int countsInSpectrum = input.size() -  rejectCnt;
    std::cout<<"Number of Counts :"<<countsInSpectrum<<endl;
    ifile.close();
    //myfile.close();
    myfile1.close();
    myfile2.close();
    //myfile3.close();
}
  
    return 0;

}
