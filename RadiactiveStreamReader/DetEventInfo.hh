#pragma once
#ifndef DETEVENTINFO_HH
#define DETEVENTINFO_HH


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
        //m_ArrCF = 0.;
        m_MeanArr = 0.;
        //m_MedianArr = 0.;
       // m_MaxPulse = 0.;
        m_globalTime = 0;
        m_nofHits = false;
        m_interactionX = 0.;
        m_interactionY =0.;
        m_interactionZ =0.;
        m_numberOfCounts = 0;
     
	}
	
	
	void AddEnergyDep(double fEdep, double globalTime, bool bHit, int id, double MeanArr, double InterActionX, double InterActionY,
    double InterActionZ, int numberOfCounts/*, double MedianArr, double MaxPulse*/) {
		m_depEnerg = fEdep;
		//m_ArrCF = ArrCF;
		m_globalTime = globalTime;
		m_nofHits = bHit;
		m_id = id;
        m_MeanArr = MeanArr;
        m_interactionX = InterActionX;
        m_interactionY = InterActionY;
        m_interactionZ = InterActionZ;
        m_numberOfCounts = numberOfCounts;
     
       // m_MedianArr = MedianArr;
       // m_MaxPulse = MaxPulse;
	}
	
  

    int id() const {
        return m_id;
    }

    double energy() const {
        return m_depEnerg;
    }
    
   /* double ArrCF() const {
        return m_ArrCF;
    }*/

    double MeanArr() const {
        return m_MeanArr;
    }
    double InterActionX() const {
        return m_interactionX;
    }
    double InterActionY() const {
        return m_interactionY;
    }
    double InterActionZ() const {
        return m_interactionZ;
    }
   int numberOfCounts() const {
       return m_numberOfCounts;
   }

/*
     double MedianArr() const {
        return m_MedianArr;
    }

    double MaxPulse() const {
        return m_MaxPulse;
    }
*/
    bool isValid() const {
        return m_nofHits;
    }

private:
    int m_id; // event-id
    double m_depEnerg; // dep. energy at detector
    //double m_ArrCF; // arrival time at detector
    double m_MeanArr; //mean arrvialValue
    double m_globalTime; //globalTime counter
    double m_MedianArr; //Median value of the pulse
    double m_MaxPulse; //Pulse Peak x value
    bool m_nofHits; // incident gammay ray?
    double m_interactionX;
    double m_interactionY;
    double m_interactionZ;
    int m_numberOfCounts;
    

};

class EventInfo {
public:
    EventInfo(){clear();}
    ~EventInfo(){}
    
    void clear() {
		m_lifetime = 0.;
		m_startTime = 0.;
        m_StartWinkelX = 0.;
        m_StartWinkelY= 0.;
        m_StartWinkelZ= 0.;
        m_StoppWinkel1X= 0.;
        m_StoppWinkel1Y= 0.;
        m_StoppWinkel1Z= 0.;
        m_StoppWinkel2X= 0.;
        m_StoppWinkel2Y= 0.;
        m_StoppWinkel2Z= 0.;
		m_detector1.clear();
		m_detector2.clear();
	}

	void setInfo(double lifetime, double startTime, float StartWinkelX, float StartWinkelY, float StartWinkelZ, float StoppWinkel1X,
    float StoppWinkel1Y, float StoppWinkel1Z, float StoppWinkel2X, float StoppWinkel2Y, float StoppWinkel2Z) {
		m_lifetime = lifetime;
		m_startTime = startTime;
        m_StartWinkelX = StartWinkelX;
        m_StartWinkelY = StartWinkelY;
        m_StartWinkelZ = StartWinkelZ;
        m_StoppWinkel1X = StoppWinkel1X;
        m_StoppWinkel1Y = StoppWinkel1Y;
        m_StoppWinkel1Z = StoppWinkel1Z;
        m_StoppWinkel2X = StoppWinkel2X;
        m_StoppWinkel2Y = StoppWinkel2Y;
        m_StoppWinkel2Z = StoppWinkel2Z;
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

    float StartWinkelX() const {
        return m_StartWinkelX;
    }

    float StartWinkelY() const {
        return m_StartWinkelY;
    }

    float StartWinkelZ() const {
        return m_StartWinkelZ;
    }

    float StoppWinkel1X() const {
        return m_StoppWinkel1X;
    }
    float StoppWinkel1Y() const {
        return m_StoppWinkel1Y;
    }
    float StoppWinkel1Z() const {
        return m_StoppWinkel1Z;
    }

     float StoppWinkel2X() const {
        return m_StoppWinkel2X;
    }
    float StoppWinkel2Y() const {
        return m_StoppWinkel2Y;
    }
    float StoppWinkel2Z() const {
        return m_StoppWinkel2Z;
    }

private:
    DetectorInfo m_detector1;
    DetectorInfo m_detector2;
	
    double m_lifetime;
    double m_startTime;
    float m_StartWinkelX;
    float m_StartWinkelY;
    float m_StartWinkelZ;
    float m_StoppWinkel1X;
    float m_StoppWinkel1Y;
    float m_StoppWinkel1Z;
    float m_StoppWinkel2X;
    float m_StoppWinkel2Y;
    float m_StoppWinkel2Z;
};

#endif
