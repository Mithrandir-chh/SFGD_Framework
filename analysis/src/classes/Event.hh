#ifndef _EVENT_H_
#define _EVENT_H_

#include <TObject.h>
#include <TClonesArray.h>
#include "Hit.hh"

class Event : public TObject {

private:

    TClonesArray*    fHits;                   //Collection of hits in the event
    Int_t            fNHits;                  //Number of hits in the event
    ULong_t          fEventID;                //Number of the event
    Double_t         fMaxCharge;              //Maximum energy in the event
    Int_t            fFEB12ch;                //FEB 12 channel
    Double_t         fFEB12LeadTime;          //lead time of FEB12 trigger
    Double_t         fFEB12hitTfromSpill;     //FEB12 hit time since spill
    Int_t            fMicropulse;             //the number of micropulses since the macropulse (spill) start
    Int_t            fOccupancyXZ[48];        //Occupancy per layer of Z in the XZ view
    Int_t            fOccupancyZY[48];        //Occupancy per layer of Z in the ZY view
    Double_t         fdEdzXZ[48];             //Energy loss along Z in XZ plane
    Double_t         fdEdzZY[48];             //Energy loss along Z in ZY plane
    Double_t         fdEdz[48];               //Energy loss along Z in both XZ and ZY planes (sum)
    Int_t            fRange;                  //last layer reached by the particle
    Int_t            fZEmpty;                 //number of Z layers with no charge deposit
    Double_t         fXStdDev;                //the standard deviation of hits along X
    Double_t         fYStdDev;                //the standard deviation of hits along Y

    Int_t            fNumberCube;
    Int_t            fNumTrueTracks;        
    Double_t	     fTrueVtx[10];
    Double_t	     fColmt[2];
    Double_t	     fTrueLocation[1000][4];
    Int_t	     fTrueVoxelNumber;
    Double_t	     fnElastic[100];

    Double_t         fPointPosition[100][4];

public:

    //constructors
    Event()
    {
      fNHits = 0;
      fHits  = new TClonesArray("Hit", 10000);
    }

    //destructor
    virtual ~Event()
    {
      Reset(0);
    }

    //functions to store data
    void SetEventID(Int_t p_fEventID) { fEventID = p_fEventID; }
    void SetMaxCharge(Int_t p_MaxCharge) {fMaxCharge = p_MaxCharge;}
    void SetOccupancyXZ(Int_t p_Z, Int_t p_N) {fOccupancyXZ[p_Z] = p_N; }
    void SetOccupancyZY(Int_t p_Z, Int_t p_N) {fOccupancyZY[p_Z] = p_N; }
    void SetdEdzXZ(Int_t p_Z, Double_t p_E) {fdEdzXZ[p_Z] = p_E; }
    void SetdEdzZY(Int_t p_Z, Double_t p_E) {fdEdzZY[p_Z] = p_E; }
    void SetdEdz(Int_t p_Z, Double_t p_E) {fdEdz[p_Z] = p_E; }
    void SetRange(Int_t p_Z) {fRange = p_Z; }
    void SetEmptyZ(Int_t p_Z) {fZEmpty = p_Z; }
    void SetXStdDev(Double_t p_Z) {fXStdDev = p_Z; }
    void SetYStdDev(Double_t p_Z) {fYStdDev = p_Z; }
    void SetFEB12ch(Double_t p_FEB12ch) {fFEB12ch = p_FEB12ch;}
    void SetFEB12LeadTime(Int_t p_FEB12LeadTime) {fFEB12LeadTime = p_FEB12LeadTime;}
    void SetFEB12hitTfromSpill(Int_t p_FEB12hitTfromSpill) { fFEB12hitTfromSpill = p_FEB12hitTfromSpill; }
    void SetMicropulse(Int_t p_Micropulse) {fMicropulse = p_Micropulse; }
    void SetNumTrueTracks(Int_t f_NTrueT) {fNumTrueTracks = f_NTrueT;}
    void SetNumberCube(Int_t f_Ncube) {fNumberCube = f_Ncube;}
    void SetTrueVtx(Int_t ii, Double_t tt) {fTrueVtx[ii] = tt;}
    void SetColmt(Int_t ii, Double_t tt){fColmt[ii] = tt;}
    void SetTrueLocation(Int_t ii, Int_t jj, Double_t kk) {fTrueLocation[ii][jj] = kk; }
    void SetTrueVoxelNumber(Int_t xx) {fTrueVoxelNumber = xx;}
    void SetnElastic(Int_t xx, Double_t yy) {fnElastic[xx] = yy;}
    void SetPointPosition(Int_t xx, Int_t yy, Double_t zz) { fPointPosition[xx][yy] = zz; }

    //functions to retrieve data
    TClonesArray *GetHits() const {return fHits;}
    Int_t GetNHits()  const { return fNHits; }
    ULong_t GetEventID() {return fEventID; }
    Double_t GetMaxCharge() {return fMaxCharge; }
    Int_t GetFEB12ch() {return fFEB12ch;}
    Double_t GetFEB12LeadTime() { return fFEB12LeadTime; }
    Double_t GetFEB12hitTfromSpill() { return fFEB12hitTfromSpill; }
    Int_t GetMicropulse() { return fMicropulse; }
    Int_t GetOccupancyXZ(Int_t p_Z) {return fOccupancyXZ[p_Z];}
    Int_t GetOccupancyZY(Int_t p_Z) {return fOccupancyZY[p_Z];}
    Double_t GetdEdzXZ(Int_t p_Z) {return fdEdzXZ[p_Z];}
    Double_t GetdEdzZY(Int_t p_Z) {return fdEdzZY[p_Z];}
    Double_t GetdEdz(Int_t p_Z) {return fdEdz[p_Z];}
    Int_t GetRange() {return fRange;}
    Int_t GetEmptyZ() {return fZEmpty;}
    Double_t GetXStdDev() {return fXStdDev;}
    Double_t GetYStdDev() {return fYStdDev;}
    Double_t GetTrueVtx(Int_t ii) {return fTrueVtx[ii];}
    Double_t GetColmt(Int_t ii) {return fColmt[ii];}
    Int_t GetTrueVoxelNumber() {return fTrueVoxelNumber;}
    Double_t GetnElastic(Int_t ii) {return fnElastic[ii];}
    Double_t GetPointPosition(Int_t ii, Int_t jj) {return fPointPosition[ii][jj];}

    Int_t GetNumTrueTracks() {return fNumTrueTracks;}
    Double_t         GetTrueLocation(Int_t xx, Int_t yy)       { return fTrueLocation[xx][yy]; }

    //methods
    void Clear(Option_t *option="");
    Hit* AddHit();
    Double_t FindMaxCharge();


    void Reset(Option_t* /*option*/="")
    {
      delete fHits;  fHits = 0;
    }

    ClassDef (Event,2);

};


#endif // _CLASSES_H_
