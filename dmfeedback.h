#include "rand.h"
#include "smooth.h"
// #include "starlifetime.h"
// #include "imf.h"


// allows us to calculate local parameters (vx, rhox) for star particles
class DMFeedbackSmoothParams : public SmoothParams
{
    double a, H; // Cosmological parameters
    double dTime;
    double dDelta;
    double dDMSigma;
    double dDMMass;
    double dMsolUnit;
    double dKpcUnit;
    double dSecUnit;
  	virtual void fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p);
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,ExternalSmoothParticle *p2);

  public:
    DMFeedbackSmoothParams() {} //empty constructor
    DMFeedbackSmoothParams(int _iType, int am, CSM csm, double _dTime, double _dDMSigma, double _dDMMass, double _dDelta, double _dMsolUnit, double _dKpcUnit, double _dSecUnit) {
    iType = _iType;
    activeRung = am;
    dTime = _dTime;
    dDelta = _dDelta;
    dDMSigma = _dDMSigma;
    dDMMass = _dDMMass;
    dMsolUnit = _dMsolUnit;
    dKpcUnit = _dKpcUnit;
    dSecUnit = _dSecUnit;
        if(csm->bComove) {
            H = csmTime2Hub(csm,dTime);
            a = csmTime2Exp(csm,dTime);
            }
        else {
            H = 0.0;
            a = 1.0;
            }

    }
    PUPable_decl(DMFeedbackSmoothParams);
    DMFeedbackSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call  base class
        p|dTime;
        p|dDelta;
        p|dDMSigma;
        p|dDMMass;
        p|dMsolUnit;
        p|dKpcUnit;
        p|a;
        p|H;
        }
    };

// Info to be passed to DARKSN to characterize Star Formation event,
// based on SFEvent class defined in feedback.h
// /* 
// class DarkSNEvent {
//   public:
//     double dDMMass; /* DM particle mass in GeV*/
//     double dDMSigma; /* DM nucleon cross section */
//     double ddMCap; /* DM capture rate in GeV per second */

//   DarkSNEvent(): dMcap(0), dNumWD(0) {};
//   DarkSNEvent(double *mcap, double *numWD) : dMcap(mcap), dNumWD(numWD) {};
// };

// double Mcrit(double dWDMass, double dWDage, double dDMMass);

// double CheckExplodes(double dWDMass, SFEvent *sfEvent, double ddMcap, Padova pdva);

// double NumDarkSNIa(double dWDMass, SFEvent *sfEvent);

// // Calculating effects of dark feedback from SNIAs. Modeled on SN 
// // class defined in supernova.h
// class DarkSN
// {
//     friend class Fdbk;
//     Padova pdva;

// public:
//     double dESN;
//     IMF *imf;

//     DARKSN() {
//         dESN = 1.e51;
//     }

//     // will calculate and propagate the dark FB effects to the code
//     void CalcDarkIaFeedback(SFEvent *sfEvent, DarkSNEvent *darksfEvent, double dTime, 
//                             double dDelta, FBEffects *fbEffects) const;
    
//     // what is energy (erg) ejected from a SNIa of sub-Chandrasekhar WD
//     double DarkSNIaE(double dWDMass) const;

//     // metals (Msun) ejected from sub-Chandrasekhar WD SNIa
//     double DarkMetalsLoss(double dWDMass) const;

//     // metal fraction of oxygen ejected from sub-Chandrasekhar WD SNia
//     double DarkMetalFracO(double dWDMass) const;

//     // metal fraction of iron ejected from sub-Chandrasekhar WD SNia
//     double DarkMetalFracFe(double dWDMass) const;

//     void pup(PUP::er& p) {
//       p|dESN;
//       p|pdva;
//       // p|*imf;  Think about this
//         };
// }; */