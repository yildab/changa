/*
DM-induced baryonic feedback module.
*/

#include <math.h>
#include <sys/stat.h>
#include "ParallelGravity.h"
#include "smooth.h"
#include "dmfeedback.h"
// #include "starlifetime.h"
// #include "supernova.h"
// #include "romberg.h"

///
///
void Main::DMFeedbackLocalCalculations(double dTime,double dDelta, int activeRung) { 
#ifdef DMFEEDBACK
        double stime = CkWallTimer();
        CkPrintf("DM-induced feedback: Calculate local DM velocity and density for each star ... ");
        DMFeedbackSmoothParams pDMFeedback(TYPE_DARK, activeRung, param.csm, dTime, param.dDMSigma, param.dDMMass,param.dDelta, param.dMsolUnit, param.dKpcUnit, param.dSecUnit );  
        treeProxy.startSmooth(&pDMFeedback, 1, param.nSmooth, 0, CkCallbackResumeThread());         
        CkPrintf("took %g seconds\n", CkWallTimer() - stime);
#endif
    }

int DMFeedbackSmoothParams::isSmoothActive(GravityParticle *p) {
    return p->isStar();
    }

void DMFeedbackSmoothParams::initSmoothParticle(GravityParticle *p) {
    }

void DMFeedbackSmoothParams::initTreeParticle(GravityParticle *p) {
   }

void DMFeedbackSmoothParams::initSmoothCache(GravityParticle *p) {
    //p->treeAcceleration = p->velocity;
    }

void DMFeedbackSmoothParams::combSmoothCache(GravityParticle *p1, ExternalSmoothParticle *p2) {
    //Vector3D<double> deltav; //vector template notation
    //deltav=p2->velocity - p2->treeAcceleration;
    //p1->velocity += deltav;
    }

void DMFeedbackSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nnList) {\
    GravityParticle *q;
    double fNorm,ih2,r2;
    double aDot, dvx, dvy, dvz, dv, dv_avg;
    double rs, ddensity;
    int i;

    aDot = a*H;
    ih2 = invH2(p);
    fNorm = M_1_PI*ih2*sqrt(ih2); 
    dv_avg = 0.0;
    ddensity = 0.0;
    for (i = 0; i < nSmooth; ++i)
    {
        q = nnList[i].p;
        if (q->iOrder != p->iOrder)  // particle should not interact with itself
        {
            // we want the relative velocity between star particle p and DM particle q
            // physical distance (how many km/s?)
                dvx = (q->velocity.x - p->velocity.x)/a - aDot*nnList[i].dx.x;
                dvy = (q->velocity.y - p->velocity.y)/a - aDot*nnList[i].dx.y;
                dvz = (q->velocity.z - p->velocity.z)/a - aDot*nnList[i].dx.z;

                dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
                dv_avg += dv;

            // we also want local DM density (method taken from DensitySmoothParams)
		        double fDist2 = nnList[i].fKey;
		        r2 = fDist2*ih2;

                rs = KERNEL(r2, nSmooth);
                ddensity += rs*nnList[i].p->mass; // units (dMsolUnit/dKpcUnit**3)
        }
    }
    double vconversion = 65.6; // got this from output files, need to find definition in code
    p->dDMVel() = dv_avg*vconversion/nSmooth; // km/s
    double dDMVelC = dv_avg*vconversion/nSmooth/(3e5); //natural units
    double KpctoPc = 0.001;
    double GeVCm3 = 37.96;
    double dDMDensityGeVcm3 = fNorm*ddensity * (dMsolUnit/(dKpcUnit * dKpcUnit * dKpcUnit)) * (KpctoPc*KpctoPc*KpctoPc) * GeVCm3; //
    p->dDMDensity() = dDMDensityGeVcm3;
    double ddMcapGeVSec = (6.0e24)*(1e8/pow(10,dDMMass))*(pow(10, dDMSigma)/(1e-42))*(dDMDensityGeVcm3/0.3)*(1e-3/dDMVelC);
    p->ddMcap() = min(3e27, ddMcapGeVSec);
    }

double Mcrit(double dWDMass, double dWDage, double dDMMass){
    // WD core temperature, as a function of WD age to incorporate WD cooling
    double dTWD = (2.3e10)*pow(dWDage, -2.5); 

    // WD central density in g/cm^3
    double mu_e = 0.5;
    double ge = 2;
    double mn = 0.93;
    double me = 0.511e-3;
    double B = (ge*pow(me, 3)*mn*mu_e)/(6*M_PI*M_PI);
    double alphaWD = 1.0033-0.3087*dWDMass - 1.1652*pow(dWDMass, 2) +2.0211*pow(dWDMass, 3) - 2.0604*pow(dWDMass, 4) +01.1687*pow(dWDMass, 5) -0.2810*pow(dWDMass, 6);
    double rhoWD = pow(pow(alphaWD, -2) - 1, 3/2)*(B/(pow(1.98e-14, 3)*5.62e23));

    // thermalization radius
    double M_pl = 1.22e19;
    double G = pow(M_pl, -2);
    double rth = (1.98e-14)*pow((9*dTWD*(8.62e-14))/(4*M_PI*G*(5.62e23)*(pow(1.98e-14, 3))*dDMMass), 1/2);

    double Mcrit = 4*M_PI*rhoWD*(5.62e23)*pow(rth, 3)/3;

    return Mcrit;
}

double CheckExplodes(double dWDMass, double dTime, double dDelta, SFEvent *sfEvent, double ddMcap, Padova pdva, double dDMMass){
    double progenitorMass = (dWDMass - 0.565)/0.0835; // Cummings et Al Parsec model
    double progenitorLifetime = pdva.Lifetime(progenitorMass, sfEvent->dMetals);
    double WDage = (dTime - sfEvent->dTimeForm) - progenitorLifetime;

    if (WDage < 0) {
        return false;
    }

    double dMcrit = Mcrit(dWDMass, WDage, dDMMass);
    double Mcap = ddMcap*WDage;
    double McapLast = ddMcap*(WDage - dDelta);

    if ((Mcap > dMcrit) && (McapLast < dMcrit)) {
        return true;
    }
    return false;
}

// double NumDarkSNIa(double dWDMass, SFEvent *sfEvent, imf *imf){
//     double progenitorMass = (dWDMass - 0.565)/0.0835; // Cummings et Al Parsec model
    
//     double numProgenitors = imf->CumMass(progenitorMass+0.05) - imf->CumMass(progenitorMass- 0.05);

//     return numProgenitors;
// }

/* 
void DarkSN::CalcDarkIaFeedback(SFEvent *sfEvent, DarkSNData *darkSNdata, double dTime, double dDelta, FBEffects *fbEffects) 
{   
    // loop over 
}
 */