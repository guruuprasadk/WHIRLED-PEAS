
#pragma once
#ifndef WPCOMPOSE_H
#define WPCOMPOSE_H

#include "WPGen.hpp"



namespace mrtud{
class WPCompose {
public:
    enum SpiralDirect {
        SPIRAL_OUT = 1,
        SPIRAL_IN = 2,
        SPIRAL_INOUT = 3,
        SPIRAL_OUTIN = 4
    };

    // Constructor for same Tau spiral
    WPCompose(WPGen* wp1_, int numEchoes_, SpiralDirect spiralType_);
    // Constructor for dual Tau spiral
    WPCompose(WPGen* wp1_, WPGen* wp2_, int centrNullPts_);

    // Destructor
    ~WPCompose();

    void composeGradients(dVector& gradX, dVector& gradY);

    void composeGradientsEPIC(int* gradX, int* gradY, int max_pg_integer);

    void composeKSPnSDC(dVector& kspX, dVector& kspY, dVector& sdc);

    void composeTimeMap(std::vector<std::vector<double> >& timeMap); 

    int getNEchoes() const { return mNEchoes; }

    int getKspPts() const { return mKspPts; }

    int getGradPts() const { return mGradPts; }

    int getGridMtxSize() const { return mWPOut->GetGridMtxSize(); }

    double getTimeEcho1() const { return mTimeEcho1; }
    double getDeltaTE() const { return mDeltaTE; }

    double getTotalDur() const { return mTotalDur; }

    int getTotalPts() const { return mTotalPts; }

    double getGradMax () const { return mGradMax; }

    double getSlewMax () const { return mSlewMax; }

    double getTrailingDur () const{ return mTrailingDur;}

    double getReadOutDur () const{ return readOutDur; }

    double getLeadingDur () const{ return mLeadingDur;}

    double getEchoToEndDur () const{ return mEchoToEndDur;}

    int getNyqArms() const { return mNyqArms; }

private:
    WPGen* mWPIn;
    WPGen* mWPOut;

    int mCentrNullPts;

    int mNEchoes;
    

    SpiralDirect mSpDirection;

    dVector mSpiralOutX, mSpiralOutY;
    dVector mRampDownX, mRampDownY;
    dVector mTrapOutX, mTrapOutY;
    dVector mKspOutX, mKspOutY, mSdcOut;
    int mKspPts, mGradPts;

    double mTimeEcho1;
    double mDeltaTE;

    double mTotalDur;

    double readOutDur;

    int mTotalPts;

    int mNyqArms;

    double mGradMax, mSlewMax;

    double mTrailingDur; // Time from end of readout to end of waveform
    double mLeadingDur; // Time from start of waveform to the start of readout
    double mEchoToEndDur; //Time from echo to the end of waveform

    double mGradRast;
    double mDwell;


    dVector mKspInX, mKspInY, mSdcIn;

    dVector     mSpiralInX, mSpiralInY, mRampUpX, mRampUpY, mTrapInX, mTrapInY;

    dVector _concatenate(const dVector& input1, const dVector& input2);

    dVector _concatenate(const dVector& input1, const dVector& input2, const dVector& input3);

    dVector _concatenate(const dVector& input1, const dVector& input2, const dVector& input3, const dVector& input4);

    dVector _repeat(const dVector& input, int nTimes);

    dVector _flip(const dVector& input);

    dVector _negative(const dVector& input);

    double _roundUpGrad(double dur);

    void _phasorConjucate(const dVector& x, const dVector& y, dVector& phCOnjX, dVector& phCOnjY);

    dVector _linearComb(const dVector& vect1, double a, const char op, const dVector& vect2, double b);
};
}
#endif // WPCOMPOSE_H
