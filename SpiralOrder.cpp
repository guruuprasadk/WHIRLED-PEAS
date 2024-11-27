#include "SpiralOrder.hpp"
namespace mrtud{
SpiralArmOrder::SpiralArmOrder(int numberArms, int numberEcho, float angularCoverage, OrderType orderType, 
bool inCoherentMech, bool debug)
    : mNArms(numberArms),
      mNEcho(numberEcho),
      mAngCoverage(angularCoverage),
      mOrderType(orderType),
      mInCohMech(inCoherentMech),
      mDebug(debug)
{
    if (mAngCoverage < 0 || mAngCoverage > 2 * M_PI)
    {
        mAngCoverage = 2 * M_PI;
    }



    mTotalAngles = mNArms;

    if(inCoherentMech)
        mTotalAngles *= mNEcho;

    armAngles.resize(mTotalAngles);
    armOrder.resize(mTotalAngles);
}

SpiralArmOrder::~SpiralArmOrder(){}

void SpiralArmOrder::_linearOrder(iVector &orderArray)
{
    for (int i = 0; i < mTotalAngles; i++)
    {
        orderArray[i] = i;
        if (mDebug){
            printf("Order: %d, index: %d\n", orderArray[i], i);
        }
    }
}

void SpiralArmOrder::_skipOrder(iVector &orderArray)
{
    int currOrder = 0;
    for (int i = 0; i < mTotalAngles; i++)
    {
        orderArray[i] = currOrder;

        currOrder += 2;
        if (currOrder >= mTotalAngles)
        {
            currOrder = 1;
        }

        if (mDebug){
            printf("Order: %d, index: %d\n", orderArray[i], i);
        }
    }
}

void SpiralArmOrder::_twowayOrder(iVector &orderArray)
{
    int currOrder = 0;
    int evenCounter = 0;
    int oddCounter = mTotalAngles - 1;
    for (int i = 0; i < mTotalAngles; i++)
    {
        if (i % 2)
        {
            currOrder = oddCounter;
            oddCounter--;
        }
        else
        {
            currOrder = evenCounter;
            evenCounter++;
        }
        orderArray[i] = currOrder;
        if (mDebug){
            printf("Order: %d, index: %d\n", orderArray[i], i);
        }
    }
}

void SpiralArmOrder::_mixedOrder(iVector &orderArray)
{
    int currOrder = 0;
    int evenCounterA = 1;
    int oddCounterA = mTotalAngles - 2;
    int evenCounterB = mTotalAngles - 1;
    int oddCounterB = 0;
    for (int i = 0; i < mTotalAngles; i++)
    {
        if (i < mTotalAngles / 2)
        {
            currOrder = i % 2 ? oddCounterA : evenCounterA;
            oddCounterA -= 2;
            evenCounterA += 2;
        }
        else
        {
            currOrder = i % 2 ? oddCounterB : evenCounterB;
            oddCounterB += 2;
            evenCounterB -= 2;
        }

        orderArray[i] = currOrder;
        if (mDebug){
            printf("Order: %d, index: %d\n", currOrder, i);
        }
    }
}

dVector SpiralArmOrder::_fillAngles(const iVector &orderArray)
{
    dVector angleArray(mTotalAngles, 0.0);

    double angleInc =  mAngCoverage / static_cast<double>(mTotalAngles);
    if (mOrderType == OrderType::GOLDEN)    
        angleInc = GOLDENAGNLE * M_PI / 180.0;

    for (int i = 0; i < mTotalAngles; i++)
    {
        angleArray[i] = fmod(angleInc * orderArray[i], 2 * M_PI);
    }

    return angleArray;
}

void SpiralArmOrder::composeAngles(std::vector<dVector> &spiralArmAngles)
{
    iVector orderArray(mTotalAngles);

    switch (mOrderType)
    {
    case OrderType::LINEAR:
        _linearOrder(orderArray);
        break;
    case OrderType::SKIPP:
        _skipOrder(orderArray);
        break;
    case OrderType::TWO_WAY:
        _twowayOrder(orderArray);
        break;
    case OrderType::MIXED:
        _mixedOrder(orderArray);
        break;
    case OrderType::GOLDEN:
        _linearOrder(orderArray);
        break;
    default:
        // Handle invalid order type
        break;
    }

    dVector angleArray = _fillAngles(orderArray);

    spiralArmAngles.resize(mNArms, dVector(mNEcho));

    for (int i = 0; i < mNArms; i++){
        for (int j = 0; j < mNEcho; j++){

            int index = i;

            if(mInCohMech)
                index = i + j * mNArms;

            spiralArmAngles[i][j] = angleArray[index];

            if (mDebug){
                printf("Arm: %d, Echo: %d, Angle: %.2f, index: %d\n", 
                i, j, 180/M_PI * angleArray[index], index);
            }
        }
    }
}

void SpiralArmOrder::composeOrder(iVector &spiralOrder)
{
    spiralOrder.resize(mTotalAngles);

    switch (mOrderType)
    {
    case OrderType::LINEAR:
        _linearOrder(spiralOrder);
        break;
    case OrderType::SKIPP:
        _skipOrder(spiralOrder);
        break;
    case OrderType::TWO_WAY:
        _twowayOrder(spiralOrder);
        break;
    case OrderType::MIXED:
        _mixedOrder(spiralOrder);
        break;
    case OrderType::GOLDEN:
        _linearOrder(spiralOrder);
        break;
    default:
        // Handle invalid order type
        break;
    }
}

}