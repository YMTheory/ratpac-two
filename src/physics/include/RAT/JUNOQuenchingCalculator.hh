/////////////
// Birks' Law used in JUNO offline
// Added by Miao Yu
/////////////


#ifndef __JUNOQuenchingCalculator__
#define __JUNOQuenchingCalculator__


class JUNOQuenchingCalculator {

    public:
        JUNOQuenchingCalculator();

        double QuenchedEnergyDeposit(const double TotalEnergyDeposit, const double delta, const double kB, const double kC);



};


#endif
