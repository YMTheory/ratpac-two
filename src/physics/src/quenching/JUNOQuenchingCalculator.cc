#include <RAT/JUNOQuenchingCalculator.hh>

JUNOQuenchingCalculator::JUNOQuenchingCalculator()
 {}


double JUNOQuenchingCalculator::QuenchedEnergyDeposit(const double TotalEnergyDeposit, const double delta,  const double kB, const double kC) {
    double QuenchedEnergyDeposit = TotalEnergyDeposit / (1 + kB * delta + kC * delta * delta);
    return QuenchedEnergyDeposit;

}
