#ifndef STATS_TOOLS_H
#define STATS_TOOLS_H 1

#include "TMath.h"

// Continuous FeldmanCousin functions computed by M. Bongrand
  // <bongrand@lal.in2p3.fr>
  double get_number_of_excluded_events(const double number_of_events_)
  {
    double number_of_excluded_events = 0.0;
    if (number_of_events_ < 29.0)
      {
        double x = number_of_events_;
        number_of_excluded_events =
          2.5617 + 0.747661 * x - 0.0666176 * std::pow(x,2)
          + 0.00432457 * std::pow(x,3) - 0.000139343 * std::pow(x,4)
          + 1.71509e-06 * std::pow(x,5);
      }
    else
      {
        number_of_excluded_events = 1.64 * std::sqrt(number_of_events_);
      }
    return number_of_excluded_events;
  }

#endif
