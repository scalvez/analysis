#ifndef conf_sensitivity_h
#define conf_sensitivity_h

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

namespace conf_sens {

  const double year2sec = 3600 * 24 * 365.25;

  const double isotope_mass_number = 82;
  const double isotope_mass = 7; //kg
  const double exposure = 2.5; //years
  const double tracker_volume = 15.3; //m3

  const double T_2nu = 9.0e19; //years

  // source
  const double eff_0nu_2e = 2.634e-1; // MM
  // const double eff_0nu_2e = 0.161204; // rhc
  // const double eff_2nu_full_range_2e = 0.10104;
  // const double eff_2nu_2e = 2.1585e-1 / 25;
  // const double eff_2nu_2e = 8.967e-2; //full
  const double eff_2nu_2e = 8.967e-2 * 2.6435e-2; //full E>2 cut
  const double eff_tl208_2e = 1.058e-3;
  const double eff_bi214_2e = 1.642e-3;
  const double eff_radon_2e = 1.951e-4;

  // topo
  // const double eff_0nu_2e = 2.363e-1; // MM
  // // const double eff_0nu_2e = 0.161204; // rhc
  // // const double eff_2nu_full_range_2e = 0.10104;
  // const double eff_2nu_2e = 1.929e-1 / 25;
  // const double eff_tl208_2e = 8.67e-4;
  // const double eff_bi214_2e = 1.36e-3;
  // const double eff_radon_2e = 5.69e-5;

  // last old
  // // const double eff_0nu_2e = 0.25462; // MM
  // const double eff_0nu_2e = 0.161204; // rhc
  // // const double eff_2nu_full_range_2e = 0.10104;
  // const double eff_2nu_2e = 0.21712 / 25;
  // const double eff_tl208_2e = 0.00121521;
  // const double eff_bi214_2e = 0.00165497;
  // const double eff_radon_2e = 0.00024013;

  //old, 0nu and 2nu for LAL talk
  // const double eff_0nu_2e = 427876./(); //uBq/kg
  // const double eff_0nu_2e = 427876./(57*30000); //uBq/kg
  // const double eff_2nu_2e = 1031781./(124*40000)/25 ; //uBq/kg  // only correct above 2MeV
  // const double eff_tl208_2e = 10875./1e7; //uBq/kg
  // const double eff_bi214_2e = 29148./1e7; //uBq/kg

  // const double eff_radon_2e = 729./1e7; //150uBq/m3
  const double tl208_activity = 2e-6; //Bq/kg
  // const double tl208_activity = 23.6 /1e6; //Bq/kg
  const double bi214_activity = 10e-6; //Bq/kg
  // const double bi214_activity = 338. /1e6; //Bq/kg
  const double radon_activity = 150e-6; //Bq/m3

  const double k_sens = log(2.) * 6.022e23  * isotope_mass * 1000 * exposure / isotope_mass_number;

  // const double eff_0nu_bdt_window = 0.1366;
  // const double eff_0nu_roi_window = 0.1502;
  // const double eff_0nu_bdt_window = 0.25022*0.70972;
  // const double eff_0nu_roi_window = 0.25022*0.69416;

  const double N_source_2nu = eff_2nu_2e * k_sens / T_2nu;
  const double N_source_tl208 = tl208_activity * isotope_mass * exposure * year2sec * eff_tl208_2e;
  const double N_source_bi214 = bi214_activity * isotope_mass * exposure * year2sec * eff_bi214_2e;
  const double N_source_radon = radon_activity * tracker_volume * exposure * year2sec * eff_radon_2e;

}

#endif
