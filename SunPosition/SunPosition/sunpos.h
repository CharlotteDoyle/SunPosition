// sunpos.h :
// 
// this header file contains the declaration of a class
// which includes all the input and output data,
// and the function that performs the calculation
//
// to calculate the sun position, follow these steps:
//
// 1.	include this file
// 2.	declare a variable of type SunCoord
// 3.	initialize the variavle giving the 9 input quantities
//		required. This can be done in the declaration,
//		listing the quantities between colons, or
//		calling the function SetCoord(). In both
//		cases, only the first 4 quantities are necessary;
//		the others are set to standard values
//		(Pressure = 1 atm, T = 20°„C£¨
//		0 all the other quantities) if omitted.
//		Longitude and Latitude must be given in radians,
//		Pressure in atm,
//		Temperature in °„C,
//		Day, Month and Year are integer,
//		UT in decimal hour from 0 to 24
//		(i.e. 3h30m p.m. becomes 15.5)
// 4.	callthe Calculate() member of the variable SunCoord
//
// example:
// the following program lines define a variable X
// of type SunCoord, initialize X with the nine 
// quantities required (UT, Day, etc. must have been
// previously defined), perform the calculation, and output 
// to the screen the Zenith and the Azimuth;
// then X is changed via X.SetCoord() and the calculation
// is repeated
// 
//		{ ...
//		
//			SunCoord X = { UT, Day, Month, Year, Delta_t, 
//				ObserverLatitude, ObserverLongitude,
//				Pressure, Temperature };
//			
//			X.Calculate();
//			cout << "Zenith = " << X.Zenith << "\tAzimuth = "
//				<< X.Azimuth << endl;
//			...
//			X.SetCoord(UT2, Day2, Month2, Year2, Delta_t2);
//			X.Calculate();
//			cout << "Zenith2 = " << x.Zenith
//				<< "\tAzimuth2 = " << X.Azimuth << endl;
//
//		... }
//
// Warning: in order to improve accessibility and efficiency,
// there is no access control in the class. The user
// is free to directly access and modify all the data,
// and there is no any control of consistency. Some caution
// in the use of the class is advisable.

#ifndef SUNPOS_H
#define SUNPOS_H

#include <cmath>

#ifndef PI
#define PI 3.14159265358979
#endif

// declaration of the class SunCoord

class SunCoord
{
public: // no access control

// input data
	double UT;
	int Day;
	int Month;
	int Year;
	double Delta_t;
	double ObserverLatitude;
	double ObserverLongitude;
	double Pressure;
	double Temperature;

// output data
	// global:
	double JD;
	double JDE;
	double HeliocLongitude;
	double GeocSolarLongitude;
	double RightAscension;
	double Declination;

	// local:

	double HourAngle;
	double TopocRightAscension;
	double TopocDeclination;
	double TopocHourAngle;
	double Elevation_no_refrac;
	double RefractionCorrection;
	double Zenith;
	double Azimuth;

// functions
	// SetCoord: inline function which allows to give
	// the input data quickly;

	void SetCoord(double dUT, int dDay, int dMonth,
		int dYear, double dDelta_t = 0,
		double dObserverLatitude = 0,
		double dObserverLongitude = 0,
		double dPressure = 1, double dTemperature = 20)
	{
		UT = dUT;
		Day = dDay;
		Month = dMonth;
		Year = dYear;
		Delta_t = dDelta_t;
		ObserverLatitude = dObserverLatitude;
		ObserverLongitude = dObserverLongitude;
		Pressure = dPressure;
		Temperature = dTemperature;
	}
	
	// Calculate: the algorithm
	void Calculate();
};

// Calculate() definition: the algorithm

inline void SunCoord::Calculate()
{
	// Calculation of JD and JDE:

	double dYear, dMonth;
	if (Month <= 2)
	{
		dYear = static_cast<double>(Year) - 1.0;
		dMonth = static_cast<double>(Month) + 12.0;
	}
	else
	{
		dYear = static_cast<double>(Year);
		dMonth = static_cast<double>(Month);
	}

	auto JD_t = static_cast<double>(trunc(365.25 * (dYear - 2000)))
		+ static_cast<double>(trunc(30.6001 * (dMonth + 1)))
		+ static_cast<double>(Day) + UT / 24.0 - 1158.5;

	auto t = JD_t + Delta_t / 86400;

	// standard JD and JDE (useless for the
	// computation, they are computed
	// for completeness)

	JDE = t + 2452640;
	JD = JD_t + 2452640;

	// HELIOCENTRIC LONGITUDE

	// linear increase + annual harmonic

	auto ang = 1.72019e-2 * t - .0563;
	HeliocLongitude = 1.740940 + 1.7202768683e-2 * t +
		3.34118e-2 * sin(ang) + 3.488e-4 * sin(2 * ang);

	// Moon perturbation

	HeliocLongitude += 3.13e-5 * sin(2.127730e-1 * t - .585);

	// Harmonic correction

	HeliocLongitude += 1.26e-5 * sin(4.243e-3 * t + 1.46) +
		2.35e-5 * sin(1.0727e-2 * t + .72) +
		2.76e-5 * sin(1.5799e-2 * t + 2.35) +
		2.75e-5 * sin(2.1551e-2 * t - 1.98) +
		1.26e-5 * sin(3.1490e-2 * t - .80);

	// Polynomial correction

	auto t2 = t / 1000;
	HeliocLongitude += (((-2.30796e-07 * t2 + 3.7976e-06) * t2
		- 2.0458e-5) * t2 + 3.976e-05) * t2 * t2;

	// to obtain the Heliocentric Longitude in the range
	// [0, 2 * PI], uncomment the following line:

	//HeliocLongitude = fmod(HeliocLongitude, 2 * PI);

	// END HELIOCENTRIC LONGITUDE CALCULATION

	// Correction to longitude due to nutation

	auto delta_psi = 8.33e-5 * sin(9.52e-4 * t - 1.173);

	// Earth axis inclination

	auto epsilon = -6.21e-9 * t + 0.409086 +
		4.46e-5 * sin(9.252e-4 * t + .397);

	// Geocentric global solar coordinates:

	GeocSolarLongitude = HeliocLongitude + PI + delta_psi - 9.932e-5;

	auto s_lambda = sin(GeocSolarLongitude);

	RightAscension = atan2(s_lambda* cos(epsilon),
		cos(GeocSolarLongitude));

	Declination = asin(sin(epsilon) * s_lambda);

	// local hour angle

	HourAngle = 6.30038809903 * JD_t + 4.8824623 +
		delta_psi * .9174 + ObserverLongitude -
		RightAscension;

	// to obtain the local hour angle in the range
	// [0, 2 * PI], uncomment the following line:

	//HourAngle = fmod(HourAngle, 2 * PI);

	auto c_lat = cos(ObserverLatitude);
	auto s_lat = sin(ObserverLatitude);
	auto c_H = cos(HourAngle);
	auto s_H = sin(HourAngle);

	// Parallax correction to Right Ascensition:

	auto d_alpha = -4.26e-5 * c_lat * s_H;
	TopocRightAscension = RightAscension + d_alpha;
	TopocHourAngle = HourAngle - d_alpha;

	// Parallax correction to Declination:

	TopocDeclination = Declination - 4.26e-5 * (s_lat -
		Declination * c_lat);
	auto s_delta_corr = sin(TopocDeclination);
	auto c_delta_corr = cos(TopocDeclination);
	auto c_H_corr = c_H + d_alpha * s_H;
	auto s_H_corr = s_H - d_alpha * c_H;

	// Solar elevation angle, without
	// refraction correction:
	Elevation_no_refrac =
		asin(s_lat * s_delta_corr +
			c_lat * c_delta_corr * c_H_corr);

	// Refraction correction: it is calculated only
	// if Elevation_no_refrac > elev_min

	const auto elev_min = -.01;

	if (Elevation_no_refrac > elev_min)
		RefractionCorrection =
		.084217 * Pressure / (273 + Temperature) /
		tan(Elevation_no_refrac + .0031376 /
			(Elevation_no_refrac + .089186));
	else
		RefractionCorrection = 0;

	// Local coordinates of the sun:

	Zenith = PI / 2 - Elevation_no_refrac -
		RefractionCorrection;
	Azimuth = atan2(s_H_corr, c_H_corr * s_lat -
		s_delta_corr / c_delta_corr * c_lat);
}

#endif
