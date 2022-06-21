#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#define PI 3.1415926536
#define n 1000					// Number of points.
#define epsilon 2.1/3.1					// The eccentricity of the orbit.
#define maxerror 0.0001				    // The maximum error allowed in the Newton's Method solution of Kepler's Equation
#define c 63239.7						// Speed of light in AU/yr

struct coordTuple {
	double x, y, z, t;
};

double eccentricAnomaly(double M);      // Function that given mean anomaly, approximates the Eccentric Anomaly using Newton's Method.
double thetaVals(double E);				// Function that given eccentric anomaly, finds the angular coordinate theta.
struct coordTuple polarToCartesian(double r, double theta, double time);			// Function that converts polar coordinates into cartesian coodinates.
struct coordTuple equatorialTransform(struct coordTuple coords);
double craftFreq(struct coordTuple velocity, struct coordTuple pulsar);


double main() {
	double g, a, angFreq, period, deltaT;
	g = 3.964 * pow(10.0, -14.0);       // Value of the product G*M in AU^3 / s^2 where M is the mass of the sun.
	a = 3.1; 							// Value of the semimajor axis in AU.
	struct coordTuple crab = { 0.7723, 0.5127, 0.3748, 187.502 };  // This gives the direction cosines of each pulsar, which are stored in x,y,z, while t contains the angular frequency(rad/s) of pulses.
	struct coordTuple pulsarTwo = { 0.6161, 0.4090, 0.6731, 2708.27 };
	struct coordTuple pulsarThree = { 0.0848, -0.9049, -0.4290, 4027.68 };

	angFreq = sqrt(g / pow(a, 3));
	period = (2 * PI) / angFreq;
	deltaT = period / ((double)n);

	double angles[n];					// Initializing an array containing values of the angular coordinate theta.
	struct coordTuple xytCoords[n];		// Initiialzing an array to contain coordinates (x,y,t) of points on the orbit.
	double radialDist[n];				// Initializing an array containing the values of the radial distance corresponding to time values of the same index in time[]
	double psi[n];						// Initializing an array containing values of the eccentric anomaly corresponding to time values with the same index in time[]
	double time[n];						// Initializing an array containing n equally spaced time values from t = 0 to t= period.
	struct coordTuple realCoords[n];
	time[0] = 0;
	int i;
	for (i = 1; i < n; i++) {
		time[i] = deltaT + time[i - 1];
	}
	int k;
	for (k = 0; k < n; k++) {
		psi[k] = eccentricAnomaly(angFreq * time[k]);
		angles[k] = thetaVals(psi[k]);
		radialDist[k] = a * (1 - epsilon * cos(psi[k]));
		xytCoords[k] = polarToCartesian(radialDist[k], angles[k], time[k]);
		realCoords[k] = equatorialTransform(xytCoords[k]);

	}
	double xDerivative[n];				// Initializing an array containing the derivatives at each point of x w.r.t time.
	double yDerivative[n];				// Initializing an array containing the derivatives at each point of y w.r.t time.
	double zDerivative[n];
	struct coordTuple velocities[n];

	xDerivative[0] = (realCoords[1].x - realCoords[0].x) / (deltaT * 3.17098 * pow(10, -8));				// The derivatives at the left and right endpoints are approximated using a
	xDerivative[n - 1] = (realCoords[n - 1].x - realCoords[n - 2].x) / (deltaT * 3.17098 * pow(10, -8));	// forward finite difference and a backward finite difference, respectively with the remaining points
	yDerivative[0] = (realCoords[1].y - realCoords[0].y) / (deltaT * 3.17098 * pow(10, -8));				// having their derivatives approximated by a central finite difference.
	yDerivative[n - 1] = (realCoords[n - 1].y - realCoords[n - 2].y) / (deltaT * 3.17098 * pow(10, -8));
	zDerivative[0] = (realCoords[1].z - realCoords[0].z) / (deltaT * 3.17098 * pow(10, -8));
	zDerivative[n - 1] = (realCoords[n - 1].z - realCoords[n - 2].z) / (deltaT * 3.17098 * pow(10, -8));

	int m;
	for (m = 1; m < n - 1; m++) {
		xDerivative[m] = (realCoords[m + 1].x - realCoords[m - 1].x) / (2 * (deltaT * 3.17098 * pow(10, -8)));
		yDerivative[m] = (realCoords[m + 1].y - realCoords[m - 1].y) / (2 * (deltaT * 3.17098 * pow(10, -8)));
		zDerivative[m] = (realCoords[m + 1].z - realCoords[m - 1].z) / (2 * (deltaT * 3.17098 * pow(10, -8)));
	}

	for (m = 0; m < n; m++) {
		velocities[m].x = xDerivative[m];
		velocities[m].y = yDerivative[m];
		velocities[m].z = zDerivative[m];
		velocities[m].t = 1;
	}


	int j;
	printf("c = [");
	for (j = 0; j < n; j++) {
		printf("%lf,", realCoords[j].z);
	}
	printf("];");


	/*	printf("Coordinates (x,y,z) of points in the orbit in AU with Components of velocity (v_x,v_y,v_z) in AU/yr, time in years, and observed frequencies in rad/s \n");
		int j;
		for (j = 0; j < n; j++) {
			printf("Position: (%lf, %lf, %lf) Velocity: (%lf, %lf, %lf) Time = %lf Relative Frequency: %lf \n", realCoords[j].x, realCoords[j].y, realCoords[j].z, velocities[j].x, velocities[j].y, velocities[j].z, realCoords[j].t * 3.17098 * pow(10, -8), craftFreq(velocities[j], crab));
		}
		printf("\n");
	*/

	system("pause");
	return 0;
}

double eccentricAnomaly(double M) {
	double E = M;		// initial guess
	double fE = E - epsilon * sin(E) - M;
	while (fabs(fE) > maxerror) {
		E = E - (E - epsilon * sin(E) - M) / (1 - epsilon * cos(E));
		fE = E - epsilon * sin(E) - M;
	}
	return E;
}

double thetaVals(double E) {
	double theta;
	if (E > PI) {
		theta = acos(-(cos(E) - epsilon) / (1 - epsilon * cos(E))) + PI;
	}
	else {
		theta = acos((cos(E) - epsilon) / (1 - epsilon * cos(E)));
	}
	return theta;
}

struct coordTuple polarToCartesian(double r, double theta, double time) {
	double x, y, z;
	x = r * cos(theta);
	y = r * sin(theta);
	z = 0;
	struct coordTuple coordinates = { x, y, z, time };
	return coordinates;

}

struct coordTuple equatorialTransform(struct coordTuple coords) {
	struct coordTuple newCoords;
	double alpha = 1.2714;
	double obliq = 23.45*PI/180;

	newCoords.t = coords.t;
	newCoords.x = coords.x * cos(alpha) - coords.y * sin(alpha);
	newCoords.y = coords.x * sin(alpha) * cos(obliq) + coords.y * cos(alpha) * cos(obliq) - coords.z * sin(obliq);
	newCoords.z = coords.x * sin(alpha) * sin(obliq) + coords.y * cos(alpha) * sin(obliq) + coords.z * cos(obliq);

	return newCoords;
}

double craftFreq(struct coordTuple velocity, struct coordTuple pulsar) {
	double relFreq = pulsar.x * velocity.x / c + pulsar.y * velocity.y / c + pulsar.z * velocity.z / c;
	return relFreq;
}