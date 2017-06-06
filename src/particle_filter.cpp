#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 10;
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];
    std::default_random_engine generator;
    std::normal_distribution<double> dist_x(x, std_x);
    std::normal_distribution<double> dist_y(y, std_y);
    std::normal_distribution<double> dist_theta(theta, std_theta);
    for (int i=0; i<num_particles; i++)
    {
        Particle p;
        p.id = i;
        p.x = dist_x(generator);
        p.y = dist_y(generator);;
        p.theta = dist_theta(generator);
        p.weight = 1.0;
        particles.push_back(p);
    }
    is_initialized = false;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];
    std::default_random_engine generator;
    for (int i=0; i<num_particles; i++)
    {
        Particle p = particles[i];
        double dtheta = yaw_rate * delta_t;
        std::normal_distribution<double> dist_theta(dtheta, std_theta);
        double theta_f = p.theta + dtheta + dist_theta(generator);
        p.theta = theta_f;

        double dx = velocity * (sin(theta_f) - sin(p.theta)) / yaw_rate;
        std::normal_distribution<double> dist_x(dx, std_x);
        p.x += dx + dist_x(generator);

        double dy = velocity * (cos(p.theta) - cos(theta_f)) / yaw_rate;
        std::normal_distribution<double> dist_y(dy, std_y);
        p.y += dy + dist_y(generator);
    }



}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void transform_coord(Particle p, std::vector<LandmarkObs> &obs)
{
    double x = p.x;
    double y = p.y;
    for (int i=0; i<obs.size(); i++)
    {
        obs[i].x = p.x + obs[i].x * cos(p.theta) - obs[i].y * sin(p.theta);
        obs[i].y = p.y + obs[i].x * sin(p.theta) + obs[i].y * cos(p.theta);
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) 
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    std::cout << "number of observations = " << observations.size() << std::endl;
    std::vector<LandmarkObs>::iterator it; 
    for (int i=0; i<num_particles; i++)
    {
        Particle p = particles[i];
        //std::cout << "particle number = " << i << "  nobs = " << observations.size() << std::endl;
        //for(it=observations.begin() ; it < observations.end(); it++ )
            //std::cout << (*it).x << "\t" << (*it).y << std::endl;

        transform_coord(p, observations);
        //std::cout << "after transform:\n";
        //for(it=observations.begin() ; it < observations.end(); it++ )
            //std::cout << (*it).x << "\t" << (*it).y << std::endl;


    }

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
