#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <assert.h>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    //mytest();
    num_particles = 100;
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
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];
    double newx, newy, theta_f;
    std::default_random_engine generator;
    for (int i=0; i<num_particles; i++)
    {
        Particle p = particles[i];
        if (fabs(yaw_rate) < 1.0e-5)
        {
            newx = p.x + velocity * cos(p.theta) * delta_t;
            newy = p.y + velocity * sin(p.theta) * delta_t;
            theta_f = p.theta;
        }
        else
        {
            double ratio = velocity / yaw_rate;

            theta_f = p.theta + yaw_rate * delta_t;
            newx = p.x + (sin(theta_f) - sin(p.theta)) * ratio;
            newy = p.y + (cos(p.theta) - cos(theta_f)) * ratio;
        }

        std::normal_distribution<double> dist_x(newx, std_x);
        std::normal_distribution<double> dist_y(newy, std_y);
        std::normal_distribution<double> dist_theta(theta_f, std_theta);

        particles[i].x = dist_x(generator);
        particles[i].y = dist_y(generator);
        particles[i].theta = dist_theta(generator);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) 
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

std::vector<LandmarkObs> transform_coord_list(Particle p, std::vector<LandmarkObs> obs)
{
    for (int i=0; i<obs.size(); i++)
    {
        double x, y;
        x = p.x + obs[i].x * cos(p.theta) - obs[i].y * sin(p.theta);
        y = p.y + obs[i].x * sin(p.theta) + obs[i].y * cos(p.theta);
        obs[i].x = x;
        obs[i].y = y;
    }
    return obs;
}

void transform_coord(Particle p, LandmarkObs &obs)
{
    double x, y;
    x = p.x + obs.x * cos(p.theta) - obs.y * sin(p.theta);
    y = p.y + obs.x * sin(p.theta) + obs.y * cos(p.theta);
    obs.x = x;
    obs.y = y;
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
    
    std::vector<LandmarkObs>::iterator it; 
    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];

    for (int i=0; i<num_particles; i++)
    {
        Particle p = particles[i];

        std::vector<int> landmark_id_list;
        std::vector<double> landmark_x_list;
        std::vector<double> landmark_y_list;

        std::vector<LandmarkObs> obs = transform_coord_list(p, observations);
        for (int j=0; j<obs.size(); j++)
        {
            double min_dist = 1.0e5;
            int closest_landmark;
            double x1 = obs[j].x;
            double y1 = obs[j].y;
            for (int l=0; l<map_landmarks.landmark_list.size(); l++)
            {
                //double x1 = p.x;
                //double y1 = p.y;
                double x2 = map_landmarks.landmark_list[l].x_f;
                double y2 = map_landmarks.landmark_list[l].y_f;
                double d = dist(x1, y1, x2, y2);
                if ( (d < min_dist) & (d < sensor_range) )
                {
                    min_dist = d;
                    closest_landmark = l;
                }
            }
            landmark_id_list.push_back(map_landmarks.landmark_list[closest_landmark].id_i);
            landmark_x_list.push_back(map_landmarks.landmark_list[closest_landmark].x_f);
            landmark_y_list.push_back(map_landmarks.landmark_list[closest_landmark].y_f);
        }
        p = SetAssociations(p, landmark_id_list, landmark_x_list, landmark_y_list);
        assert (observations.size() == landmark_id_list.size());
        assert (observations.size() == p.associations.size());

        double sum = 0;
        for (int j=0; j<obs.size(); j++)
        {
            LandmarkObs iobs = obs[j];
            sum += pow(p.sense_x[j] - iobs.x, 2) / (2 * sigma_x * sigma_x);
            sum += pow(p.sense_y[j] - iobs.y, 2) / (2 * sigma_y * sigma_y);
        }
        p.weight = exp(-sum) / pow(2 * M_PI * sigma_x * sigma_y, obs.size());
        particles[i] = p;
    }

}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    std::vector<double> w_dist;
    std::vector<int> id_dist;

    for (int i=0; i<num_particles; i++)
        w_dist.push_back(particles[i].weight);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> distribution(w_dist.begin(), w_dist.end());
    std::map<int, int> m;
    for(int n=0; n<num_particles; ++n)
        ++m[distribution(gen)];

    std::vector<Particle> refresh_p;
    for (int i=0; i<num_particles; i++)
    {
        int c = m[i];
        while (c > 0)
        {
            refresh_p.push_back(particles[i]);
            c--;
        }
    }
    assert (refresh_p.size() == particles.size());
    particles = refresh_p;

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

void ParticleFilter::mytest()
{
    Particle p;
    p.id = 1;
    p.x = 4.0;
    p.y = 5.0;
    p.theta = -M_PI/2;
    p.weight = 0;

    std::vector<LandmarkObs> obs;
    LandmarkObs o;
    o.id = 1;
    o.x = 2;
    o.y = 2;
    obs.push_back(o);
    o.id = 1;
    o.x = 3;
    o.y = -2;
    obs.push_back(o);
    o.id = 1;
    o.x = 0;
    o.y = -4;
    obs.push_back(o);

    std::vector<LandmarkObs> lms;
    LandmarkObs l;
    l.id = 1;
    l.x = 5;
    l.y = 3;
    lms.push_back(l);
    l.id = 2;
    l.x = 2;
    l.y = 1;
    lms.push_back(l);
    l.id = 3;
    l.x = 6;
    l.y = 1;
    lms.push_back(l);
    l.id = 4;
    l.x = 7;
    l.y = 4;
    lms.push_back(l);
    l.id = 5;
    l.x = 4;
    l.y = 7;
    lms.push_back(l);
    double sensor_range = 50;
    double sigma_x = 0.3;
    double sigma_y = 0.3;


    obs = transform_coord_list(p, obs);
    for (int i=0; i<obs.size(); i++)
    {
        std::cout << i << ": " << obs[i].x << ", " << obs[i].y << std::endl;
    }



        std::vector<int> landmark_id_list;
        std::vector<double> landmark_x_list;
        std::vector<double> landmark_y_list;

        for (int j=0; j<obs.size(); j++)
        {
            double min_dist = 1.0e5;
            int closest_landmark;
            double x1 = obs[j].x;
            double y1 = obs[j].y;
            for (int l=0; l<lms.size(); l++)
            {
                //double x1 = p.x;
                //double y1 = p.y;
                double x2 = lms[l].x;
                double y2 = lms[l].y;
                double d = dist(x1, y1, x2, y2);
                if ( (d < min_dist) & (d < sensor_range) )
                {
                    min_dist = d;
                    closest_landmark = l;
                }
            }
            landmark_id_list.push_back(lms[closest_landmark].id);
            landmark_x_list.push_back(lms[closest_landmark].x);
            landmark_y_list.push_back(lms[closest_landmark].y);
        }
        p = SetAssociations(p, landmark_id_list, landmark_x_list, landmark_y_list);
        assert (obs.size() == landmark_id_list.size());
        assert (obs.size() == p.associations.size());
        std::cout << "assosiations: \n";

        double sum = 0;
        for (int j=0; j<obs.size(); j++)
        {
            LandmarkObs iobs = obs[j];
            std::cout << j << ": " << p.sense_x[j] << ", " << p.sense_y[j] << std::endl;
            sum += pow(p.sense_x[j] - iobs.x, 2) / (2 * sigma_x * sigma_x);
            sum += pow(p.sense_y[j] - iobs.y, 2) / (2 * sigma_y * sigma_y);
        }
        p.weight = exp(-sum) / pow(2 * M_PI * sigma_x * sigma_y, obs.size());
        std::cout << "weight = " << p.weight << std::endl;
    assert (1==2);
}

