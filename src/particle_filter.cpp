/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>

#include "particle_filter.h"

//calculates a multivariate Gaussian
double mv_gauss(double x, double y, double mu_x, double mu_y, double sig_x, double sig_y) {
  return exp(-( (x-mu_x) * (x-mu_x) / (2.0*sig_x*sig_x) + (y-mu_y) * (y-mu_y) / (2.0*sig_y*sig_y))) / (2.0 * M_PI *sig_x*sig_y);
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 200;
  particles.resize(num_particles);
	weights.resize(num_particles);

	std::default_random_engine generator;

	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	for(int i=0; i < num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(generator);
		p.y = dist_y(generator);
		p.theta = dist_theta(generator);
    if(i == 0) {
      p.weight = 1.0;
      weights[i] = 1.0;
    } else {
		    p.weight = 0.0;
		    weights[i] = 0.0;
  }

		particles[i] = p;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	// 'for reference:
	// xf=x0+v/θ˙[sin(θ0+θ˙(dt))−sin(θ0)]
	// yf = y0 + v/θ˙[cos(θ0) - cos((θ0 + θ˙(dt))]
	// θf = θ0 + θ˙(dt)'
	std::default_random_engine generator;
	double x,y,theta;
		for(int i; i < particles.size(); i++) {

			//check and make sure yaw_rate is not zero and calc accordingly
			if (fabs(yaw_rate) == 0.0) {
        x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
        y =  particles[i].y + velocity*delta_t*sin(particles[i].theta);
        theta = particles[i].theta;
		} else {
        x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
        y = particles[i].y + -velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
        theta = yaw_rate * delta_t + particles[i].theta;
			}


			// add gaussian noise
			std::normal_distribution<double> dist_x(x, std_pos[0]);
			std::normal_distribution<double> dist_y(y, std_pos[1]);
			std::normal_distribution<double> dist_theta(theta, std_pos[2]);

			particles[i].x = dist_x(generator);
			particles[i].y = dist_y(generator);
			particles[i].theta = dist_theta(generator);
		}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	//std::cout << "Associate DATA!" << std::endl;
	for(int i = 0; i < observations.size(); i ++) {
		//set old_distance to keep track of smallest distance
		double old_distance = 10000000.0; //set with a big number at first to initialize
		for(int j = 0; j < predicted.size(); j++) {
			double calc_distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);//sqrt(pow((predicted[j].x - observations[i].x),2) + pow((predicted[j].y - observations[i].y) ,2) );
			if(calc_distance < old_distance) {
				old_distance = calc_distance;
				observations[i].id = predicted[j].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html


	//std::cout << "updating weights" << std::endl;

  //clear old weights
	weights.clear();

	for(int i = 0; i < particles.size(); i++) {
		//transform the coordinates to map coordinates
    std::vector<LandmarkObs> obs_map;
		for (int j=0; j < observations.size(); j++) {
      LandmarkObs observation;
			observation.x = 0.0; //make sure to initialize to zero
			observation.x += observations[j].x * cos(particles[i].theta);
			observation.x += -observations[j].y * sin(particles[i].theta);
			observation.x += particles[i].x;
			observation.y = 0.0; //init to zero
	    observation.y += observations[j].x * sin(particles[i].theta);
			observation.y += observations[j].y * cos(particles[i].theta);
			observation.y += particles[i].y;

			observation.id = -1; //needs an id yet???

			obs_map.push_back(observation);
	}

	//get predicted measurements
	std::vector<LandmarkObs> predicted;
	for(int j = 0; j < map_landmarks.landmark_list.size() ; j++) {
		double distance =  dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y);//sqrt(pow( (particles[i].x - map_landmarks.landmark_list[j].x_f) ,2) + pow( (particles[i].y - map_landmarks.landmark_list[j].y_f) ,2));
		if(distance <= sensor_range) {
				//std::cout << "getting predictedOb!" << map_landmarks.landmark_list[j].id_i << std::endl;
				LandmarkObs predictedOb;
				predictedOb.id = map_landmarks.landmark_list[j].id_i;
				predictedOb.x = map_landmarks.landmark_list[j].x_f;
				predictedOb.y = map_landmarks.landmark_list[j].y_f;

			 	predicted.push_back(predictedOb);
		}
	}

  //perform data association
	dataAssociation(predicted, obs_map);
	// std::cout << "Finished Association!" << std::endl;

	//update the weights
	double newWeight = 1.0; //initialize to 1
	double weight;

	for(int j = 0; j < predicted.size(); j++) {
		int closest_index = -1;
		double min_distance = 10000000.0; //set to large distance to  initialize

		for(int k = 0; k < obs_map.size(); k++) {

      if(predicted[j].id == obs_map[k].id){
				double distance = dist(predicted[j].x,predicted[j].y,obs_map[k].x,obs_map[k].y);//sqrt(pow((predicted[j].x - obs_map[k].x),2) + pow((predicted[j].y - obs_map[k].y) ,2) );

				if (distance < min_distance){
						closest_index = k;
						min_distance = distance;
				}
			}

		}
		if(closest_index != -1) {
			//std::cout << "closest_index: " << closest_index << std::endl;
			//calc multivariate gaussian terms
			// double firstTerm = ( 2.0 * M_PI * std_landmark[0]* std_landmark[1]);
			// std::cout << "First Term: " << firstTerm << std::endl;
			// double xTerm = pow( (predicted[j].x - obs_map[closest_index].x), 2) / ( 2*(pow(std_landmark[0], 2) ));
			// std::cout << "X Term: " << xTerm << std::endl;
			// double yTerm = (predicted[j].y - obs_map[closest_index].y)* / ( 2*(pow(std_landmark[1], 2) ));
			// std::cout << "Y Term: " << yTerm << std::endl;
			weight = mv_gauss(predicted[j].x, predicted[j].y, obs_map[closest_index].x, obs_map[closest_index].y, std_landmark[0], std_landmark[1]); //exp(-1.0*(xTerm + yTerm) ) / firstTerm;
			//std::cout << "weight i " << weight << std::endl;
			newWeight *= weight;
		}
	}

		// std::cout << "update with new weights" << std::endl;
		particles[i].weight = newWeight; //(1/(2*M_PI*std_landmark[0]*std_landmark[1]))*exp(-( (pow((predicted[i].x - obs_map[i].x),2)/(2*(pow(std_landmark[0],2)))) + (pow((predicted[i].y - obs_map[i].y),2)/(2*(pow(std_landmark[1],2)))) ) );
		weights.push_back(newWeight);
}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	//std::cout << "resample!" << std::endl;

  //From Udacity Forum suggestions on resampling with discrete_distribution fx:
	std::random_device rd;
  std::mt19937 generator_wts(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  std::vector<Particle> resampledParticles;

  resampledParticles.resize(num_particles);
  for (int i=0; i < num_particles; i++) {
      Particle p = particles[d(generator_wts)];
      resampledParticles.push_back(p);
  }

    particles = resampledParticles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
