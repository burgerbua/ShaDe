//
//  ransac.hxx
//  ShaDe
//
//  Created by Matthias Messner on 4/21/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#ifndef ransac_h
#define ransac_h

#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <memory>
#include <algorithm>
#include <vector>
//#include <omp.h>

namespace shade
{
	// abstract model type for generic RANSAC model fitting
	template <typename T, int N> // minimum number of parameters required to define this model
	class Model
	{
	protected:
		std::array<T, N> min_model_params;

		virtual double compute_distance_measure(const T& param) const = 0;

	public:
		using param_type = T;
		const static int num_params = N;
		using ptr = std::shared_ptr<Model>;

		virtual std::pair<double, std::vector<T>> evaluate(
			const std::vector<T>& params,
			double threshold)
		{
			std::vector<T> inliers;
			for (auto& param : params) 
			{
				if (compute_distance_measure(param) < threshold)
				{
					inliers.push_back(param);
				}
			}
			const double inlier_fraction = static_cast<double>(inliers.size()) / static_cast<double>(params.size()); // This is the inlier fraction
			return std::make_pair(inlier_fraction, inliers);
		};
	};

	// M - AbstractModel
	template <class M>
	class RANSAC
	{
		std::vector<typename M::ptr> sampled_models; // vector of all sampled models
		typename M::ptr best_model; // pointer to the best model, valid only after estimate() is called
		std::vector<typename M::param_type> best_inliers;

		int max_iterations; // number of iterations before termination
		double threshold; // the threshold for computing model consensus
		double best_model_score; // the score of the best model
		int best_model_idx;

		std::mt19937 rand_engine; // mersenne twister high quality RNG that support *OpenMP* multi-threading

	public:
		RANSAC() : rand_engine(std::random_device())
		{
			reset();
		}

		virtual ~RANSAC() {}

		void reset()
		{
			// clear sampled models, etc. and prepare for next call.
			// reset RANSAC estimator state
			sampled_models.clear();

			best_model_idx = -1;
			best_model_score = 0.0;
		}

		void initialize(
			double _threshold, 
			int _max_iterations = 1000)
		{
			threshold = _threshold;
			max_iterations = _max_iterations;
		}

		typename M::ptr get_best_model() 
		{
			return best_model;
		}

		const std::vector<typename M::param_type>& get_best_inliers()
		{
			return best_inliers;
		}

		bool estimate(const std::vector<typename M::param_type>& data)
		{
			if (data.size() <= M::num_params)
			{
				std::cout << "[ WARN ]: RANSAC - Number of data points is too less. Not doing anything." << std::endl;
				return false;
			}

			std::vector<double> inlier_fraction_accum(max_iterations);
			std::vector<std::vector<typename M::ptr>> inliers_accum(max_iterations);
			sampled_models.resize(max_iterations);

			for (int i = 0; i < max_iterations; ++i)
			{
				// select M::num_params random samples
				std::vector<typename M::param_type> random_samples(M::num_params);
				std::vector<typename M::param_type> remainder_samples = data; // without the chosen random samples

				std::shuffle(remainder_samples.begin(), remainder_samples.end(), rand_engine); // to avoid picking the same element more than once
				std::copy(remainder_samples.begin(), remainder_samples.begin() + M::num_params, random_samples.begin());
				remainder_samples.erase(remainder_samples.begin(), remainder_samples.begin() + M::num_params);

				//typename M::ptr random_model = std::make_shared<M>(random_samples);
				typename M::ptr random_model = std::make_shared<M>(M::make_model(random_samples));

				// check if the sampled model is the best so far
				std::pair<double, std::vector<typename M::param_type>> eval_pair = random_model->evaluate(remainder_samples, threshold);
				inlier_fraction_accum[i] = eval_pair.first;
				inliers_accum[i] = eval_pair.second;

				// push back into history, could be removed later
				sampled_models[i] = random_model;
			}

			for (int i = 0; i < max_iterations; ++i)
			{
				if (inlier_fraction_accum[i] > best_model_score)
				{
					best_model_score = inlier_fraction_accum[i];
					best_model_idx = sampled_models.size() - 1;
					best_model = sampled_models[i];
					best_inliers = inliers_accum[i];
				}
			}

			reset();

			return true;
		}
	};
} // namespace shade

#endif /* ransac_h */
