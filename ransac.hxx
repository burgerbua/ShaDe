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
    template <typename T, int N>
    class Model {
    protected:
        std::array<T, N> min_model_params;
        
        virtual double compute_distance_measure(const T& param) const = 0;
        
    public:
        typedef T param_type;
        const static int num_params = N;
        typedef std::shared_ptr<Model> ptr;

        std::pair<double, std::vector<T>> evaluate(typename std::vector<T>::iterator itb, typename std::vector<T> ite, const double threshold) {
            std::vector<T> inliers;
            std::copy_if(itb, ite, std::back_inserter(inliers),
                         [&](const T& param){return this->compute_distance_measure(param)<threshold;});
            
            const double inlier_fraction = static_cast<double>(inliers.size()) / static_cast<double>(std::distance(itb, ite));
            return std::make_pair(inlier_fraction, inliers);
        }
    };

    template <class M>
    class RANSAC {
        typename M::ptr best_model; // pointer to the best model, valid only after estimate() is called
        std::vector<typename M::param_type> best_inliers;
        double best_model_score; // the score of the best model
        
        const double threshold; // the threshold for computing model consensus
        const int max_iterations; // number of iterations before termination
        
        void reset() {
            best_model = nullptr;
            best_inliers.clear();
            best_model_score = 0.0;
        }
        
    public:
        RANSAC(const double _threshold, const int _max_iterations = 1000)
        : threshold(_threshold), max_iterations(_max_iterations) {
            reset();
        }
        
        typename M::ptr get_best_model() const {
            return best_model;
        }
        
        const std::vector<typename M::param_type>& get_best_inliers() const {
            return best_inliers;
        }
        
        const double get_best_model_score() const {
            return best_model_score;
        }
        
        bool estimate(const std::vector<typename M::param_type>& params) {
            
            reset();

            if(params.size() <= M::num_params)
            {
                std::cout << "[ WARN ]: RANSAC - Number of data points is too less. Not doing anything." << std::endl;
                return false;
            }
            
            // for creating random models
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<size_t> dis(0, params.size()-1);

            // find best of the random models
            best_model_score = 0.0;
            for (int i=0; i<max_iterations; ++i) {
                std::vector<typename M::param_type> random_params;
                for (size_t i=0; i<M::num_params; ++i)
                    random_params.emplace_back(params.at(dis(gen)));
                typename M::ptr random_model = std::make_shared<M>(random_params);
                
                std::pair<double, std::vector<typename M::param_type>> eval_pair = random_model->evaluate(std::begin(params), std::end(params), threshold);
                
                if (eval_pair.first>best_model_score) {
                    best_model_score = eval_pair.first;
                    best_inliers = std::move(eval_pair.second);
                    best_model = random_model;
                }
            }
            
            return true;
        }
    };
} // namespace shade

#endif /* ransac_h */
