/*
 * Cluster.hpp
 *
 *  Created on: Dec 30, 2011
 *      Author: yasir
 *      Given a set of loop closuers, clusters them
 */

#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_

#include <iostream>
#include <vector>
#include <cstdlib>

struct cluster
{
	int startLow, startHigh;
	int endLow, endHigh;
	int size;

	cluster(): startLow(-1), startHigh(-1), endLow(-1), endHigh(-1), size(0) {}
	cluster(int start, int end) : startLow(start), startHigh(start), endLow(end), endHigh(end), size(1){}
};

class Clusterizer
{
	std::vector<cluster> _clustersFound;
public:
	// Assumes that in the vector the values are passed as (start_1,end_1), (start_2,end_2), ...
	void clusterize( const std::vector<int> loops , const int threshold, std::vector<int>& membership, int& clusterCount)
	{
		if(loops.size() < 2)
		{
			std::cerr<<"clusterize(): "<<__LINE__<<" no loops to make clusters"<<std::endl;
			clusterCount = 0;
			membership = std::vector<int>();
			return;
		}
		_clustersFound.clear();
		membership.clear();

		for(size_t i=0; i<loops.size();i+=2)
		{
			int start 	= loops[i];
			int end 	= loops[i+1];

			if(_clustersFound.empty())
			{
				cluster s(start,end);
				_clustersFound.push_back(s);
				membership.push_back(_clustersFound.size()-1);

			}
			else
			{
				cluster& currentCluster = _clustersFound.back();
				if(
						(
							std::abs(start-currentCluster.startHigh)<threshold or
							std::abs(start-currentCluster.startLow)<threshold
						)
						and
						(
							std::abs(end-currentCluster.endHigh)<threshold or
							std::abs(end-currentCluster.endLow)<threshold
						)
				)
				{
					currentCluster.size++;
					membership.push_back(_clustersFound.size()-1);
					if(start<currentCluster.startLow)	currentCluster.startLow = start;
					if(start>currentCluster.startHigh)	currentCluster.startHigh = start;

					if(end<currentCluster.endLow) currentCluster.endLow = end;
					if(end>currentCluster.endHigh) currentCluster.endHigh = end;
				}
				else
				{
					cluster s(start,end);
					_clustersFound.push_back(s);
					membership.push_back(_clustersFound.size()-1);
				}

			}
		}
		if(0)
		{
			std::cout<<" \% Clusters formed "<<_clustersFound.size()<<std::endl;
			std::cout<<"limits = [ "<<std::endl;
			for(size_t i=0 ; i< _clustersFound.size() ; i++)
			{
				//std::cout<<_clustersFound[i].size<<" :"<<std::endl;
				std::cout<<" "<<_clustersFound[i].startLow<<" "<<_clustersFound[i].startHigh<<" ";
				std::cout<<" "<<_clustersFound[i].endLow<<" "<<_clustersFound[i].endHigh<<std::endl;;

			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;


			std::cout<<"membership =[ ";
			for(size_t i=0; i<membership.size();i++)
			{
				std::cout<<membership[i]<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;
		}
		clusterCount = _clustersFound.size();
	}

};

#endif /* CLUSTER_HPP_ */
