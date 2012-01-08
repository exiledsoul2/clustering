/*
 * OptCluster.hpp
 *
 *  Created on: Jan 1, 2012
 *      Author: yasir
 */

#ifndef OPTCLUSTER_HPP_
#define OPTCLUSTER_HPP_

#include <g2o/core/graph_optimizer_sparse.h>
#include <g2o/core/block_solver.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/types/slam2d/vertex_se2.h>

#include <fstream>

#include <boost/math/distributions/chi_squared.hpp>

#include <Cluster.hpp>

using namespace g2o;



class GraphStorage
{
	private:

	struct _edge{
		int from;
		int to;
		g2o::SE2 transform;
		Eigen::Matrix3d information;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	std::vector<g2o::SE2, Eigen::aligned_allocator<g2o::SE2> > _vertices;
	std::vector< _edge , Eigen::aligned_allocator<_edge> > _odom, _loops;

	public:

	void store(g2o::SparseOptimizer& optimizer)
	{
		_vertices.reserve(optimizer.vertices().size());

		for(size_t i=0; i<optimizer.vertices().size(); i++)
		{
			_vertices.push_back(dynamic_cast<g2o::VertexSE2*>(optimizer.vertex(i))->estimate());
		}

		g2o::OptimizableGraph::EdgeSet::iterator it = optimizer.edges().begin(),
				end = optimizer.edges().end();

		for( ; it!=end ; it++)
		{
			g2o::EdgeSE2* thisEdge = dynamic_cast<g2o::EdgeSE2*>(*it);
			_edge e;
			e.from = thisEdge->vertices()[0]->id();
			e.to = thisEdge->vertices()[1]->id();
			e.transform = thisEdge->measurement();
			e.information = thisEdge->information();

			std::cout<<e.from<<" > "<<e.to<<std::endl;

			if((e.to - e.from)==1)
			{
				_odom.push_back(e);
			}
			else
				_loops.push_back(e);
		}

		std::cout<<_vertices.size()<<" vertices stored"<<std::endl;
		std::cout<<_odom.size()<<" odom stored"<<std::endl;
		std::cout<<_loops.size()<<" loops stored"<<std::endl;

	}

	void get(SparseOptimizer& optimizer)
	{
		optimizer.clear();

		//std::vector<g2o::VertexSE2*> newVertices(_vertices.size(), NULL);
		std::vector<g2o::EdgeSE2*> newOdom(_odom.size(),NULL), newLoops(_loops.size(),NULL);

		std::cerr<<"Getting .. "<<std::endl;
		for(size_t i=0; i< _vertices.size(); i++)
		{
			g2o::VertexSE2* newVertex = new VertexSE2;
			newVertex->setId(i);
			newVertex->setEstimate(_vertices[i]);
			if(!optimizer.addVertex(newVertex)) std::cerr<<"Could not add vertex "<<i<<std::endl;
		}

		for(size_t i=0; i< _odom.size(); i++)
		{
			g2o::EdgeSE2* thisEdge = new EdgeSE2;
			thisEdge->vertices()[0] = dynamic_cast<VertexSE2*>(optimizer.vertices().find(_odom[i].from)->second);
			thisEdge->vertices()[1] = dynamic_cast<VertexSE2*>(optimizer.vertices().find(_odom[i].to)->second);
			thisEdge->setMeasurement(_odom[i].transform);
			thisEdge->setInformation(_odom[i].information);
			optimizer.addEdge(thisEdge);
		}


		for(size_t i=0; i< _loops.size(); i++)
		{
			g2o::EdgeSE2* thisEdge  = new EdgeSE2;
			thisEdge->vertices()[0] = dynamic_cast<VertexSE2*>(optimizer.vertices().find(_loops[i].from)->second);;
			thisEdge->vertices()[1] = dynamic_cast<VertexSE2*>(optimizer.vertices().find(_loops[i].to)->second);;
			thisEdge->measurement() = _loops[i].transform;
			thisEdge->information() = _loops[i].information;
			optimizer.addEdge(thisEdge);
		}

	}

};


class ClusterOptimizer
{
	typedef BlockSolver< BlockSolverTraits<-1, -1> >  SlamBlockSolver;
	typedef LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
	typedef std::vector< std::vector< g2o::EdgeSE2*> > EdgeMap;

	// allocating the optimizer
	SparseOptimizer optimizer, _backUp;
	SlamLinearSolver* linearSolver ; //= new SlamLinearSolver();
	SlamBlockSolver* solver ;
	unsigned int vertexCount;
	bool initialized;

	std::vector<g2o::OptimizableGraph::Edge*> _edgeSet, _odomSet;
	std::vector<bool> _edgeSetEnabled;

	std::vector<int> _outliers;

	std::vector<g2o::OptimizableGraph::Vertex*> _vertices;

	Clusterizer clusterizer;

	char * _filename;

	std::vector<int> _bestH;
	int _bestHPairings;
	double _bestChi2, previousChi2;
	// For graph "refreshing" ..

	GraphStorage graphStorage;

public:
	ClusterOptimizer()
	{
		linearSolver = new SlamLinearSolver();
		linearSolver->setBlockOrdering(false);
		solver = new SlamBlockSolver(&optimizer, linearSolver);
		optimizer.setSolver(solver);
		vertexCount = 0;
		initialized = false;
		_bestHPairings = 0;
		previousChi2 = 0;
	}
	bool write(std::string filename)
	{
		std::ofstream out(filename.c_str());
		if(!out.good()) std::cerr<<"Unable to open file "<<std::endl;

		for(size_t i=0; i<_vertices.size();i++)
		{
			out<<"VERTEX_SE2 "<<dynamic_cast<g2o::VertexSE2*>(_vertices[i])->id()<<" "; dynamic_cast<g2o::VertexSE2*>(_vertices[i])->write(out); out<<std::endl;
		}



		for(size_t i=0; i<_odomSet.size(); i++)
		{
			g2o::EdgeSE2* thisEdge = dynamic_cast<g2o::EdgeSE2*>(_odomSet[i]);
			out<<"EDGE_SE2 "<<thisEdge->vertices()[0]->id()<<" "<<thisEdge->vertices()[1]->id()<<" "; thisEdge->write(out); out<<std::endl;
		}



		for(size_t i=0; i<_edgeSet.size(); i++)
		{
			if(i<_edgeSetEnabled.size() and _edgeSetEnabled[i])
			{
				g2o::EdgeSE2* thisEdge = dynamic_cast<g2o::EdgeSE2*>(_edgeSet[i]);
				out<<"EDGE_SE2 "<<thisEdge->vertices()[0]->id()<<" "<<thisEdge->vertices()[1]->id()<<" "; thisEdge->write(out); out<<std::endl;
			}
		}



		//optimizer.save(out);
		if(out.good()) return true;
		return false;


	}

	bool reload()
	{
		_odomSet.clear(); _edgeSet.clear(); _vertices.clear();
		graphStorage.get(optimizer);

		int numberOfVertices = optimizer.vertices().size();
		if(numberOfVertices == 0)
		{
			std::cerr<<"No vertices in the graph "<<std::endl;
			return false;
		}

		for(int i=0 ; i<numberOfVertices ; i++)
			_vertices.insert(_vertices.end(),optimizer.vertex(i));

		g2o::OptimizableGraph::EdgeSet::iterator edgeIterator = optimizer.edges().begin(), end = optimizer.edges().end();


		for(; edgeIterator!=end; edgeIterator++)
		{
			g2o::EdgeSE2 * thisEdge = dynamic_cast<g2o::EdgeSE2*>(*edgeIterator);
			//thisEdge->setRobustKernel(true);
			//thisEdge->setHuberWidth(5.0);
			if(
					((thisEdge->vertices()[1])->id() - (thisEdge->vertices()[0])->id())
					==1
			)
			{
				//std::cout<<"O "<<thisEdge->vertices()[0]->id()<<" "<<thisEdge->vertices()[1]->id()<<std::endl;
				_odomSet.push_back(thisEdge);
			}
			else
			{
				//std::cout<<"L "<<thisEdge->vertices()[0]->id()<<" "<<thisEdge->vertices()[1]->id()<<std::endl;
				_edgeSet.push_back(thisEdge);
			}
			//_edgeSetEnabled = std::vector<bool>(_edgeSet.size(),true);
		}
		std::sort(_edgeSet.begin(), _edgeSet.end(), loopSort);

//		for(size_t i=0 ; i<_edgeSet.size(); i++)
//		{
//			std::cout<<_edgeSet[i]->vertices()[0]->id()<<" ";
//		}
//		std::cout<<std::endl;


	}

	static bool loopSort(g2o::OptimizableGraph::Edge* e1, g2o::OptimizableGraph::Edge* e2)
	{
		return (e1->vertices()[0]->id() < e2->vertices()[0]->id());
	}

	bool read(const char* filename)
	{
		_odomSet.clear(); _edgeSet.clear(); _vertices.clear();

		bool readOk = optimizer.load(filename,true);

		_filename = (char*)filename;

		if(!readOk)
		{
			std::cerr<<"Could not read "<<filename<<std::endl;
			return false;
		}

		graphStorage.store(optimizer);

		int numberOfVertices = optimizer.vertices().size();
		if(numberOfVertices == 0)
		{
			std::cerr<<"No vertices in the graph "<<std::endl;
			return false;
		}

		for(int i=0 ; i<numberOfVertices ; i++)
			_vertices.insert(_vertices.end(),optimizer.vertex(i));

		g2o::OptimizableGraph::EdgeSet::iterator edgeIterator = optimizer.edges().begin(), end = optimizer.edges().end();


		for(; edgeIterator!=end; edgeIterator++)
		{
			g2o::EdgeSE2 * thisEdge = dynamic_cast<g2o::EdgeSE2*>(*edgeIterator);
			//thisEdge->setRobustKernel(true);
			//thisEdge->setHuberWidth(5.0);
			if(
					((thisEdge->vertices()[1])->id() - (thisEdge->vertices()[0])->id())
					==1
			)
			{
				//std::cout<<"O "<<thisEdge->vertices()[0]->id()<<" "<<thisEdge->vertices()[1]->id()<<std::endl;
				_odomSet.push_back(thisEdge);
			}
			else
			{
				//std::cout<<"L "<<thisEdge->vertices()[0]->id()<<" "<<thisEdge->vertices()[1]->id()<<std::endl;
					_edgeSet.push_back(thisEdge);
			}
			//_edgeSetEnabled = std::vector<bool>(_edgeSet.size(),true);
		}
		std::sort(_edgeSet.begin(),_edgeSet.end(),loopSort);
		std::cerr<<"Read Done!"<<std::endl;
		return true;
	}

	int initializeCluster(int clusterID, const std::vector<int>& membership, g2o::OptimizableGraph::EdgeSet& activeEdges,
			std::vector<g2o::EdgeSE2*>& loopEdges)
	{

		int added = 0;
		//std::cout<<"ClusterInit()"<<std::endl;
		for(size_t i=0; i<membership.size();i++)
		{
			if(membership[i] == clusterID)
			{
				activeEdges.insert(_edgeSet[i]);
				loopEdges.push_back(dynamic_cast<EdgeSE2*>(_edgeSet[i]));
				//std::cout<<clusterID<<" "<<_edgeSet[i]->vertices()[0]->id()<<" "<<_edgeSet[i]->vertices()[1]->id()<<std::endl;
				added++;
			}
		}
		//std::cout<<"/ClusterInit()"<<std::endl;
		return added;
	}

	void refineOnResiduals(const std::vector<g2o::EdgeSE2*>& loopEdges, std::vector<bool>& isOkay)
	{
		std::vector<g2o::EdgeSE2*>::const_iterator it = loopEdges.begin() , end = loopEdges.end();

		for( ; it!=end ; it++)
		{
			(*it)->computeError();
			if( (*it)->chi2()<5.99){ std::cout<< "[ ]"; isOkay.push_back(true); }
			else { std::cout<<"[X]"; isOkay.push_back(false); }
		}
		std::cout<<std::endl;
	}

	void IC(int clusterID, std::vector<int>& membership, const int iterations)
	{
		//optimizer.clear();
		//read(_filename);

		reload();

		g2o::OptimizableGraph::EdgeSet activeEdges;
		activeEdges.insert(_odomSet.begin(),_odomSet.end());

		std::vector<g2o::EdgeSE2*> loopEdges;

		initializeCluster(clusterID,membership,activeEdges, loopEdges);

		//activeEdges.insert(edgeMap[clusterID].begin(), edgeMap[clusterID].end());
		//loopEdges.insert(loopEdges.end(),edgeMap[clusterID].begin(), edgeMap[clusterID].end());

		optimizer.setVerbose(false);
		optimizer.vertex(0)->setFixed(true);
		optimizer.initializeOptimization(activeEdges);
		optimizer.optimize(iterations);
		optimizer.computeActiveErrors();

		std::vector<g2o::EdgeSE2*>::const_iterator it = loopEdges.begin() , end = loopEdges.end();

		std::vector<bool> goodOnes(loopEdges.size(),false);

		double chiSquare = 0;
		int correctOnes = 0;

		for(int i=0 ; it!=end ; it++, i++)
		{
			(*it)->computeError();
			//std::cout<<(*it)->chi2()<<" ";
			if( (*it)->chi2()<5.99)
			{
				std::cout<< "[ ]";
				goodOnes[i]=true;
				correctOnes ++;
				chiSquare+= (*it)->chi2();
			}
			else std::cout<<"[X]";
		}
		std::cout<<std::endl;

		if(correctOnes > 0)
		{
			double threshold = boost::math::quantile(boost::math::chi_squared(3*correctOnes-1),0.95);
			double threshold2 = boost::math::quantile(boost::math::chi_squared(3*optimizer.activeEdges().size()-1),0.95);
			std::cout<<" Chi2 "<<chiSquare<<"/"<<threshold<<" "<<optimizer.activeChi2()<<"/"<<threshold2<<std::endl;
		}


		for(size_t i=0, j=0; i< membership.size() ; i++)
		{
			if(membership[i]==clusterID)
			{
				if(goodOnes[j]==false)
				{
					membership[i]= -2;
				}
				j++;
			}
		}

	}

	bool JC( std::vector<int>& H, const std::vector<int>& membership, int iterations)//, std::vector<int>& rejected, int checkLast = 0)
	{

		if(H.empty())
		{
			std::cout<<" H is empty! "<<std::endl;
			return true;
		}

		if(H.back()==-1)
			return true;

		//optimizer.clear();
		//read(_filename);

		reload();

		g2o::OptimizableGraph::EdgeSet activeEdges;
		activeEdges.insert(_odomSet.begin(),_odomSet.end());

		std::vector<g2o::EdgeSE2*> loopEdges;

		std::vector<int> membershipNow;

		std::vector<int> count(H.size(),0);

		std::cout<<"JC : [ ";
		for(size_t i=0; i<H.size(); i++)
		{
			std::cout<<H[i]<<" ";
			if(H[i]< 0 ) continue;
			int added = initializeCluster(H[i],membership,activeEdges,loopEdges);
			count[i]=added;
			for(int j=0; j<added ; j++)
			{
				membershipNow.push_back(H[i]);
			}
		}
		std::cout<<std::endl;

		optimizer.setVerbose(false);
		optimizer.vertex(0)->setFixed(true);

		optimizer.initializeOptimization(activeEdges);
		optimizer.optimize(iterations);
		optimizer.computeActiveErrors();

		std::vector<g2o::EdgeSE2*>::const_iterator it = loopEdges.begin() , end = loopEdges.end();

		double chiSquare = 0;

		for( int i=0 ; it!=end ; it++, i++)
		{
			(*it)->computeError();
			chiSquare += (*it)->chi2();
		}

		previousChi2 = chiSquare;

		double allThreshold =boost::math::quantile(boost::math::chi_squared(3* optimizer.edges().size()-1),0.95);
		double threshold = boost::math::quantile(boost::math::chi_squared(3*loopEdges.size()-1),0.95);
		std::cout<<"ChiSqaure "<<chiSquare<<"/"<<threshold<<" "
				<<optimizer.activeChi2()<<"/"<<allThreshold<<std::endl;

		if (chiSquare < threshold and optimizer.activeChi2()< allThreshold) //  // and
			return true;
		else
		{
			return false;
		}

	}

	void JointSelectionClusters(std::vector<int>& H, const std::vector<int>& membership, int iterations, int numClusters, int depth = 0)
	{
		std::cout<<" --- ";
		for(size_t i=0 ; i< H.size() ; i++)
		{
			std::cout<<H[i]<<" ";
		}
		std::cout<<std::endl;

		if(depth >= numClusters)
		{
			std::cout<<" Leaf node "<<std::endl;
			int currentPairing = 0;
			for(size_t i=0 ; i< H.size() ; i++)
			{
				std::cout<<H[i]<<" ";
				if(H[i]>-1) currentPairing++;
			}

			if( (currentPairing > _bestHPairings)
					or
					(currentPairing == _bestHPairings and _bestChi2 > previousChi2))
			{
				_bestH = H;
				_bestHPairings = currentPairing;
				_bestChi2 = previousChi2;
			}

			std::cout<<std::endl;
			std::cerr<<"Done leaf noce"<<std::endl;



			return;
		}

		std::vector<int> _H = H;
		_H.push_back(depth);
		if(JC(_H,membership,iterations))
		{
			JointSelectionClusters(_H,membership,iterations,numClusters,depth+1);
		}

		int currentPairing = 0;

		for(size_t i=0; i<H.size(); i++)
		{
			if(H[i]> 0) currentPairing++;
		}

		_H = H;
		if(currentPairing + (numClusters-depth-1) > _bestHPairings )
		{
			_H.push_back(-1);
			JointSelectionClusters(_H,membership,iterations,numClusters,depth+1);
		}
		return;
	}


	void display(std::set<int>&s )
	{
		std::set<int>::iterator it = s.begin(), end = s.end();
		for( ; it!=end ; it++)
		{
			std::cout<<*it<<" ";
		}
		std::cout<<std::endl;
	}

	void run(int iterations, int clusterThreshold)
	{
		std::cerr<<"run"<<std::endl;
		std::vector<int> loops;

		for(size_t i =0 ; i<_edgeSet.size(); i++)
		{
			std::cout<<_edgeSet[i]->vertices()[0]->id()<<" -> "<<_edgeSet[i]->vertices()[1]->id()<<std::endl;
			loops.push_back(_edgeSet[i]->vertices()[0]->id());
			loops.push_back(_edgeSet[i]->vertices()[1]->id());
		}

		std::cerr<<"clusterize "<<std::endl;
		// Form clusters //
		std::vector<int> membership;
		int clusterCount;
		clusterizer.clusterize(loops,clusterThreshold, membership, clusterCount);

		for(size_t i =0 ; i<_edgeSet.size(); i++)
		{
			std::cout<<membership[i]<<" ";
			std::cout<<_edgeSet[i]->vertices()[0]->id()<<" -> "<<_edgeSet[i]->vertices()[1]->id()<<std::endl;
			loops.push_back(_edgeSet[i]->vertices()[0]->id());
			loops.push_back(_edgeSet[i]->vertices()[1]->id());
		}

		std::cout<<" done clusters "<<std::endl;

		_edgeSetEnabled = std::vector<bool>(membership.size(),false);

		for(int i=0; i< clusterCount; i++)
		{
			IC(i,membership,iterations);

		}

		for(size_t i=0 ; i<membership.size(); i++)
		{
			std::cout<<membership[i]<<" ";
		}
		std::cout<<std::endl;

		std::set<int> goodClusters,			// Clusters selected in every iteration
					  selectedClusters, 	// Overall set of JC sets that are selected
					  rejectedClusters;		// Clusters that have been rejected (NOT JC) in the prev iteration.

		g2o::OptimizableGraph::EdgeSet activeEdges;
		std::vector<g2o::EdgeSE2*> loopEdges;

		std::vector<int> H;
		_bestChi2 = 1e10;

		optimizer.setVerbose(false);

		JointSelectionClusters(H,membership,iterations,clusterCount,0);

		std::cerr<<"Return .. "<<std::endl;

		_edgeSetEnabled = std::vector<bool> (membership.size(),false);

		ofstream out("memebership.txt");

		for(size_t j=0; j<membership.size(); j++)
		{
			for(size_t i=0; i<_bestH.size(); i++)
			{
				if(_bestH[i]==-1) continue;
				if(membership[j]==_bestH[i]) _edgeSetEnabled[j] = true;
			}
		}

		reload();

		write("second.g2o");

		for(size_t j=0; j<_edgeSetEnabled.size(); j++)
		{
			if(_edgeSetEnabled[j]) out<<1<<" ";
			else out<<0<<" ";
		}

		out.close();

		return;
		/*** JCBB *

		_edgeSetEnabled = std::vector<bool>(membership.size());
		for(int i=0; i<membership.size();i++)
		{
			_edgeSetEnabled[i] = (membership[i]>-1);
		}
		write("IC.g2o");



		std::vector<int> H;

		_bestH.clear();
		_bestHPairings = 0;

		JointSelectionClusters(H,membership,iterations,clusterCount,0);

		optimizer.clear();
		read(_filename);

		_edgeSetEnabled = std::vector<bool>(_edgeSet.size(),false);

		std::cerr<<_edgeSet.size()<<" "<<membership.size()<<std::endl;
		for(int i=0; i<membership.size(); i++)
		{
			for(int j=0; j< _bestH.size(); j++)
			{
				if(membership[i]==_bestH[j])
				{
					_edgeSetEnabled[i] = true;
					break;
				}
			}
		}

	/*** JCBB */


	}


};
#endif /* OPTCLUSTER_HPP_ */
