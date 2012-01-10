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


class ClusterManager{
	typedef g2o::OptimizableGraph::Edge* EdgePtrType;
	typedef std::pair< EdgePtrType, int > EdgePairType;
	std::map<EdgePtrType, int > _membershipLookUp;
	std::vector< std::vector< EdgePtrType > > _clusters;

public:
	ClusterManager(){}
	void init(std::vector<int>& membership, std::vector<EdgePtrType>& edges, const int& clusterCount)
	{
		//_clusters.reserve(clusterCount);
		for(size_t i=0, runTill = membership.size(); i<runTill;  i++)
		{
			if( membership[i] > (int)(_clusters.size()-1))
			{
				_clusters.push_back(std::vector<EdgePtrType>());
			}
			if(membership[i]>=0)
			{
				_clusters[membership[i]].push_back(edges[i]);
				_membershipLookUp.insert(EdgePairType(edges[i],membership[i]));
			}
		}
	}

	void update(std::vector<int>& membership, std::vector<EdgePtrType>& edges, const int& clusterCount)
	{
		for(size_t i=0 ; i< _clusters.size(); i++)
		{
			_clusters[i].clear();
		}
		_membershipLookUp.clear();
		init(membership,edges,clusterCount);
	}

	void getCluster(int clusterID, g2o::OptimizableGraph::EdgeSet& cluster)
	{
		cluster.insert(_clusters[clusterID].begin(), _clusters[clusterID].end());
	}

	int getCluster(int clusterID, std::vector<g2o::EdgeSE2*>& cluster)
	{
		std::vector<EdgePtrType>::iterator it = _clusters[clusterID].begin(), end = _clusters[clusterID].end();
		for ( ; it!=end ; it++)
			cluster.push_back(dynamic_cast<g2o::EdgeSE2*>(*it));
		return (int)_clusters[clusterID].size();
	}

	int getMembership(const EdgePtrType& e)
	{
		if(_membershipLookUp.find(e)!=_membershipLookUp.end())
			return _membershipLookUp.find(e)->second;
		else
			return -2;
	}
};


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
	std::vector< _edge , Eigen::aligned_allocator<_edge> > _edges;

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

			_edges.push_back(e);

		}

		std::cout<<_vertices.size()<<" vertices stored"<<std::endl;
		std::cout<<_edges.size()<<" odom stored"<<std::endl;

	}

	void get(SparseOptimizer& optimizer)
	{
		//optimizer.clear();

		//std::vector<g2o::VertexSE2*> newVertices(_vertices.size(), NULL);
		//std::vector<g2o::EdgeSE2*> newOdom(_odom.size(),NULL), newLoops(_loops.size(),NULL);

		std::cerr<<"Getting .. "<<std::endl;

		for(size_t i=0, verticesSize = _vertices.size() ; i< verticesSize ; i++)
		{
			dynamic_cast<g2o::VertexSE2*>(optimizer.vertex(i))->setEstimate(_vertices[i]);
			/*
			g2o::VertexSE2* newVertex = new VertexSE2;
			newVertex->setId(i);
			newVertex->setEstimate(_vertices[i]);
			if(!optimizer.addVertex(newVertex)) std::cerr<<"Could not add vertex "<<i<<std::endl;
			*/
		}

		g2o::OptimizableGraph::EdgeSet::iterator it = optimizer.edges().begin(), end = optimizer.edges().end();

		for(size_t i=0 ; it!=end ; i++ , it++)
		{
			g2o::EdgeSE2* thisEdge = dynamic_cast<g2o::EdgeSE2*>(*it);
			//g2o::EdgeSE2* thisEdge = new EdgeSE2;
			thisEdge->vertices()[0] = dynamic_cast<VertexSE2*>(optimizer.vertices().find(_edges[i].from)->second);
			thisEdge->vertices()[1] = dynamic_cast<VertexSE2*>(optimizer.vertices().find(_edges[i].to)->second);
			thisEdge->setMeasurement(_edges[i].transform);
			thisEdge->setInformation(_edges[i].information);
			//optimizer.addEdge(thisEdge);
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

	// For graph "refreshing" ..

	GraphStorage graphStorage;
	ClusterManager clusterManager;

public:
	ClusterOptimizer()
	{
		linearSolver = new SlamLinearSolver();
		linearSolver->setBlockOrdering(false);
		solver = new SlamBlockSolver(&optimizer, linearSolver);
		optimizer.setSolver(solver);
		vertexCount = 0;
		initialized = false;
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
		//_odomSet.clear(); _edgeSet.clear(); _vertices.clear();
		graphStorage.get(optimizer);

		return true;
		/*

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
	*/
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

		_vertices.reserve(optimizer.vertices().size());
		_odomSet.reserve(optimizer.edges().size());
		_edgeSet.reserve(optimizer.edges().size());

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
					std::abs((thisEdge->vertices()[1])->id() - (thisEdge->vertices()[0])->id())
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
		//std::cout<<"ClusterInit()"<<std::endl;
		clusterManager.getCluster(clusterID,activeEdges);
		int added = clusterManager.getCluster(clusterID, loopEdges);
		return added;
	}

	void refineOnResiduals(const std::vector<g2o::EdgeSE2*>& loopEdges, std::vector<bool>& isOkay)
	{
		isOkay.reserve(loopEdges.size());
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


		for(size_t i=0, j=0 , runTill = membership.size(); i< runTill ; i++)
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

	bool JC( std::vector<int>& H, const std::vector<int>& membership, int iterations, std::vector<int>& rejected, int checkLast = 0)
	{

		if(H.empty() or checkLast == 0)
		{
			std::cout<<" H is empty! "<<std::endl;
			return true;
		}

		reload();

		g2o::OptimizableGraph::EdgeSet activeEdges;
		activeEdges.insert(_odomSet.begin(),_odomSet.end());

		std::vector<g2o::EdgeSE2*> loopEdges;

		std::vector<int> membershipNow;

		std::vector<int> count(H.size(),0);

		for(size_t i=0 , runTill = H.size(); i< runTill; i++)
		{
			//std::cout<<"Adding Cluster with iD "<<H[i]<<std::endl;
			int added = initializeCluster(H[i],membership,activeEdges,loopEdges);
			count[i]=added;
			for(int j=0; j<added ; j++)
			{
				membershipNow.push_back(H[i]);
			}
		}

		//std::cout<<"# Loops "<<loopEdges.size()<<std::endl;
		//for(std::vector<g2o::EdgeSE2*>::iterator it = loopEdges.begin(); it!=loopEdges.end(); it++)
		//{
		//	std::cout<<(*it)->vertices()[0]->id()<<" "<<(*it)->vertices()[1]->id()<<std::endl;
		//}

		optimizer.setVerbose(false);
		optimizer.vertex(0)->setFixed(true);

		optimizer.initializeOptimization(activeEdges);
		optimizer.optimize(iterations);
		optimizer.computeActiveErrors();

		std::map<int,int> reverse_map;

		for(size_t i=0 , runTill = H.size() ; i< runTill; i++)
		{
			reverse_map.insert(std::pair<int,int>(H[i],i));
		}

		std::vector<double> errors(H.size(),0);


		double chiSquare = 0;

		double maxLinkError = 0;
		int maxLinkErrorID = -1;

		std::vector<g2o::EdgeSE2*>::const_iterator it = loopEdges.begin() , end = loopEdges.end();

		for( int i=0 ; it!=end ; it++, i++)
		{
			(*it)->computeError();
			chiSquare += (*it)->chi2();
			//std::cout<<"["<<membershipNow[i]<<" "<<(*it)->chi2()<<" ]"<<std::endl;;
			if ( (*it)->chi2()>maxLinkError )
			{
				maxLinkError = (*it)->chi2();
				maxLinkErrorID = reverse_map[membershipNow[i]];
			}
			errors[reverse_map[membershipNow[i]]] +=(*it)->chi2();
		}
		//std::cout<<std::endl;

		//std::cout<<"Avg. Residual per cluster .. "<<std::endl;
		for(size_t i=0; i<H.size(); i++)
		{
			std::cout<<H[i]<<" ";
			std::cout<<errors[i]/count[i]<<std::endl;
		}

		double allThreshold =boost::math::quantile(boost::math::chi_squared(3* optimizer.edges().size()-1),0.95);
		double threshold = boost::math::quantile(boost::math::chi_squared(3*loopEdges.size()-1),0.95);
		std::cout<<"ChiSqaure "<<chiSquare<<"/"<<threshold<<" "
				<<optimizer.activeChi2()<<"/"<<allThreshold<<std::endl;

		if (chiSquare < threshold and optimizer.activeChi2()< allThreshold) //  // and
			return true;
		else
		{
			double maxAvg = 0; int maxIndex = -1;
			for(size_t i= H.size() - checkLast ; i<H.size(); i++)
			//for(size_t i= 0 ; i<H.size(); i++)
			{
				if(errors[i] > maxAvg)//if(errors[i]/count[i] > maxAvg)
				{
					maxAvg = errors[i];//maxAvg = errors[i]/count[i];
					maxIndex = i;
				}
			}

			rejected.push_back(H[maxIndex]);
			H.erase(H.begin()+maxIndex);
			//rejected.push_back(H[maxLinkErrorID]);
			//H.erase(H.begin()+maxLinkErrorID);

			return JC(H,membership,iterations,rejected, checkLast-1);
		}

	}
/*
	void JointSelectionClusters(std::vector<int>& H, const std::vector<int>& membership, int iterations, int numClusters, int depth = 0)
	{
		std::cout<<" --- ";
		for(size_t i=0 ; i< H.size() ; i++)
		{
			std::cout<<H[i]<<" ";
		}
		std::cout<<std::endl;

		if(depth >= numClusters-1)
		{
			std::cout<<" Leaf node "<<std::endl;
			int currentPairing = 0;
			for(size_t i=0 ; i< H.size() ; i++)
			{
				std::cout<<H[i]<<" ";
				if(H[i]>-1) currentPairing++;
			}
			if(currentPairing > _bestHPairings)
			{
				_bestH = H;
				_bestHPairings = currentPairing;
			}
			std::cout<<std::endl;

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
	*/

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
		std::ofstream mFile("membership.txt");

		std::cerr<<"run"<<std::endl;
		std::vector<int> loops;

		for(size_t i =0 , runTill = _edgeSet.size() ; i<runTill; i++)
		{
			std::cout<<_edgeSet[i]->vertices()[0]->id()<<" -> "<<_edgeSet[i]->vertices()[1]->id()<<std::endl;
			loops.push_back(_edgeSet[i]->vertices()[0]->id());
			loops.push_back(_edgeSet[i]->vertices()[1]->id());
		}

		// Form clusters //
		std::vector<int> membership;
		int clusterCount;
		clusterizer.clusterize(loops,clusterThreshold, membership, clusterCount);


#if 1
		for(size_t i =0 ; i<_edgeSet.size(); i++)
		{
			mFile<<membership[i]<<" ";
		}
		mFile<<std::endl;
#endif
		std::cout<<" done clusters "<<std::endl;

		clusterManager.init(membership,_edgeSet,clusterCount);

		for(int i=0; i< clusterCount; i++)
		{
			IC(i,membership,iterations);

		}
		clusterManager.update(membership,_edgeSet,clusterCount);

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

		std::vector<int> H(clusterCount);

		bool done = false;

		optimizer.setVerbose(false);

		while(!done)
		{
			reload();

			activeEdges.clear();
			loopEdges.clear();

			activeEdges.insert(_odomSet.begin(), _odomSet.end());

			for(int i=0 ; i<clusterCount; i++)
			{
				if(selectedClusters.find(i)==selectedClusters.end()){ // Not in the selected Cluster
					initializeCluster(i,membership,activeEdges,loopEdges);
				}
			}

			optimizer.vertex(0)->setFixed(true);
			optimizer.initializeOptimization(activeEdges);
			optimizer.optimize(iterations);
			optimizer.computeActiveErrors();

			done = true;
			for(size_t i=0; i<_edgeSet.size(); i++)
			{
				_edgeSet[i]->computeError();
				if(	membership[i]>=0
					and selectedClusters.find(membership[i])==selectedClusters.end()
					and rejectedClusters.find(membership[i])==rejectedClusters.end()
					and	_edgeSet[i]->chi2() < 5.99 )
				{
					//std::cout<<"["<<membership[i]<<"]";
					goodClusters.insert(membership[i]);
					done = false;
				}
				//else{	std::cout<<"[X]";}
			}


			//reload();

			H.clear();
			H.insert(H.end(),selectedClusters.begin(), selectedClusters.end());
			H.insert(H.end(), goodClusters.begin(), goodClusters.end());

			std::vector<int> rejected;

			int oldSize = selectedClusters.size();

			JC(H,membership,iterations, rejected, goodClusters.size());

			goodClusters.clear();
			selectedClusters.clear();

			selectedClusters.insert(H.begin(),H.end());

			if(selectedClusters.size() > oldSize)
			{
				rejectedClusters.clear();
			}
			rejectedClusters.insert(rejected.begin(), rejected.end());

			display(selectedClusters);
			std::cout<<"Rejected List" << std::endl ; display(rejectedClusters);

		}


		_edgeSetEnabled = std::vector<bool>(membership.size(),false);
		for(size_t i=0 ; i<membership.size() ; i++)
		{
			if(selectedClusters.find(membership[i])!=selectedClusters.end()) _edgeSetEnabled[i]=true;
			else _edgeSetEnabled[i]=false;
		}

		reload();
		write("second.g2o");


		H.clear();
		H.insert(H.end(),selectedClusters.begin(), selectedClusters.end());

		activeEdges.clear();
		loopEdges.clear();

		activeEdges.insert(_odomSet.begin(), _odomSet.end());

		for(int i=0 ; i<clusterCount; i++)
		{
			if(selectedClusters.find(i)!=selectedClusters.end()){
				initializeCluster(i,membership,activeEdges,loopEdges);
			}
		}

		optimizer.vertex(0)->setFixed(true);
		optimizer.initializeOptimization(activeEdges);
		optimizer.optimize(iterations*2);
		optimizer.computeActiveErrors();

		for(size_t i=0; i<_edgeSet.size(); i++)
		{
			_edgeSet[i]->computeError();
			if(	_edgeSet[i]->chi2() < 5.99 )
			{
				std::cout<<"["<<membership[i]<<"]";
			}
			else{	std::cout<<"[X]";}
		}

		reload();
		write("check.g2o");

		for(size_t i=0 ; i<membership.size() ; i++)
		{
			if(selectedClusters.find(membership[i])!=selectedClusters.end())
				mFile<<1<<" ";
			else
				mFile<<0<<" ";
		}
		mFile<<std::endl;

		mFile.close();

		return;

	}


};
#endif /* OPTCLUSTER_HPP_ */
