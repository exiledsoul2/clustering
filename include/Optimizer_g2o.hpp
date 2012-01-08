
#include <g2o/core/graph_optimizer_sparse.h>
#include <g2o/core/block_solver.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/types/slam2d/edge_se2.h>
#include <g2o/types/slam2d/vertex_se2.h>

#include <fstream>

#include <boost/math/distributions/chi_squared.hpp>

#include <Cluster.hpp>

using namespace g2o;
class Optimizer
{
	typedef BlockSolver< BlockSolverTraits<-1, -1> >  SlamBlockSolver;
	typedef LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

  // allocating the optimizer
	SparseOptimizer optimizer;
	SlamLinearSolver* linearSolver ; //= new SlamLinearSolver();
	SlamBlockSolver* solver ; 
	unsigned int vertexCount;
	bool initialized;

	std::vector<g2o::OptimizableGraph::Edge*> _edgeSet, _odomSet;
	std::vector<bool> _edgeSetEnabled;

	std::vector<int> _outliers;

	std::vector<g2o::OptimizableGraph::Vertex*> _vertices;

	Clusterizer clusterizer;

	int _selected;

	char * _filename;
	public:
	Optimizer()
	{
		linearSolver = new SlamLinearSolver();
		linearSolver->setBlockOrdering(false);
		solver = new SlamBlockSolver(&optimizer, linearSolver);
		optimizer.setSolver(solver);
		vertexCount = 0;
		initialized = false;
		_selected = 0;
	}
	bool addVertex(const Eigen::Vector3d& v)
	{
		VertexSE2* vertex = new VertexSE2;
		vertex->setId(vertexCount);
		vertex->setEstimate(SE2(v(0),v(1),v(2)));
		optimizer.addVertex(vertex);
		_vertices.push_back(vertex);
		vertexCount++;
		return true;
	}
	bool addEdgeLoop(
		unsigned int from,
		unsigned int to,
		const Eigen::Vector3d& t, 
		const Eigen::Matrix3d& information)
	{
		SE2 transf(t(0),t(1),t(2));
		EdgeSE2* e = new EdgeSE2;
		e->setMeasurement(transf);
		e->setInverseMeasurement(transf.inverse());
		e->setInformation(information);
		e->vertices()[0] = optimizer.vertex(from);
		e->vertices()[1] = optimizer.vertex(to);

		optimizer.addEdge(e);

		_edgeSet.push_back(e);
		_edgeSetEnabled.push_back(true);

		return true;
	}
	
	bool addEdgeOdometry(unsigned int from, unsigned int to, const Eigen::Matrix3d& information)
	{
		EdgeSE2* e = new EdgeSE2;
		VertexSE2* fromV = dynamic_cast<VertexSE2*>(optimizer.vertices().find(from)->second);
		VertexSE2* toV = dynamic_cast<VertexSE2*>(optimizer.vertices().find(to)->second);
		
		SE2 transf = fromV->estimate().inverse() * toV->estimate();
		
		e->setMeasurement(transf);
		e->setInverseMeasurement(transf.inverse());
		e->setInformation(information);
		e->vertices()[0] = optimizer.vertex(from);
		e->vertices()[1] = optimizer.vertex(to);
		
		optimizer.addEdge(e);
		_odomSet.push_back(e);

		return true;
	}
	
	bool initializeEdges(int n = -1)
	{
		//for(size_t i=0; i<_vertices.size();i++) optimizer.addVertex(_vertices[i]);
		//for(g2o::OptimizableGraph::EdgeSet::iterator it = _odomSet.begin(); it!=_odomSet.end() ; it++)
		//		optimizer.addEdge(dynamic_cast<g2o::EdgeSE2*>(*it));
		//for(std::vector<g2o::OptimizableGraph::Edge*>::iterator it = _edgeSet.begin(); it!=_edgeSet.end() ; it++)
		//		optimizer.addEdge(dynamic_cast<g2o::EdgeSE2*>(*it));

		optimizer.clear();
		read(_filename);

		//std::cerr<<"Done "<<std::endl;

		g2o::OptimizableGraph::EdgeSet _edges;
		_edges.insert(_odomSet.begin(),_odomSet.end());

		if(n<0)
		{
			for(size_t i=0; i<_edgeSet.size(); i++)
			{
				if(_edgeSetEnabled[i])
					_edges.insert(_edgeSet[i]);
			}
		}
		else
		{
			_edges.insert(_edgeSet[n]);
		}


		optimizer.vertex(0)->setFixed(true);
		optimizer.initializeOptimization(_edges);
		initialized = true;
		return true;
	}
	
	bool initializeEdges(std::vector<int>& edgeList)
	{
		g2o::OptimizableGraph::EdgeSet _edges;
		_edges.insert(_odomSet.begin(),_odomSet.end());

		for(size_t i=0; i<edgeList.size();i++)
		{
			_edges.insert(_edgeSet[edgeList[i]]);
		}

		optimizer.vertex(0)->setFixed(true);
		optimizer.initializeOptimization(_edges);
		return true;
	}

	bool initializeEdgesAll() // include everything even those that are rejected
	{
		g2o::OptimizableGraph::EdgeSet _edges;
		_edges.insert(_odomSet.begin(),_odomSet.end());
		_edges.insert(_edgeSet.begin(),_edgeSet.end());

		optimizer.vertex(0)->setFixed(true);
		//std::cout<<"Graph is fixed by vertex 0"<<std::endl;
		optimizer.initializeOptimization(_edges);
		return true;
	}

	bool initializeCluster(int clusterID, const std::vector<int>& membership)
	{
		optimizer.clear();
		read(_filename);

		g2o::OptimizableGraph::EdgeSet _edges;
		_edges.insert(_odomSet.begin(),_odomSet.end());

		for(size_t i=0; i<membership.size();i++)
		{
			if(membership[i] == clusterID)
				_edges.insert(_edgeSet[i]);
		}
		optimizer.vertex(0)->setFixed(true);
		optimizer.initializeOptimization(_edges);
		return true;
	}

	double activeChi2()
	{
		return optimizer.activeChi2();
	}
	
	double optimize(int iterations = 1)
	{
		//if(!initialized)
		optimizer.optimize(iterations);
		return optimizer.activeChi2();
	}
	
	void oneFunction(int iterations)
	{
		std::vector<int> loops;
		for(size_t i =0 ; i<_edgeSet.size(); i++)
		{
			std::cout<<_edgeSet[i]->vertices()[0]->id()<< " -> "<<_edgeSet[i]->vertices()[1]->id() <<std::endl;
			loops.push_back(_edgeSet[i]->vertices()[0]->id());
			loops.push_back(_edgeSet[i]->vertices()[1]->id());
		}

		// Form clusters //
		std::vector<int> membership;
		int clusterCount;
		clusterizer.clusterize(loops,50, membership, clusterCount);

		for(int i=0 ; i<clusterCount;i++)
		{
			srand(time(0));

			int randomClusterID = rand()%clusterCount;
			randomClusterID = i;

			std::cerr<<" Cluster Selected : "<<randomClusterID<<std::endl;

			initializeCluster(i, membership);

			optimizer.optimize(iterations);
			optimizer.computeActiveErrors();
			refineOnResidualErrors();

			for(size_t j=0; j< _edgeSetEnabled.size() ; j++)
			{
				if(!_edgeSetEnabled[j] and membership[j]==randomClusterID)
					membership[j] = -1;
			}

		}



		//refineOnResidualErrors((_selected+1)*3 -1,true);

		write("firstGroupDisabled.g2o");

		// Now _edgeListEnable Contains good links.

		std::vector<bool> goodSet = _edgeSetEnabled;

		for(size_t i=0; i<_edgeSetEnabled.size(); i ++ )
		{
			_edgeSetEnabled[i]=!goodSet[i];
		}

		bool done = false;

		while(!done){

			initializeEdges();
			optimizer.optimize(iterations);
			done = !refineOnResidualErrors(2,true);

			for(size_t i=0; i<_edgeSetEnabled.size(); i ++ )
			{
				goodSet[i] = goodSet[i] or _edgeSetEnabled[i];
			}

			for(size_t i=0; i<_edgeSetEnabled.size(); i ++ )
			{
				_edgeSetEnabled[i]=!goodSet[i];
			}
		}
		write("secondGroupDisabled.g2o");


		std::cout<<" ----------------------------------------------- "<<std::endl;

		for(size_t i=0 ; i< goodSet.size() ; i++)
		{
			(goodSet[i])? std::cout<<"[ ]":std::cout<<"[X]";
		}
		std::cout<<std::endl;

		std::cin.get();

		_edgeSetEnabled = goodSet;

		_selected = 1;
		int _previousSelected = 0;
		//for(int i=0 ; i<10 ; i++)
		int i=0;
		while(_previousSelected != _selected)
		{
			initializeEdges();

			char fname[256];;
			sprintf(fname,"output/Opt-%02d.g2o",i);
			write(fname);

			optimizer.optimize(3);
			std::cout<<"After second optimization " << optimizer.activeChi2()<<std::endl;
			std::cout<<_selected<<" links selected"<<std::endl;
			_previousSelected = _selected;
			int p = _selected ;
			do
			{
				p = _selected;
				refineOnResidualErrors((_selected+1)*3 -1);
				std::cout<<_selected<<" after checking links selected"<<std::endl;
			}while(_selected>p);
			std::cout<<"Run :"<<++i<<std::endl;
		}





	}

	void jointCompatibility(std::vector<int> includedClusterIDs)
	{


	}
	//bool getNormalizedResiduals(std::vector<double>& res)
	bool refineOnResidualErrors(double dof=2, bool value = true)
	{
		std::vector<g2o::OptimizableGraph::Edge*>::iterator it;
		
		_edgeSetEnabled = std::vector<bool>(_edgeSet.size(),!value);
		_selected = 0;

		double threshold = boost::math::quantile(boost::math::chi_squared(dof),0.95);
		std::cout<<" Chi : dof : "<<dof<<" threshold : "<<threshold<<std::endl;

		double residualSum = 0;

		for(size_t i=0; i< _edgeSet.size(); i++)
		{
			EdgeSE2* thisEdge = dynamic_cast<EdgeSE2*>(_edgeSet[i]);

			thisEdge->computeError();
			double normalizedResidual = thisEdge->chi2();

			(normalizedResidual<threshold)?std::cout<<"[ ]":std::cout<<"[X]";

			if(normalizedResidual < threshold)
			{
				residualSum += normalizedResidual;
				_edgeSetEnabled[i] = value;
				_selected ++ ;
			}
		}

		if(_selected == 0)
			return false;

		double chiSquareBound = boost::math::quantile(boost::math::chi_squared(3*_selected-1),0.95);
		std::cout<<std::endl;
		std::cout<<"Sum of Residuals : "<<residualSum<<" (bound : "<<chiSquareBound<<") ";
		if(residualSum < chiSquareBound)
		{
			std::cout<<" (in bound) "<<std::endl;
			return true;
		}
		else
		{
			std::cout<<" (out of bound) "<<std::endl;
			return false;
		}
		return true;
	}
	
	bool write(std::string filename)
	{
		std::ofstream out(filename.c_str());
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
		return true;
	}


};
