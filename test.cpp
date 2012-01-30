#include "Simulator/Simulator.hpp"

//#include <Optimizer_g2o.hpp>
//#include <OptCluster.hpp>
//#include <ClusterJCBB.hpp>
#include<PersisCluster.hpp>

//#include<ClusterRJC.hpp>

using namespace Simulator;
#define STEP_LENGTH 1.0

//#define MINIMIZE_ORIGINAL
#if 0
void lawn_mower_sim(ClusterOptimizer& optimizer)
{
	Simulator2D sim;

			std::cerr<<"Setting Noise models ...";

			sim.OdometryNoise()[0]= (0.1*STEP_LENGTH)*(0.1*STEP_LENGTH);
			sim.OdometryNoise()[1]= (0.1*STEP_LENGTH)*(0.1*STEP_LENGTH);
			sim.OdometryNoise()[2]= 0.0174532925 * 0.0174532925;

			std::cerr<<"Setting motion ...";
			MotionPlan& plan = sim.motionPlan();
			plan.stepLength() = STEP_LENGTH;

			sim.OdometryNoise()[0]= (0.01*STEP_LENGTH)*(0.01*STEP_LENGTH) ;
			sim.OdometryNoise()[1]= (0.01*STEP_LENGTH)*(0.01*STEP_LENGTH) ;
			sim.OdometryNoise()[2]= 0.0174532925 * 0.0174532925 ;

			int length = 100; int width = 20;

			for(int i=0 ; i<10;i++)
			{
				plan.addMotionAhead(length);
				plan.addMotionRotateLeft();
				plan.addMotionAhead(width);
				plan.addMotionRotateLeft();
				plan.addMotionAhead(length);
				plan.addMotionRotateRight();
				plan.addMotionAhead(width);
				plan.addMotionRotateRight();

			}


			sim.run();

			optimizer.addVertex(sim.RobotPoses()[0].poseSimulated);
			for(unsigned int i=1; i<sim.RobotPoses().size(); i++)
			{
				optimizer.addVertex(sim.RobotPoses()[i].poseSimulated);
				optimizer.addEdgeOdometry(i-1,i,sim.RobotPoses()[i].covariance.inverse());

			}
			// Wrong loop closures

			for(int i=0; i<19; i++)
			{
				int startNode = rand()%sim.RobotPoses().size();
				int endNode = rand()%sim.RobotPoses().size()/2;
				if(startNode - endNode < 20 ) i--;
				else
				{
					optimizer.addEdgeLoop(startNode,endNode,Eigen::Vector3d(0.,0.,0.),sim.RobotPoses()[2].covariance.inverse());
				}
			}

			// Correct Loop Clousers //

			for(int i=2; i<19; i++)
			{
					optimizer.addEdgeLoop(i*(length+1)+(i-1)*(width+1)-1,(length+width+2)*(i-2),Eigen::Vector3d(0,(i%2)?-20:20.0,M_PI),sim.RobotPoses()[2].covariance.inverse());
					optimizer.addEdgeLoop(50+(i-1)*122,50+(i-2)*122,Eigen::Vector3d(0,(i%2)?-20:20.0,M_PI),sim.RobotPoses()[2].covariance.inverse());
			}

			//optimizer.addEdgeLoop(sim.RobotPoses().size()-1,0,Eigen::Vector3d(0.,0.,0.),Eigen::Matrix3d::Identity());
			//optimizer.addEdgeLoop(sim.RobotPoses().size()-1,1,Eigen::Vector3d(0.,0.,0.),Eigen::Matrix3d::Identity());
			//optimizer.addEdgeLoop(220,0,Eigen::Vector3d(0.,0.,0.),Eigen::Matrix3d::Identity());
			//optimizer.addEdgeLoop(50,0,Eigen::Vector3d(0.,0.,0.),Eigen::Matrix3d::Identity());
			//optimizer.addEdgeLoop(100,0,Eigen::Vector3d(0.,0.,0.),Eigen::Matrix3d::Identity());

}

#endif
int main(int argc, char** argv)
{
	char *filename = NULL;
	bool hasInput = false;
	int numIterations = 3, clusterThreshold = 50;
	if(argc == 4)
	{
		filename = argv[1];
		numIterations = atoi(argv[2]);
		clusterThreshold = atoi(argv[3]);
		std::cout<<"Will read from "<<filename<<std::endl;
		hasInput = true;
	}
	ClusterOptimizer optimizer;
	
	if(!hasInput)
	{
		//lawn_mower_sim(optimizer);
	}
	else
	{
		if(!optimizer.read(filename))
		{
			std::cout<<"Could not read "<<filename<<std::endl;
			return 1;
		}
	}
	
#if 0
	optimizer.run(numIterations, clusterThreshold);

#else

//	6973
//	15331
//	21809
//	30144
//	optimizer.run(numIterations, clusterThreshold, 1761);
//	optimizer.run(numIterations, clusterThreshold, 8341);
//	optimizer.run(numIterations, clusterThreshold, 14754);
//	optimizer.run(numIterations, clusterThreshold, 22942);
//	optimizer.run(numIterations, clusterThreshold, 30896);

//		optimizer.run(numIterations, clusterThreshold, 6973);
//		optimizer.run(numIterations, clusterThreshold, 15331);
//		optimizer.run(numIterations, clusterThreshold, 21809);
//		optimizer.run(numIterations, clusterThreshold, 30144);
//		optimizer.run(numIterations, clusterThreshold, 40928);
// Laser Odometry Bicocca
	system("date >> icremental.txt");
	optimizer.run(numIterations, clusterThreshold, 7350);
	system("date >> icremental.txt");
	optimizer.run(numIterations, clusterThreshold, 16126);
	system("date >> icremental.txt");
	optimizer.run(numIterations, clusterThreshold, 22987);
	system("date >> icremental.txt");
	optimizer.run(numIterations, clusterThreshold, 31777);
	system("date >> icremental.txt");
	optimizer.run(numIterations, clusterThreshold, 43116);
	system("date >> icremental.txt");



#endif


	//std::vector<double> res;
	//std::cout<<optimizer.optimize(30)<<std::endl;
	//optimizer.refineOnResidualErrors(res);

	//std::cout<<optimizer.optimize(30)<<std::endl;
	//optimizer.refineOnResidualErrors(res);
	
	
	return 0;
	
}
