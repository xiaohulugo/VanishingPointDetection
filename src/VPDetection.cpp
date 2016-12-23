#include "VPDetection.h"
#include "time.h"
#include <iostream>

using namespace std;
using namespace cv;


VPDetection::VPDetection(void)
{
}


VPDetection::~VPDetection(void)
{
}

void VPDetection::run( std::vector<std::vector<double> > &lines, cv::Point2d pp, double f, std::vector<cv::Point3d> &vps, std::vector<std::vector<int> > &clusters )
{
	this->lines = lines;
	this->pp = pp;
	this->f = f;
	this->noiseRatio = 0.5; 

	cout<<"get vp hypotheses . . ."<<endl;
	std::vector<std::vector<cv::Point3d> > vpHypo;
	getVPHypVia2Lines( vpHypo );

	cout<<"get sphere grid . . ."<<endl;
	std::vector<std::vector<double> > sphereGrid;
	getSphereGrids( sphereGrid );

	cout<<"test vp hypotheses . . ."<<endl;
	getBestVpsHyp( sphereGrid, vpHypo, vps );

	cout<<"get final line clusters . . ."<<endl;
	double thAngle = 6.0 / 180.0 * CV_PI;
	lines2Vps( thAngle, vps, clusters );
	int clusteredNum = 0;
	for ( int i=0; i<3; ++i )
	{
		clusteredNum += clusters[i].size();
	}

	cout<<"total: " <<lines.size()<<"  clusered: "<<clusteredNum;
	cout<<"   X: "<<clusters[0].size()<<"   Y: "<<clusters[1].size()<<"   Z: "<<clusters[2].size()<<endl;
}

void VPDetection::getVPHypVia2Lines( std::vector<std::vector<cv::Point3d> > &vpHypo )
{
	int num = lines.size();

	double p = 1.0 / 3.0 * pow( 1.0 - noiseRatio, 2 );

	double confEfficience = 0.9999;
	int it = log( 1 - confEfficience ) / log( 1.0 - p );
	
	int numVp2 = 360;
	double stepVp2 = 2.0 * CV_PI / numVp2;

	// get the parameters of each line
	lineInfos.resize( num );
	for ( int i=0; i<num; ++i )
	{
		cv::Mat_<double> p1 = ( cv::Mat_<double>(3, 1) << lines[i][0], lines[i][1], 1.0 );
		cv::Mat_<double> p2 = ( cv::Mat_<double>(3, 1) << lines[i][2], lines[i][3], 1.0 );

		lineInfos[i].para = p1.cross( p2 );

		double dx = lines[i][0] - lines[i][2];
		double dy = lines[i][1] - lines[i][3];
		lineInfos[i].length = sqrt( dx * dx + dy * dy );

		lineInfos[i].orientation = atan2( dy, dx );
		if ( lineInfos[i].orientation < 0 )
		{
			lineInfos[i].orientation += CV_PI;
		}
	}

	// get vp hypothesis for each iteration
	vpHypo = std::vector<std::vector<cv::Point3d> > ( it * numVp2, 3 );
	int count = 0;
	srand((unsigned)time(NULL));  
	for ( int i = 0; i < it; ++ i )
	{
		int idx1 = rand() % num;
		int idx2 = rand() % num;
		while ( idx2 == idx1 )
		{
			idx2 = rand() % num;
		}

		// get the vp1
		cv::Mat_<double> vp1_Img = lineInfos[idx1].para.cross( lineInfos[idx2].para );
		if ( vp1_Img(2) == 0 )
		{
			i --;
			continue;
		}
		cv::Mat_<double> vp1 = ( cv::Mat_<double>(3, 1) << vp1_Img(0) / vp1_Img(2) - pp.x, vp1_Img(1) / vp1_Img(2) - pp.y, f );
		if ( vp1(2) == 0 ) { vp1(2) = 0.0011; }
		double N = sqrt( vp1(0) * vp1(0) + vp1(1) * vp1(1) + vp1(2) * vp1(2) );
		vp1 *= 1.0 / N;

		// get the vp2 and vp3
		cv::Mat_<double> vp2 = ( cv::Mat_<double>(3, 1) << 0.0, 0.0, 0.0 );
		cv::Mat_<double> vp3 = ( cv::Mat_<double>(3, 1) << 0.0, 0.0, 0.0 );
		for ( int j = 0; j < numVp2; ++ j )
		{
			// vp2
			double lambda = j * stepVp2;

			double k1 = vp1(0) * sin( lambda ) + vp1(1) * cos( lambda );
			double k2 = vp1(2);
			double phi = atan( - k2 / k1 );

			double Z = cos( phi );
			double X = sin( phi ) * sin( lambda );
			double Y = sin( phi ) * cos( lambda );

			vp2(0) = X;  vp2(1) = Y;  vp2(2) = Z;
			if ( vp2(2) == 0.0 ) { vp2(2) = 0.0011; }
			N = sqrt( vp2(0) * vp2(0) + vp2(1) * vp2(1) + vp2(2) * vp2(2) );
			vp2 *= 1.0 / N;
			if ( vp2(2) < 0 ) { vp2 *= -1.0; }		

			// vp3
			vp3 = vp1.cross( vp2 );
			if ( vp3(2) == 0.0 ) { vp3(2) = 0.0011; }
			N = sqrt( vp3(0) * vp3(0) + vp3(1) * vp3(1) + vp3(2) * vp3(2) );
			vp3 *= 1.0 / N;
			if ( vp3(2) < 0 ) { vp3 *= -1.0; }		

			//
			vpHypo[count][0] = cv::Point3d( vp1(0), vp1(1), vp1(2) );
			vpHypo[count][1] = cv::Point3d( vp2(0), vp2(1), vp2(2) );
			vpHypo[count][2] = cv::Point3d( vp3(0), vp3(1), vp3(2) );

			count ++;
		}
	}
}


void VPDetection::getSphereGrids( std::vector<std::vector<double> > &sphereGrid )
{	
	// build sphere grid with 1 degree accuracy
	double angelAccuracy = 1.0 / 180.0 * CV_PI;
	double angleSpanLA = CV_PI / 2.0;
	double angleSpanLO = CV_PI * 2.0;
	int gridLA = angleSpanLA / angelAccuracy;
	int gridLO = angleSpanLO / angelAccuracy;

	sphereGrid = std::vector<std::vector<double> >( gridLA, gridLO );
	for ( int i=0; i<gridLA; ++i )
	{
		for ( int j=0; j<gridLO; ++j )
		{
			sphereGrid[i][j] = 0.0;
		}
	}

	// put intersection points into the grid
	double angelTolerance = 60.0 / 180.0 * CV_PI;
	cv::Mat_<double> ptIntersect;
	double x = 0.0, y = 0.0;
	double X = 0.0, Y = 0.0, Z = 0.0, N = 0.0;
	double latitude = 0.0, longitude = 0.0;
	int LA = 0, LO = 0;
	double angleDev = 0.0;
	for ( int i=0; i<lines.size()-1; ++i )
	{
		for ( int j=i+1; j<lines.size(); ++j )
		{
			ptIntersect = lineInfos[i].para.cross( lineInfos[j].para );

			if ( ptIntersect(2,0) == 0 )
			{
				continue;
			}

			x = ptIntersect(0,0) / ptIntersect(2,0);
			y = ptIntersect(1,0) / ptIntersect(2,0);

			X = x - pp.x;
			Y = y - pp.y;
			Z = f;
			N = sqrt( X * X + Y * Y + Z * Z );

			latitude = acos( Z / N );
			longitude = atan2( X, Y ) + CV_PI;

			LA = int( latitude / angelAccuracy );
			if ( LA >= gridLA ) 
			{
				LA = gridLA - 1;
			}

			LO = int( longitude / angelAccuracy );
			if ( LO >= gridLO ) 
			{
				LO = gridLO - 1;
			}

			// 
			angleDev = abs( lineInfos[i].orientation - lineInfos[j].orientation );
			angleDev = min( CV_PI - angleDev, angleDev );
			if ( angleDev > angelTolerance )
			{
				continue;
			}

			sphereGrid[LA][LO] += sqrt( lineInfos[i].length * lineInfos[j].length ) * ( sin( 2.0 * angleDev ) + 0.05 );
		}
	}

	// 
	int halfSize = 1;
	int winSize = halfSize * 2 + 1;
	int neighNum = winSize * winSize;

	// get the weighted line length of each grid
	std::vector<std::vector<double> > sphereGridNew( gridLA, gridLO );
	for ( int i=halfSize; i<gridLA-halfSize; ++i )
	{
		for ( int j=halfSize; j<gridLO-halfSize; ++j )
		{
			double neighborTotal = 0.0;
			for ( int m=0; m<winSize; ++m )
			{
				for ( int n=0; n<winSize; ++n )
				{
					neighborTotal += sphereGrid[i-halfSize+m][j-halfSize+n];
				}
			}

			sphereGridNew[i][j] = sphereGrid[i][j] + neighborTotal / neighNum;
		}
	}
	sphereGrid = sphereGridNew;
}

void VPDetection::getBestVpsHyp( std::vector<std::vector<double> > &sphereGrid, std::vector<std::vector<cv::Point3d> > &vpHypo, std::vector<cv::Point3d> &vps )
{
	int num = vpHypo.size();
	double oneDegree = 1.0 / 180.0 * CV_PI;

	// get the corresponding line length of every hypotheses
	std::vector<double> lineLength( num, 0.0 );
	for ( int i = 0; i < num; ++ i )
	{
		std::vector<cv::Point2d> vpLALO( 3 ); 
		for ( int j = 0; j < 3; ++ j )
		{
			if ( vpHypo[i][j].z == 0.0 )
			{
				continue;
			}

			if ( vpHypo[i][j].z > 1.0 || vpHypo[i][j].z < -1.0 )
			{
				cout<<1.0000<<endl;
			}
			double latitude = acos( vpHypo[i][j].z );
			double longitude = atan2( vpHypo[i][j].x, vpHypo[i][j].y ) + CV_PI;

			int gridLA = int( latitude / oneDegree );
			if ( gridLA == 90 ) 
			{
				gridLA = 89;
			}
			
			int gridLO = int( longitude / oneDegree );
			if ( gridLO == 360 ) 
			{
				gridLO = 359;
			}

			lineLength[i] += sphereGrid[gridLA][gridLO];
		}
	}

	// get the best hypotheses
	int bestIdx = 0;
	double maxLength = 0.0;
	for ( int i = 0; i < num; ++ i )
	{
		if ( lineLength[i] > maxLength )
		{
			maxLength = lineLength[i];
			bestIdx = i;
		}
	}

	vps = vpHypo[bestIdx];
}


void VPDetection::lines2Vps( double thAngle, std::vector<cv::Point3d> &vps, std::vector<std::vector<int> > &clusters )
{
	clusters.clear();
	clusters.resize( 3 );

	//get the corresponding vanish points on the image plane
	std::vector<cv::Point2d> vp2D( 3 ); 
	for ( int i = 0; i < 3; ++ i )
	{
		vp2D[i].x =  vps[i].x * f / vps[i].z + pp.x;
		vp2D[i].y =  vps[i].y * f / vps[i].z + pp.y;
	}

	for ( int i = 0; i < lines.size(); ++ i )
	{
		double x1 = lines[i][0];
		double y1 = lines[i][1];
		double x2 = lines[i][2];
		double y2 = lines[i][3];
		double xm = ( x1 + x2 ) / 2.0;
		double ym = ( y1 + y2 ) / 2.0;

		double v1x = x1 - x2;
		double v1y = y1 - y2;
		double N1 = sqrt( v1x * v1x + v1y * v1y );
		v1x /= N1;   v1y /= N1;

		double minAngle = 1000.0;
		int bestIdx = 0;
		for ( int j = 0; j < 3; ++ j )
		{
			double v2x = vp2D[j].x - xm;
			double v2y = vp2D[j].y - ym;
			double N2 = sqrt( v2x * v2x + v2y * v2y );
			v2x /= N2;  v2y /= N2;

			double crossValue = v1x * v2x + v1y * v2y;
			if ( crossValue > 1.0 )
			{
				crossValue = 1.0;
			}
			if ( crossValue < -1.0 )
			{
				crossValue = -1.0;
			}
			double angle = acos( crossValue );
			angle = min( CV_PI - angle, angle );

			if ( angle < minAngle )
			{
				minAngle = angle;
				bestIdx = j;
			}
		}

		//
		if ( minAngle < thAngle )
		{
			clusters[bestIdx].push_back( i );
		}
	}
}